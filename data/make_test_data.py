#!/usr/bin/env python3
"""Generate the workflow's test inputs: data/pheno.csv and data/mvals.csv.gz.

The data are simulated but laid out on real EPIC probes, so every stage of the
workflow has something meaningful to do:

* EWAS / BACON: 120 samples (60 F, 60 M) -- enough per stratum for a stable
  fit -- and a null background, so lambda sits near 1 and the QQ plot has a
  real null line to depart from.
* Known answers: a set of CpGs carries a planted BMI effect. Every planted CpG
  and its true effect is written to data/test_truth.tsv, so a run can be
  checked for what it recovers, not just for whether it finishes.
* DMR (comb-p): regions of 6-15 EPIC probes spaced within the comb-p window,
  half with a coordinated BMI effect and half as correlated null regions (the
  negative control comb-p should not call).
* Enrichment: most single-CpG hits are drawn from ChromHMM active-TSS (TssA)
  probes, so KYCG enrichment has a real, known signal to find.
* Plots: CpGs on every chromosome including X and Y, with realistic X
  inactivation in females; effects are moderate (roughly p 1e-4 to 1e-20),
  so nothing underflows and the figures look like a real study.

Row order is identical in the two files. The combined (unstratified) EWAS pairs
methylation rows with phenotype rows by position, so this matters.

Usage (from the repo root, after one workflow run has filled the annotation
cache -- or point --manifest at any copy of the Zhou gene manifest):

    python data/make_test_data.py

Needs numpy and pandas. The ChromHMM bias needs `yame` on PATH and the KYCG
cache; without them hits are drawn uniformly and a note is printed.
"""
import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path("resources/annotation/zhou/EPIC/v8.1")
STD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", default=str(CACHE / "EPIC.hg38.manifest.gencode.v41.tsv.gz"),
                    help="Zhou EPIC gene manifest (CpG_chrm, CpG_beg, probeID columns)")
    ap.add_argument("--kycg-dir", default=str(CACHE / "KYCG"),
                    help="KYCG cache holding ChromHMM.cm (optional)")
    ap.add_argument("--ordering", default=str(CACHE / "EPIC.ordering.tsv.gz"),
                    help="EPIC probe ordering the .cm files are aligned to (optional)")
    ap.add_argument("--out-dir", default="data")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--n-per-sex", type=int, default=60)
    ap.add_argument("--n-background", type=int, default=4800,
                    help="randomly sampled CpGs genome-wide")
    ap.add_argument("--min-chrY", type=int, default=20,
                    help="top up chrY so the Y block is visible on the Manhattan plot")
    ap.add_argument("--n-dmr", type=int, default=10, help="regions with a planted effect")
    ap.add_argument("--n-null-regions", type=int, default=10, help="correlated regions, no effect")
    ap.add_argument("--n-hits", type=int, default=60, help="single CpGs with a planted effect")
    ap.add_argument("--tssa-share", type=float, default=2 / 3,
                    help="share of single-CpG hits drawn from TssA probes")
    # 0 by default: the README asks for cleaned input with no missing values.
    # Raise it only to test how the workflow copes with input that breaks that.
    ap.add_argument("--na-rate", type=float, default=0.0, help="missing M-values, at random")
    return ap.parse_args()


def chromhmm_states(kycg_dir, ordering):
    """Probe -> ChromHMM state via `yame unpack`, or None if unavailable."""
    cm = Path(kycg_dir) / "ChromHMM.cm"
    if not (cm.exists() and Path(ordering).exists() and shutil.which("yame")):
        print("note: ChromHMM.cm, the ordering file or yame not found; "
              "single-CpG hits will be drawn uniformly", file=sys.stderr)
        return None
    labels = subprocess.run(["yame", "unpack", "-a", str(cm)], check=True,
                            capture_output=True, text=True).stdout.splitlines()
    with gzip.open(ordering, "rt") as f:
        next(f)
        ids = [line.split("\t", 1)[0] for line in f]
    if len(ids) != len(labels):
        sys.exit(f"ordering has {len(ids)} rows but ChromHMM.cm has {len(labels)}")
    return pd.Series(labels, index=ids, name="chromhmm")


def dense_regions(probes, max_step=200, n_min=6, n_max=15):
    """Runs of autosomal probes each within max_step bp of the previous one."""
    q = probes[probes.CpG_chrm.isin(STD_CHROMS[:22])].sort_values(["CpG_chrm", "CpG_beg"])
    gap = q.CpG_beg.diff().where(q.CpG_chrm.eq(q.CpG_chrm.shift()))
    cl = (gap.isna() | (gap > max_step)).cumsum()
    size = cl.map(cl.value_counts())
    q = q.assign(cluster=cl.values)[(size >= n_min) & (size <= n_max)]
    return q


def main():
    a = parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out_dir)

    # ---- probes -----------------------------------------------------------
    man = pd.read_csv(a.manifest, sep="\t", usecols=["CpG_chrm", "CpG_beg", "probeID"],
                      low_memory=False)
    man = man[man.probeID.str.startswith("cg") & man.CpG_chrm.notna()]
    man = man.drop_duplicates("probeID")
    st = chromhmm_states(a.kycg_dir, a.ordering)
    man["chromhmm"] = man.probeID.map(st) if st is not None else "NA"
    man["chromhmm"] = man.chromhmm.fillna("NA")

    # Regions first, so background sampling cannot split one.
    reg = dense_regions(man)
    ids = rng.choice(reg.cluster.unique(), a.n_dmr + a.n_null_regions, replace=False)
    reg = reg[reg.cluster.isin(ids)].copy()
    rank = {c: i for i, c in enumerate(ids)}
    reg["region"] = reg.cluster.map(rank)
    reg["role"] = np.where(reg.region < a.n_dmr, "dmr", "null_region")
    reg["region"] = [f"{'DMR' if r == 'dmr' else 'NULL'}{k + 1:02d}"
                     for r, k in zip(reg.role, reg.region)]

    rest = man[~man.probeID.isin(reg.probeID)]
    bg = rest.sample(a.n_background, random_state=a.seed)
    n_y = (bg.CpG_chrm == "chrY").sum()
    if n_y < a.min_chrY:
        extra = rest[(rest.CpG_chrm == "chrY") & ~rest.probeID.isin(bg.probeID)]
        bg = pd.concat([bg, extra.sample(min(len(extra), a.min_chrY - n_y), random_state=a.seed)])
    bg = bg.assign(role="background", region="")

    # Single-CpG hits, on standard chromosomes, mostly from active TSSs.
    pool = bg[bg.CpG_chrm.isin(STD_CHROMS[:23])]
    n_tssa = int(round(a.n_hits * a.tssa_share)) if st is not None else 0
    tssa = pool[pool.chromhmm == "TssA"].sample(n_tssa, random_state=a.seed)
    other = pool.drop(tssa.index).sample(a.n_hits - n_tssa, random_state=a.seed + 1)
    bg.loc[tssa.index.union(other.index), "role"] = "single_hit"

    cpg = pd.concat([bg, reg.drop(columns="cluster")]).reset_index(drop=True)
    p = len(cpg)

    # ---- samples ----------------------------------------------------------
    n = 2 * a.n_per_sex
    sex = np.array(["F"] * a.n_per_sex + ["M"] * a.n_per_sex)
    order = rng.permutation(n)
    sex = sex[order]
    pheno = pd.DataFrame({
        "sampleID": [f"S{i + 1:03d}" for i in range(n)],
        "sex": sex,
        "re": rng.choice([1, 2, 3], size=n, p=[0.5, 0.3, 0.2]),
        "BMI": np.clip(rng.normal(27, 5, n), 17, 45).round(1),
    })
    female = (pheno.sex == "F").to_numpy()
    bmi_c = (pheno.BMI - pheno.BMI.mean()).to_numpy()

    # ---- methylation, in M-values ------------------------------------------
    # Baseline by chromatin state: active promoters unmethylated, repressed
    # and bivalent intermediate, everything else largely methylated.
    state = cpg.chromhmm.to_numpy()
    low = np.isin(state, ["TssA", "TssFlnk", "TssFlnkU", "TssFlnkD"])
    mid = np.isin(state, ["TssBiv", "EnhBiv", "ReprPC", "ReprPCWk"])
    beta0 = np.where(low, rng.beta(1.2, 12, p),
                     np.where(mid, rng.beta(3, 3, p), rng.beta(10, 2, p)))
    beta0 = beta0.clip(0.01, 0.99)
    m0 = np.log2(beta0 / (1 - beta0))
    sigma = rng.uniform(0.2, 0.5, p)
    M = m0[None, :] + rng.normal(0, 1, (n, p)) * sigma[None, :]

    # Small per-CpG differences by `re`.
    re_eff = rng.normal(0, 0.05, (3, p))
    M += re_eff[pheno.re.to_numpy() - 1]

    # Sex chromosomes. X inactivation pulls unmethylated X promoters in
    # females towards 50%; chrY in females is background noise.
    onX = (cpg.CpG_chrm == "chrX").to_numpy()
    onY = (cpg.CpG_chrm == "chrY").to_numpy()
    xi = onX & low
    M[np.ix_(female, xi)] += 0.8 * (0.0 - m0[xi])[None, :]
    M[np.ix_(female, onY)] = rng.normal(-3.5, 1.0, (female.sum(), onY.sum()))

    # Planted BMI effects, in M-value units per kg/m^2.
    beta_bmi = np.zeros(p)
    hit = (cpg.role == "single_hit").to_numpy()
    beta_bmi[hit] = rng.uniform(0.02, 0.06, hit.sum()) * rng.choice([-1, 1], hit.sum())

    # Regions: a shared per-sample factor makes neighbouring CpGs correlated,
    # as they are in real data; DMRs add a coordinated, same-sign effect.
    for name, rows in cpg[cpg.role.isin(["dmr", "null_region"])].groupby("region").groups.items():
        rows = np.asarray(rows)
        f = rng.normal(0, 0.3, n)
        M[:, rows] += f[:, None] * rng.uniform(0.5, 1.0, len(rows))[None, :]
        if name.startswith("DMR"):
            b = rng.uniform(0.025, 0.05) * rng.choice([-1, 1])
            beta_bmi[rows] = b * rng.uniform(0.6, 1.2, len(rows))
    M += bmi_c[:, None] * beta_bmi[None, :]

    if a.na_rate > 0:
        M[rng.random((n, p)) < a.na_rate] = np.nan

    # ---- write ------------------------------------------------------------
    out.mkdir(parents=True, exist_ok=True)
    pheno.to_csv(out / "pheno.csv", index=False)
    mv = pd.DataFrame(M.round(3), columns=cpg.probeID)
    mv.insert(0, "sampleID", pheno.sampleID)
    mv.to_csv(out / "mvals.csv.gz", index=False, compression="gzip")

    truth = cpg.assign(true_beta_per_BMI=beta_bmi.round(4), baseline_beta=beta0.round(3))
    truth = truth[truth.role != "background"].rename(
        columns={"CpG_chrm": "chrom", "CpG_beg": "pos"})
    truth = truth[["probeID", "chrom", "pos", "role", "region", "chromhmm",
                   "true_beta_per_BMI", "baseline_beta"]].sort_values(["role", "chrom", "pos"])
    truth.to_csv(out / "test_truth.tsv", sep="\t", index=False)

    print(f"samples {n} ({female.sum()} F) | CpGs {p}: "
          f"{(cpg.role == 'background').sum()} background, {hit.sum()} single hits "
          f"({n_tssa} TssA), {a.n_dmr} DMRs + {a.n_null_regions} null regions "
          f"({cpg.role.isin(['dmr', 'null_region']).sum()} CpGs) | chrX {onX.sum()}, "
          f"chrY {onY.sum()}, other {(~cpg.CpG_chrm.isin(STD_CHROMS)).sum()} | "
          f"NA {np.isnan(M).mean():.2%}")


if __name__ == "__main__":
    main()
