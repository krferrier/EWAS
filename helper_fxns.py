from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Iterable, Optional, Union
import os
import pandas as pd

def _to_bool(x: Union[str, bool]) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, str):
        return x.strip().lower() in {"y", "yes", "true", "1"}
    return bool(x)

def _norm_path(p: Union[str, Path]) -> Path:
    # Expand env vars and ~, then resolve relative paths against CWD
    if isinstance(p, Path):
        s = str(p)
    else:
        s = p
    s = os.path.expandvars(os.path.expanduser(s))
    return Path(s)

class ConfigWizard(object):
    __slots__ = (
        "config", "mvals", "pheno", "assoc_var", "stratified", "strat_vars",
        "dmr", "genome_build", "min_pval", "win_size", "region_filter",
        "chunk_size", "processing_type", "n_workers", "out_dir", "out_type",
        "_groups", "_bacon_plot_kinds", "dmr_anno_cache_dir", "dmr_anno_cache_tag",
        "gene_table", "ucsc_database_base", "ucsc_gbdb_base", "hgnc_complete_set_url"
    )

    def __init__(self, cfg: Dict):
        self.config = cfg

        # --- Required keys (raise early with a good message if missing) ---
        required = [
            "mvals", "pheno", "association_variable", "stratified_ewas",
            "stratify_variables", "dmr_analysis", "genome_build", "min_pvalue",
            "window_size", "region_filter", "chunk_size", "processing_type",
            "workers", "out_directory", "out_type",
        ]
        missing = [k for k in required if k not in cfg]
        if missing:
            raise KeyError(f"Missing required config keys: {', '.join(missing)}")

        # --- Basic fields ---
        self.mvals = _norm_path(cfg["mvals"])
        self.pheno = _norm_path(cfg["pheno"])

        self.assoc_var: str = str(cfg["association_variable"])
        self.stratified: bool = _to_bool(cfg["stratified_ewas"])
        self.strat_vars: List[str] = list(cfg.get("stratify_variables") or [])

        self.dmr: bool = _to_bool(cfg["dmr_analysis"])
        self.genome_build: str = str(cfg["genome_build"])
        self.min_pval: float = float(cfg["min_pvalue"])
        self.win_size: int = int(cfg["window_size"])
        self.region_filter: float = float(cfg["region_filter"])

        self.chunk_size: int = int(cfg["chunk_size"])
        self.processing_type: str = str(cfg["processing_type"])
        self.n_workers: int = int(cfg["workers"])

        # Ensure the output directory is a Path, and does NOT need a trailing slash
        self.out_dir: Path = _norm_path(cfg["out_directory"])
        self.out_type: str = str(cfg["out_type"])  # e.g. ".csv" or ".csv.gz"

        # Optional DMR annotation-cache settings
        anno_cfg = cfg.get("dmr_annotation", {}) or {}

        self.dmr_anno_cache_dir: Path = _norm_path(
            anno_cfg.get("cache_dir", cfg.get("dmr_anno_cache_dir", "resources/ucsc"))
        )
        self.dmr_anno_cache_tag: str = str(
            anno_cfg.get("cache_tag", cfg.get("dmr_anno_cache_tag", "latest"))
        )
        self.gene_table: str = str(
            anno_cfg.get("gene_table", cfg.get("gene_table", "refGene"))
        )
        self.ucsc_database_base: str = str(
            anno_cfg.get(
                "ucsc_database_base",
                cfg.get("ucsc_database_base", "https://hgdownload.soe.ucsc.edu/goldenPath"),
            )
        ).rstrip("/")
        self.ucsc_gbdb_base: str = str(
            anno_cfg.get(
                "ucsc_gbdb_base",
                cfg.get("ucsc_gbdb_base", "https://hgdownload.soe.ucsc.edu/gbdb"),
            )
        ).rstrip("/")
        self.hgnc_complete_set_url: str = str(
            anno_cfg.get(
                "hgnc_complete_set_url",
                cfg.get(
                    "hgnc_complete_set_url",
                    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                ),
            )
        )


        # Keep plot kinds centralized
        self._bacon_plot_kinds: List[str] = ["traces", "posteriors", "fit", "qqs"]

        # Defer computing groups until requested (but you can force with CW.groups)
        self._groups: Optional[List[str]] = None

    # ---------- Commonly-used simple properties ----------
    @property
    def bacon_plot_kinds(self) -> List[str]:
        return list(self._bacon_plot_kinds)

    @property
    def groups(self) -> List[str]:
        """Observed strat groups from phenotype file. Returns ["all"] if not stratified."""
        if self._groups is not None:
            return self._groups

        if not self.stratified or not self.strat_vars:
            self._groups = ["all"]
            return self._groups

        # Load phenotype and build observed combos
        df = pd.read_csv(self.pheno)
        # Convert everything used for strat to string to avoid 0 vs "0" issues
        for col in self.strat_vars:
            if col not in df.columns:
                raise KeyError(f"Stratify variable '{col}' not found in phenotype file.")
        df[self.strat_vars] = df[self.strat_vars].astype(str)

        combos = (
            df.groupby(self.strat_vars)
              .size()
              .reset_index()
        )
        combos["combination"] = combos[self.strat_vars].agg("_".join, axis=1)
        observed = combos["combination"].tolist()
        # Guarantee deterministic sort
        self._groups = sorted(set(observed))
        return self._groups

    # ---------- Path helpers ----------
    def _prefix(self, group: Optional[str] = None) -> str:
        """
        "BMI"                       (unstratified)
        "female_BMI" or "F_1_BMI"  (stratified with observed group name)
        """
        if not group or group == "all":
            return f"{self.assoc_var}"
        return f"{group}_{self.assoc_var}"

    def _out(self, *parts: Union[str, Path]) -> Path:
        return self.out_dir.joinpath(*map(lambda p: str(p), parts))

    # ---------- Unstratified outputs ----------
    @property
    def raw_results(self) -> Path:
        return self._out(f"{self._prefix()}_ewas_results{self.out_type}")

    @property
    def bacon_results(self) -> Path:
        return self._out(f"{self._prefix()}_ewas_bacon_results{self.out_type}")

    @property
    def manhattan_qq_plot(self) -> Path:
        # Keep jpg since that’s what your Snakefile showed
        return self._out(f"{self._prefix()}_ewas_manhattan_qq_plots.jpg")

    @property
    def annotated_results(self) -> Path:
        return self._out(f"{self._prefix()}_ewas_annotated_results{self.out_type}")

    @property
    def meta_analysis_results(self) -> Path:
        return self._out(f"{self._prefix()}_ewas_meta_analysis_results_1.txt")

    def bacon_plot_files(self) -> List[str]:
        # Return strings for Snakemake expand friendliness
        return [str(self._out("bacon_plots", f"{self._prefix()}_{k}.jpg"))
                for k in self._bacon_plot_kinds]

    # ---------- Stratified outputs ----------
    def strat_raw_results(self) -> List[str]:
        if self.groups == ["all"]:
            return []
        return [str(self._out(g, f"{self._prefix(g)}_ewas_results{self.out_type}"))
                for g in self.groups]

    def strat_bacon_results(self) -> List[str]:
        if self.groups == ["all"]:
            return []
        return [str(self._out(g, f"{self._prefix(g)}_ewas_bacon_results{self.out_type}"))
                for g in self.groups]

    def strat_bacon_plot_files(self) -> List[str]:
        if self.groups == ["all"]:
            return []
        files = []
        for g in self.groups:
            files.extend(
                str(self._out(g, "bacon_plots", f"{self._prefix(g)}_{k}.jpg"))
                for k in self._bacon_plot_kinds
            )
        return files

    # ---------- DMR outputs ----------
    @property
    def dmr_results_bed(self) -> Path:
        return self._out(f"{self._prefix()}_ewas_annotated_results.bed")

    @property
    def dmr_acf(self) -> Path:
        return self._out("dmr", f"{self._prefix()}_ewas.acf.txt")

    @property
    def dmr_args(self) -> Path:
        return self._out("dmr", f"{self._prefix()}_ewas.args.txt")

    @property
    def dmr_fdr(self) -> Path:
        # Keep gz suffix per your Snakefile
        return self._out("dmr", f"{self._prefix()}_ewas.fdr.bed.gz")

    @property
    def dmr_regions(self) -> Path:
        return self._out("dmr", f"{self._prefix()}_ewas.regions-p.bed.gz")

    @property
    def dmr_slk(self) -> Path:
        return self._out("dmr", f"{self._prefix()}_ewas.slk.bed.gz")

    @property
    def dmr_anno(self) -> Path:
        return self._out("dmr", f"{self._prefix()}_ewas.anno.{self.genome_build}.bed")

    @property
    def dmr_cpg_anno(self):
        return self._out("dmr", f"{self._prefix()}_ewas.anno.{self.genome_build}.with_cpgs.bed")

    # ---------- Local DMR annotation-cache resources ----------
    @property
    def dmr_anno_resource_dir(self) -> Path:
        return self.dmr_anno_cache_dir.joinpath(self.genome_build, self.dmr_anno_cache_tag)

    @property
    def dmr_refgene_txt(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath(f"{self.gene_table}.txt.gz")

    @property
    def dmr_cpg_island_txt(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("cpgIslandExt.txt.gz")

    @property
    def dmr_hgnc_bb(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("hgnc.bb")

    @property
    def dmr_hgnc_complete_set(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("hgnc_complete_set.txt")

    @property
    def dmr_refgene_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath(f"{self.gene_table}.bed.gz")

    @property
    def dmr_cpg_island_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("cpgIslandExt.bed.gz")

    @property
    def dmr_hgnc_raw_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("hgnc.raw.bed.gz")

    @property
    def dmr_hgnc_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("hgnc.approved.bed.gz")

    @property
    def dmr_annotation_manifest(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("annotation_manifest.tsv")