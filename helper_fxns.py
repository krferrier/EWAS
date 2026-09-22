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
        "_groups", "_bacon_plot_kinds", "anno_cache_dir", "dmr_anno_cache_tag",
        "gene_table", "ucsc_database_base", "ucsc_gbdb_base",
        "array_platform", "zhou_release", "gencode_release", "zhou_raw_base",
        "zhou_anno_raw_base", "zhou_anno_api_base",
        "bios_eqtm_url", "hgnc_complete_set_url",
        "enrichment", "enrich_significance", "enrich_threshold",
        "enrich_min_set_size", "ewas_atlas_url", "_enrich_kbs",
        "_enrich_qc_kbs",
        "enrichment_plots", "enrich_plot_top_n",
        "dmr_make_zoom", "dmr_make_combined", "dmr_plot_min_probes",
        "dmr_plot_max_y", "dmr_zoom_padding", "dmr_cluster_gap",
        "dmr_combined_formats",
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

        # --- Annotation settings (one block, one cache) ---
        # `annotation` covers every annotation resource the workflow downloads:
        # the Zhou Infinium manifest used for CpG-level annotation, the UCSC
        # tracks used for CpG islands and for DMR gene annotation, and the
        # eQTM/HGNC tables. They share a single cache root so there is one
        # place to inspect, archive, or delete.
        anno_cfg = cfg.get("annotation", {}) or {}

        self.anno_cache_dir: Path = _norm_path(
            anno_cfg.get("cache_dir", "resources/annotation")
        )
        self.dmr_anno_cache_tag: str = str(anno_cfg.get("cache_tag", "latest"))
        self.gene_table: str = str(anno_cfg.get("gene_table", "refGene"))
        self.ucsc_database_base: str = str(
            anno_cfg.get(
                "ucsc_database_base", "https://hgdownload.soe.ucsc.edu/goldenPath"
            )
        ).rstrip("/")
        self.ucsc_gbdb_base: str = str(
            anno_cfg.get("ucsc_gbdb_base", "https://hgdownload.soe.ucsc.edu/gbdb")
        ).rstrip("/")

        # Zhou lab Infinium annotation moved off zhouserver.research.chop.edu to
        # GitHub. Gene/promoter annotation now comes from the "coherent tables"
        # repo zhou-lab/InfiniumAnnotationData, pinned by git tag.
        self.array_platform: str = str(anno_cfg.get("array_platform", "EPIC"))
        # Git tag in zhou-lab/InfiniumAnnotationData. "main" tracks the latest
        # release; a tag such as "v8.1" freezes it for a reproducible run.
        self.zhou_release: str = str(anno_cfg.get("zhou_release", "v8.1"))
        # GENCODE release embedded in the manifest filename. This is genome
        # specific: hg38 uses v41, hg19 uses v26lift37.
        self.gencode_release: str = str(anno_cfg.get("gencode_release", "v41"))
        self.zhou_raw_base: str = str(
            anno_cfg.get(
                "zhou_raw_base",
                "https://github.com/zhou-lab/InfiniumAnnotationData/raw",
            )
        ).rstrip("/")
        # The KYCG feature sets and the probe ordering file live in the OTHER
        # Zhou repo (InfiniumAnnotation, the versioned "coherent" release),
        # not in InfiniumAnnotationData.
        self.zhou_anno_raw_base: str = str(
            anno_cfg.get(
                "zhou_anno_raw_base",
                "https://github.com/zhou-lab/InfiniumAnnotation/raw",
            )
        ).rstrip("/")
        self.zhou_anno_api_base: str = str(
            anno_cfg.get(
                "zhou_anno_api_base",
                "https://api.github.com/repos/zhou-lab/InfiniumAnnotation/contents",
            )
        ).rstrip("/")
        self.bios_eqtm_url: str = str(
            anno_cfg.get(
                "bios_eqtm_url",
                "https://molgenis26.gcc.rug.nl/downloads/biosqtlbrowser/"
                "2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt",
            )
        )
        self.hgnc_complete_set_url: str = str(
            anno_cfg.get(
                "hgnc_complete_set_url",
                "https://storage.googleapis.com/public-download-files/hgnc/"
                "archive/archive/quarterly/tsv/hgnc_complete_set_2025-07-01.txt",
            )
        )

        # --- Functional enrichment of the significant CpG set ---
        enrich_cfg = cfg.get("enrichment", {}) or {}
        self.enrichment: bool = _to_bool(enrich_cfg.get("run", "yes"))
        # How "significant" is defined when building the foreground set.
        # The background is always the set of CpGs actually tested, never the
        # whole array -- using the array would inflate every enrichment by the
        # coverage bias of the probes that failed QC.
        self.enrich_significance: str = str(
            enrich_cfg.get("significance", "fdr")
        ).lower()
        self.enrich_threshold: float = float(enrich_cfg.get("threshold", 0.05))
        self.enrich_min_set_size: int = int(enrich_cfg.get("min_set_size", 20))
        # Which KYCG knowledgebases to test. A list of set names, or "all".
        self._enrich_kbs = enrich_cfg.get("knowledgebases")
        self._enrich_qc_kbs = enrich_cfg.get("qc_knowledgebases")
        self.ewas_atlas_url: str = str(
            enrich_cfg.get(
                "ewas_atlas_url",
                "https://ngdc.cncb.ac.cn/ewas/downloads/batch"
                "?file=EWAS_Atlas_associations.tsv",
            )
        )
        # One summary plot per enrichment analysis. Always written when
        # enabled, including when nothing reached significance -- the plot then
        # states that, matching how the enrichment tables write an empty table
        # naming the reason. That keeps the outputs trackable by Snakemake.
        self.enrichment_plots: bool = _to_bool(enrich_cfg.get("make_plots", "yes"))
        self.enrich_plot_top_n: int = int(enrich_cfg.get("plot_top_n", 10))

        # --- DMR plotting ---
        # scripts/dmr_plot.R exposes more knobs than are surfaced here; the
        # rest keep the script's own defaults. Only settings that change which
        # FILES appear, or that depend on the data, are configurable.
        dmr_plot_cfg = cfg.get("dmr_plots", {}) or {}
        # Zoom plots and the combined figure are off by default because the
        # number of files is unknown until comb-p has run -- one zoom plot and
        # one gene table per significant region cluster.
        self.dmr_make_zoom: str = "yes" if _to_bool(
            dmr_plot_cfg.get("make_zoom", "no")) else "no"
        self.dmr_make_combined: str = "yes" if _to_bool(
            dmr_plot_cfg.get("make_combined", "no")) else "no"
        self.dmr_plot_min_probes: int = int(dmr_plot_cfg.get("min_probes", 2))
        self.dmr_plot_max_y: float = float(dmr_plot_cfg.get("max_y", -1))
        self.dmr_zoom_padding: int = int(dmr_plot_cfg.get("zoom_padding", 2000))
        self.dmr_cluster_gap: int = int(dmr_plot_cfg.get("cluster_gap", 3000))
        self.dmr_combined_formats: str = str(
            dmr_plot_cfg.get("combined_formats", "pdf"))

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

    # ---------- Group-specific paths ----------
    def group_dir(self, group: str) -> Path:
        """Directory containing outputs for one stratum."""
        return self._out(group)

    def group_pheno(self, group: str) -> Path:
        return self.group_dir(group) / f"{group}_pheno.fst"

    def group_mvals(self, group: str) -> Path:
        return self.group_dir(group) / f"{group}_mvals.fst"

    def group_ewas_results(self, group: str) -> Path:
        return self.group_dir(group) / (
            f"{group}_{self.assoc_var}_ewas_results{self.out_type}"
        )

    def group_bacon_results(self, group: str) -> Path:
        return self.group_dir(group) / (
            f"{group}_{self.assoc_var}_ewas_bacon_results{self.out_type}"
        )

    def group_bacon_plot(self, group: str, kind: str) -> Path:
        if kind not in self._bacon_plot_kinds:
            raise ValueError(
                f"Unknown bacon plot kind '{kind}'. "
                f"Allowed values: {', '.join(self._bacon_plot_kinds)}"
            )

        return self.group_dir(group) / "bacon_plots" / (
            f"{group}_{self.assoc_var}_{kind}.jpg"
        )

    @property
    def stratified_pheno_files(self) -> List[str]:
        """Concrete phenotype outputs for every observed group."""
        if self.groups == ["all"]:
            return []

        return [str(self.group_pheno(group)) for group in self.groups]

    @property
    def stratified_mvals_files(self) -> List[str]:
        """Concrete methylation outputs for every observed group."""
        if self.groups == ["all"]:
            return []

        return [str(self.group_mvals(group)) for group in self.groups]

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

    # ---------- Stratified EWAS outputs ----------
    def strat_raw_results(self) -> List[str]:
        if self.groups == ["all"]:
            return []

        return [
            str(self.group_ewas_results(group))
            for group in self.groups
        ]

    def strat_bacon_results(self) -> List[str]:
        if self.groups == ["all"]:
            return []

        return [
            str(self.group_bacon_results(group))
            for group in self.groups
        ]

    def strat_bacon_plot_files(self) -> List[str]:
        if self.groups == ["all"]:
            return []
        return [
            str(self.group_bacon_plot(group, kind))
            for group in self.groups
            for kind in self._bacon_plot_kinds
        ]

    # ---------- METAL meta-analysis outputs ----------
    @property
    def metal_out_prefix(self) -> Path:
        """
        Prefix passed to METAL's OUTFILE command.

        METAL appends its output index and extension:
        <prefix>1.txt
        """
        return self._out(f"{self.assoc_var}_ewas_meta_analysis_results_")

    @property
    def meta_analysis_results(self) -> Path:
        """Expected first METAL meta-analysis output."""
        return Path(f"{self.metal_out_prefix}1.txt")

    @property
    def metal_command_script(self) -> Path:
        """Generated METAL command file; not a source script."""
        return self._out(
            "meta_analysis",
            f"{self.assoc_var}_metal_commands.sh",
        )

    # ---------- DMR outputs ----------
    @property
    def dmr_out_dir(self) -> Path:
        return self._out("dmr")

    @property
    def dmr_out_prefix(self) -> Path:
        return self._out("dmr", f"{self.assoc_var}_ewas")

    @property
    def dmr_results_bed(self) -> Path:
        """Combined EWAS BED file supplied to comb-p."""
        return self._out(f"{self.assoc_var}_ewas_annotated_results.bed")

    @property
    def dmr_acf(self) -> Path:
        return Path(f"{self.dmr_out_prefix}.acf.txt")

    @property
    def dmr_args(self) -> Path:
        return Path(f"{self.dmr_out_prefix}.args.txt")

    @property
    def dmr_fdr(self) -> Path:
        return Path(f"{self.dmr_out_prefix}.fdr.bed.gz")

    @property
    def dmr_regions(self) -> Path:
        return Path(f"{self.dmr_out_prefix}.regions.bed.gz")

    @property
    def dmr_regions_p(self) -> Path:
        return Path(f"{self.dmr_out_prefix}.regions-p.bed.gz")

    @property
    def dmr_slk(self) -> Path:
        return Path(f"{self.dmr_out_prefix}.slk.bed.gz")

    @property
    def dmr_anno_final(self) -> Path:
        return self._out("dmr", f"{self.assoc_var}_dmr_annotated_results.tsv")

    @property
    def dmr_manhattan_plot(self) -> Path:
        return self._out("dmr", f"{self.assoc_var}_dmr_manhattan.jpg")

    # ---------- Annotation cache ----------
    # One cache root (config: annotation.cache_dir), subdivided by source
    # because each source is versioned differently:
    #
    #   <cache>/zhou/<platform>/<zhou_release>/   Zhou Infinium manifest
    #   <cache>/ucsc/<genome>/<cache_tag>/        UCSC tracks
    #   <cache>/eqtm/                             BIOS eQTM + HGNC
    @property
    def ewas_zhou_dir(self) -> Path:
        """Zhou annotation cache, keyed by platform and release tag."""
        return self.anno_cache_dir.joinpath(
            "zhou", self.array_platform, self.zhou_release
        )

    @property
    def eqtm_dir(self) -> Path:
        return self.anno_cache_dir.joinpath("eqtm")

    # KYCG feature sets carried through to the annotated results as one column
    # each. Keyed by the filename prefix in <platform>/KYCG/ (the published
    # names are date-stamped and differ per platform, so they are matched by
    # prefix and saved under the normalised name).
    #
    # Deliberately excluded: CGI (already derived from the UCSC track),
    # InfiniumChemistry / ProbeType / Tetranuc2 / nFlankCG / MetagenePC
    # (probe design rather than biology), and the large multi-record sets
    # HM, TFBSrm, rmsk2 and REMCChromHMM, which are tested for enrichment
    # instead of being flattened into a column.
    KYCG_FEATURE_SETS = {
        "ChromHMM": "chromHMM_state",
        "PMD": "PMD",
        "ABCompartment": "AB_compartment",
        "rmsk1": "repeat_class",
        "ImprintingDMR": "imprinting_DMR",
        "CTCFbind": "CTCF_binding",
        "Blacklist": "ENCODE_blacklist",
    }

    # Knowledgebases tested for enrichment. Deliberately a short list, and
    # much shorter than what a platform publishes (EPIC 17 sets, MSA 32).
    #
    # Chosen to be biological, non-redundant, and not already derivable from
    # the annotated results:
    #   ChromHMM, REMCChromHMM  chromatin state, two references
    #   HM                      histone marks and variants
    #   TFBSrm                  transcription factor binding, 1188 motifs
    #   CTCFbind                methylation-sensitive CTCF binding
    #   CGI                     island / shore / shelf / open sea
    #   MetagenePC              position relative to the gene
    #   ABCompartment           Hi-C A/B compartments
    #   PMD                     partially methylated domains
    #   rmsk1, rmsk2            repeat classes and families
    #   Tetranuc2               WCGW / SCGS sequence context
    #   ImprintingDMR           positive control: allele-specific methylation
    #
    # Large and overlapping sets are all kept because the FDR is computed
    # within each knowledgebase, not across all of them, so no set spends
    # anything the others need.
    #
    # On EPIC these 13 plus the four QC sets below are everything the platform
    # publishes. Other platforms publish more -- MSA has 32 sets, including
    # tissue signatures, CoRSIVs and evolutionary conservation -- and any of
    # them can be named here. Several of these overlap by construction
    # (REMCChromHMM restates ChromHMM from another reference; HM is the histone
    # data ChromHMM states are called from; rmsk2 is rmsk1 at finer resolution;
    # Tetranuc2 shares its context partition with nFlankCG; CGI covers what the
    # CGIposition column already carries), which is a reason to read them as
    # corroborating rather than independent -- not a reason to drop them.
    #
    # Widen with `enrichment: knowledgebases:` in the config -- a list of set
    # names, or "all" for everything the platform publishes.
    KYCG_ENRICHMENT_SETS = (
        "ChromHMM", "REMCChromHMM", "HM", "TFBSrm", "CTCFbind",
        "CGI", "MetagenePC", "ABCompartment", "PMD",
        "rmsk1", "rmsk2", "Tetranuc2", "ImprintingDMR",
    )

    # Design and QC sets. Tested and reported in the same table (with
    # role = "qc"), but plotted separately, because enrichment here is not a
    # finding about biology: it says the hit list tracks probe design,
    # artefact-prone regions or CpG density. Zhou's registry frames the first
    # three as exactly this kind of post-hoc check.
    KYCG_QC_SETS = ("ProbeType", "InfiniumChemistry", "Blacklist", "nFlankCG")

    @property
    def kycg_enrichment_sets(self) -> tuple:
        """Set names tested by enrich_features, or ("all",)."""
        cfg = getattr(self, "_enrich_kbs", None)
        if cfg is None:
            return self.KYCG_ENRICHMENT_SETS
        if isinstance(cfg, str):
            if cfg.strip().lower() == "all":
                return ("all",)
            return tuple(x.strip() for x in cfg.split(",") if x.strip())
        return tuple(str(x).strip() for x in cfg if str(x).strip())

    @property
    def kycg_qc_sets(self) -> tuple:
        """Design/QC set names tested and reported under role = "qc"."""
        cfg = getattr(self, "_enrich_qc_kbs", None)
        if cfg is None:
            return self.KYCG_QC_SETS
        if isinstance(cfg, str):
            if cfg.strip().lower() in ("", "none"):
                return ()
            return tuple(x.strip() for x in cfg.split(",") if x.strip())
        return tuple(str(x).strip() for x in cfg if str(x).strip())

    @property
    def kycg_tested_sets(self) -> tuple:
        """Everything enrich_features tests: biological plus QC."""
        bio = self.kycg_enrichment_sets
        if bio == ("all",):
            return ("all",)
        return tuple(dict.fromkeys(bio + self.kycg_qc_sets))

    @property
    def kycg_download_sets(self) -> tuple:
        """Sets the fetch rule needs: annotation columns plus tested sets.

        Downloading only what is used matters on MSA, which publishes 32 sets;
        TFBSrm and HM alone are hundreds of megabytes.
        """
        tested = self.kycg_tested_sets
        if tested == ("all",):
            return ("all",)
        return tuple(dict.fromkeys(tuple(self.KYCG_FEATURE_SETS) + tested))

    @property
    def kycg_dir(self) -> Path:
        """Normalised KYCG feature sets for the configured platform."""
        return self.ewas_zhou_dir.joinpath("KYCG")

    @property
    def kycg_manifest(self) -> Path:
        return self.ewas_zhou_dir.joinpath("kycg_manifest.tsv")

    @property
    def kycg_registry(self) -> Path:
        """Zhou's knowledgebase registry: what each KYCG set IS.

        One row per set with title, biology, upstream source, citation and
        processing notes. Cached beside the sets so a results directory carries
        the definitions of the features it reports -- a feature label such as
        ABCompartment's "B4" is not self-explanatory.
        """
        return self.kycg_dir.joinpath("knowledgebases.tsv")

    # Last zhou-lab/kycg commit that still carries data/knowledgebases.tsv.
    # Upstream deleted the file in abb75cea (2026-09-20) and now compiles the
    # definitions into the kycg binary (src/kbinfo.h), so there is no longer a
    # TSV on any branch or tag -- v0.5 and v0.6 postdate the removal. The URL
    # used to track `main`, which is how a routine upstream refactor broke the
    # fetch; a commit SHA cannot move. It ships with the kycg CLI rather than
    # the annotation release, so it is independent of zhou_release. Covers
    # every knowledgebase this workflow tests; bump deliberately, if ever.
    KYCG_REGISTRY_COMMIT = "d6df6f36b81234e23c176c72a2fcc287471ea13e"

    @property
    def kycg_registry_url(self) -> str:
        return (
            "https://raw.githubusercontent.com/zhou-lab/kycg/"
            f"{self.KYCG_REGISTRY_COMMIT}/data/knowledgebases.tsv"
        )

    @property
    def probe_ordering(self) -> Path:
        """Row order the KYCG .cm files are aligned to."""
        return self.ewas_zhou_dir.joinpath(f"{self.array_platform}.ordering.tsv.gz")

    @property
    def probe_ordering_url(self) -> str:
        return (
            f"{self.zhou_anno_raw_base}/{self.zhou_release}/"
            f"{self.array_platform}/{self.array_platform}.ordering.tsv.gz"
        )

    @property
    def kycg_api_url(self) -> str:
        return (
            f"{self.zhou_anno_api_base}/{self.array_platform}/KYCG"
            f"?ref={self.zhou_release}"
        )

    @property
    def kycg_raw_base(self) -> str:
        return (
            f"{self.zhou_anno_raw_base}/{self.zhou_release}/"
            f"{self.array_platform}/KYCG"
        )

    @property
    def ewas_gene_manifest_name(self) -> str:
        return (
            f"{self.array_platform}.{self.genome_build}"
            f".manifest.gencode.{self.gencode_release}.tsv.gz"
        )

    @property
    def ewas_gene_manifest(self) -> Path:
        """Zhou gene/promoter annotation table (probeID -> genes, distToTSS)."""
        return self.ewas_zhou_dir.joinpath(self.ewas_gene_manifest_name)

    @property
    def ewas_gene_manifest_url(self) -> str:
        return (
            f"{self.zhou_raw_base}/{self.zhou_release}/Anno/"
            f"{self.array_platform}/{self.ewas_gene_manifest_name}"
        )

    @property
    def bios_eqtm_txt(self) -> Path:
        return self.eqtm_dir.joinpath(os.path.basename(self.bios_eqtm_url))

    @property
    def hgnc_complete_set_txt(self) -> Path:
        return self.eqtm_dir.joinpath(
            os.path.basename(self.hgnc_complete_set_url)
        )

    @property
    def bios_eqtm_annotation(self) -> Path:
        """BIOS eQTM table with HGNC-resolved GRCh38 gene symbols."""
        stem = self.hgnc_complete_set_txt.name
        stem = stem.replace("hgnc_complete_set", "").strip("_").removesuffix(".txt")
        suffix = f"_{stem}" if stem else ""
        return self.eqtm_dir.joinpath(f"eQTM_annotations_BIOS_HGNC{suffix}.tsv")

    @property
    def ewas_annotation_manifest(self) -> Path:
        """Provenance record for the CpG-level annotation sources."""
        return self.ewas_zhou_dir.joinpath("annotation_manifest.tsv")

    # ---------- Enrichment ----------
    @property
    def ewas_atlas_txt(self) -> Path:
        """Cached EWAS Atlas association table (CpG to trait)."""
        return self.anno_cache_dir.joinpath("ewas_atlas", "EWAS_Atlas_associations.tsv")

    @property
    def enrichment_out_dir(self) -> Path:
        return self._out("enrichment")

    @property
    def enrichment_feature_results(self) -> Path:
        return self._out("enrichment", f"{self.assoc_var}_enrichment_features.tsv")

    @property
    def enrichment_pathway_results(self) -> Path:
        return self._out("enrichment", f"{self.assoc_var}_enrichment_pathways.tsv")

    @property
    def enrichment_trait_results(self) -> Path:
        return self._out("enrichment", f"{self.assoc_var}_enrichment_traits.tsv")

    # One summary plot per enrichment analysis, named by kind so the
    # plot_enrichment rule can carry a {kind} wildcard.
    # "features" and "features_qc" plot two slices of the same table, split on
    # its role column: the biological knowledgebases and the design/QC ones.
    # They are kept apart because enrichment in a QC set is a warning about the
    # other results, not a result of its own.
    ENRICHMENT_KINDS = ("features", "features_qc", "pathways", "traits")

    def enrichment_plot(self, kind: str) -> Path:
        return self._out("enrichment", f"{self.assoc_var}_enrichment_{kind}.jpg")

    def enrichment_plot_files(self) -> List[Path]:
        if not self.enrichment_plots:
            return []
        return [self.enrichment_plot(k) for k in self.ENRICHMENT_KINDS]

    def enrichment_table(self, kind: str) -> Path:
        return {
            "features": self.enrichment_feature_results,
            "features_qc": self.enrichment_feature_results,
            "pathways": self.enrichment_pathway_results,
            "traits": self.enrichment_trait_results,
        }[kind]

    @property
    def enrichment_cpg_set(self) -> Path:
        """The foreground/background split actually tested, kept for the record."""
        return self._out("enrichment", f"{self.assoc_var}_enrichment_cpg_sets.tsv")

    # ---------- UCSC track cache (CpG islands, refGene, HGNC) ----------
    @property
    def dmr_anno_resource_dir(self) -> Path:
        return self.anno_cache_dir.joinpath(
            "ucsc", self.genome_build, self.dmr_anno_cache_tag
        )

    @property
    def cpg_island_manifest(self) -> Path:
        """Provenance for the CpG-island track, which an EWAS-only run needs."""
        return self.dmr_anno_resource_dir.joinpath("cpg_island_manifest.tsv")

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
    def dmr_refgene_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath(f"{self.gene_table}.bed.gz")

    @property
    def dmr_cpg_island_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("cpgIslandExt.bed.gz")

    @property
    def dmr_hgnc_bed(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("hgnc.bed.gz")

    @property
    def dmr_annotation_manifest(self) -> Path:
        return self.dmr_anno_resource_dir.joinpath("annotation_manifest.tsv")

    # ---------- Run provenance ----------
    @property
    def provenance_dir(self) -> Path:
        """Per-run record of the command and configuration that produced the results."""
        return self._out("provenance")

    @property
    def provenance_log(self) -> Path:
        """One line per run: run_id, times, status, command."""
        return self.provenance_dir.joinpath("runs.tsv")


# ---------------------------------------------------------------------------
# Run provenance
#
# Written from the Snakefile's onstart/onsuccess/onerror handlers rather than
# from a rule. Three reasons a rule does not work here:
#
#   * A rule's output is cached. On a second invocation the file is already
#     present and up to date, so the rule does not re-run and the recorded
#     command stays stale from the first run -- the opposite of what a
#     provenance record is for.
#   * Making the output unique per run means a timestamp in the path, and the
#     Snakefile is re-parsed by every job subprocess. Each would compute a
#     different timestamp and the target would stop matching.
#   * A rule body runs in a job subprocess, where sys.argv is Snakemake's
#     internal re-invocation ("--target-jobs ... --mode subprocess"), not the
#     command the user typed. onstart runs in the main process, where sys.argv
#     is the real command line.
#
# Consequence worth knowing: these handlers fire only when Snakemake actually
# executes jobs. A dry run, or a run that reports "Nothing to be done", writes
# no record -- correct, since no results were produced, and the record from the
# run that did produce them is already on disk.
# ---------------------------------------------------------------------------

def _git_state(repo_dir: Union[str, Path]) -> Dict[str, str]:
    """Commit, branch and dirty flag for the workflow checkout, if it is one."""
    import subprocess

    def _git(*args: str) -> Optional[str]:
        try:
            out = subprocess.run(
                ["git", "-C", str(repo_dir), *args],
                capture_output=True, text=True, timeout=10,
            )
            return out.stdout.strip() if out.returncode == 0 else None
        except Exception:
            return None

    if _git("rev-parse", "--is-inside-work-tree") != "true":
        return {"commit": "not_a_git_checkout"}
    status = _git("status", "--porcelain")
    return {
        "commit": _git("rev-parse", "HEAD") or "unknown",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD") or "unknown",
        "describe": _git("describe", "--tags", "--always", "--dirty") or "unknown",
        "dirty": "yes" if status else "no",
        "uncommitted_files": str(len(status.splitlines())) if status else "0",
    }


def record_run_provenance(CW, workflow, config, argv) -> Optional[Path]:
    """Snapshot the command and configuration for this run.

    Creates <out_directory>/provenance/<run_id>/ containing:
      command.txt          the exact command line, and each argument on its own
                           line so long invocations stay readable
      config_snapshot.yml  verbatim copy of every configuration file used
      config_resolved.yml  the merged configuration Snakemake actually ran with,
                           including any --config overrides
      run_info.yml         Snakemake and Python versions, workflow git state,
                           host, user, working directory and start time

    Never raises: a provenance failure warns and returns None rather than
    taking down an analysis run.
    """
    import getpass
    import platform
    import shutil
    import socket
    import sys
    from datetime import datetime, timezone

    try:
        import yaml

        started = datetime.now(timezone.utc).astimezone()
        base = CW.provenance_dir
        run_id = started.strftime("%Y%m%dT%H%M%S")
        run_dir = base.joinpath(run_id)
        suffix = 2
        while run_dir.exists():
            run_dir = base.joinpath(f"{run_id}_{suffix}")
            suffix += 1
        run_dir.mkdir(parents=True)

        # --- the command ---
        command = " ".join(shlex_quote(a) for a in argv)
        with open(run_dir.joinpath("command.txt"), "w") as fh:
            fh.write("# Command that produced the results in this directory.\n")
            fh.write(f"# Run {run_dir.name}, started {started.isoformat()}\n")
            fh.write(f"# Working directory: {os.getcwd()}\n\n")
            fh.write(command + "\n\n")
            fh.write("# One argument per line:\n")
            for a in argv:
                fh.write(f"#   {a}\n")

        # --- the configuration files, verbatim ---
        seen, copied = set(), []
        for cf in (workflow.configfiles or []):
            src = Path(str(cf))
            if not src.is_file():
                continue
            key = src.resolve()
            if key in seen:
                continue
            seen.add(key)
            dest_name = ("config_snapshot.yml" if not copied
                         else f"config_snapshot.{len(copied) + 1}.{src.name}")
            shutil.copyfile(src, run_dir.joinpath(dest_name))
            copied.append({"source": str(key), "saved_as": dest_name})

        # --- the configuration actually used, after --config overrides ---
        with open(run_dir.joinpath("config_resolved.yml"), "w") as fh:
            fh.write("# Merged configuration Snakemake ran with, including any\n")
            fh.write("# --config overrides. This, not config_snapshot.yml, is\n")
            fh.write("# what the workflow actually saw.\n")
            yaml.safe_dump(_plain(dict(config)), fh, default_flow_style=False,
                           sort_keys=False)

        # --- everything else ---
        try:
            import snakemake
            smk_version = getattr(snakemake, "__version__", "unknown")
        except Exception:
            smk_version = "unknown"

        info = {
            "run_id": run_dir.name,
            "started": started.isoformat(),
            "status": "running",
            "command": command,
            "working_directory": os.getcwd(),
            "snakefile": str(getattr(workflow, "main_snakefile", "unknown")),
            "config_files": copied,
            "snakemake_version": smk_version,
            "python_version": sys.version.split()[0],
            "host": socket.gethostname(),
            "user": _safe(getpass.getuser),
            "platform": platform.platform(),
            "workflow_git": _git_state(
                Path(str(getattr(workflow, "main_snakefile", "."))).parent),
        }
        with open(run_dir.joinpath("run_info.yml"), "w") as fh:
            yaml.safe_dump(info, fh, default_flow_style=False, sort_keys=False)

        # --- append-only index across runs ---
        log = CW.provenance_log
        if not log.exists():
            with open(log, "w") as fh:
                fh.write("run_id\tstarted\tended\tstatus\tcommand\n")
        with open(log, "a") as fh:
            fh.write(f"{run_dir.name}\t{started.isoformat()}\t\trunning\t{command}\n")

        print(f"[provenance] recording this run in {run_dir}")
        return run_dir
    except Exception as exc:  # never fail a run over bookkeeping
        print(f"[provenance] WARNING: could not record run provenance: {exc}")
        return None


def finalize_run_provenance(run_dir: Optional[Path], status: str) -> None:
    """Record the outcome once the workflow finishes."""
    if run_dir is None:
        return
    try:
        import yaml
        from datetime import datetime, timezone

        ended = datetime.now(timezone.utc).astimezone()
        info_path = Path(run_dir).joinpath("run_info.yml")
        info = {}
        if info_path.is_file():
            with open(info_path) as fh:
                info = yaml.safe_load(fh) or {}
        info["status"] = status
        info["ended"] = ended.isoformat()
        started = info.get("started")
        if started:
            try:
                info["duration_seconds"] = round(
                    (ended - datetime.fromisoformat(started)).total_seconds(), 1)
            except Exception:
                pass
        with open(info_path, "w") as fh:
            yaml.safe_dump(info, fh, default_flow_style=False, sort_keys=False)

        # rewrite this run's line in the index
        log = Path(run_dir).parent.joinpath("runs.tsv")
        if log.is_file():
            lines = log.read_text().splitlines(keepends=True)
            run_id = Path(run_dir).name
            for i, line in enumerate(lines):
                if line.startswith(run_id + "\t"):
                    parts = line.rstrip("\n").split("\t")
                    while len(parts) < 5:
                        parts.append("")
                    parts[2], parts[3] = ended.isoformat(), status
                    lines[i] = "\t".join(parts) + "\n"
                    break
            log.write_text("".join(lines))
    except Exception as exc:
        print(f"[provenance] WARNING: could not finalize run provenance: {exc}")


def shlex_quote(s: str) -> str:
    import shlex
    return shlex.quote(str(s))


def _safe(fn):
    try:
        return fn()
    except Exception:
        return "unknown"


def _plain(obj):
    """Make a config tree safe for yaml.safe_dump (Paths, sets, numpy scalars)."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_plain(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)
