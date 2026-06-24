def get_file(wildcards):
    if config["stratified_ewas"] == "yes":
        in_file = CW.meta_analysis_results
    else:
        in_file = CW.bacon_results
    return(in_file)

rule get_annotation_data:
    output:
        epic_hg38 = protected("resources/ewas_annotation/EPIC_hg38.tsv.gz"),
        epic_snp = protected("resources/ewas_annotation/EPIC_snp_key.tsv.gz"),
        bios_eqtm = protected("resources/ewas_annotation/2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt"),
        hgnc = protected("resources/ewas_annotation/hgnc_complete_set_2025-07-01.txt")
    shell:
        """
        wget https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/EPIC/EPIC.hg38.manifest.gencode.v36.tsv.gz \
        -O {output.epic_hg38}
        wget https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.commonsnp.tsv.gz \
        -O {output.epic_snp}
        wget https://molgenis26.gcc.rug.nl/downloads/biosqtlbrowser/2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt \
        -O {output.bios_eqtm}
        wget https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2025-07-01.txt \
        -O {output.hgnc}
        """

rule prep_bios_eqtm_annotation:
    input:
        bios_eqtm = rules.get_annotation_data.output.bios_eqtm,
        hgnc = rules.get_annotation_data.output.hgnc
    output:
        annotation = "resources/ewas_annotation/eQTM_annotations_BIOS_HGNC_2025-07-01.tsv"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript scripts/prepare_bios_eqtm_annotation.R \
        --bios-eqtm {input.bios_eqtm} \
        --hgnc-complete-set {input.hgnc} \
        --output {output.annotation}
        """

rule add_annotation:
    input:
        in_file = get_file,
        hg38_file = rules.get_annotation_data.output.epic_hg38,
        snp_file = rules.get_annotation_data.output.epic_snp,
        eqtm_file = rules.prep_bios_eqtm_annotation.output.annotation,
        script = "scripts/annotation.R"
    params:
        o_dir = CW.out_dir,
        strat = CW.stratified,
        assoc = CW.assoc_var,
        o_type = CW.out_type
    output: 
        CW.annotated_results
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --EPIC-anno {input.hg38_file} \
        --EPIC-snp {input.snp_file} \
        --eQTM-anno {input.eqtm_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} \
        --out-type {params.o_type}
        """