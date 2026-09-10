localrules: get_annotation_data

def get_file(wildcards):
    if config["stratified_ewas"] == "yes":
        in_file = CW.meta_analysis_results
    else:
        in_file = CW.bacon_results
    return(in_file)

rule get_annotation_data:
    # Zhou lab Infinium annotation. The zhouserver.research.chop.edu paths this
    # rule used to wget are retired; the tables now live in the GitHub repo
    # zhou-lab/InfiniumAnnotationData and are pinned here by release tag
    # (config: ewas_annotation.zhou_release) rather than tracking main.
    #
    # Only the gene/promoter manifest is fetched. The companion
    # EPIC.hg38.commonsnp.tsv.gz is a frozen 1000 Genomes Phase 3 / dbSNP 151
    # probe-to-variant lookup and is deliberately not used: SNP-affected probes
    # are expected to be masked before the EWAS is run.
    output:
        gene_manifest = protected(CW.ewas_gene_manifest),
        bios_eqtm = protected(CW.bios_eqtm_txt),
        hgnc = protected(CW.hgnc_complete_set_txt),
        manifest = CW.ewas_annotation_manifest
    params:
        gene_url = CW.ewas_gene_manifest_url,
        eqtm_url = CW.bios_eqtm_url,
        hgnc_url = CW.hgnc_complete_set_url,
        platform = CW.ewas_anno_platform,
        genome = CW.genome_build,
        release = CW.zhou_release,
        gencode = CW.gencode_release
    shell:
        """
        wget -O {output.gene_manifest}.tmp {params.gene_url}
        mv {output.gene_manifest}.tmp {output.gene_manifest}

        wget -O {output.bios_eqtm}.tmp {params.eqtm_url}
        mv {output.bios_eqtm}.tmp {output.bios_eqtm}

        wget -O {output.hgnc}.tmp {params.hgnc_url}
        mv {output.hgnc}.tmp {output.hgnc}

        {{
          echo -e "resource\tplatform\tgenome_build\trelease\tsource\tcreated"
          echo -e "zhou_gene_manifest_gencode_{params.gencode}\t{params.platform}\t{params.genome}\t{params.release}\t{params.gene_url}\t$(date -Iseconds)"
          echo -e "bios_eqtm\tNA\tNA\tNA\t{params.eqtm_url}\t$(date -Iseconds)"
          echo -e "hgnc_complete_set\tNA\tNA\tNA\t{params.hgnc_url}\t$(date -Iseconds)"
        }} > {output.manifest}
        """

rule prep_bios_eqtm_annotation:
    input:
        bios_eqtm = rules.get_annotation_data.output.bios_eqtm,
        hgnc = rules.get_annotation_data.output.hgnc
    output:
        annotation = CW.bios_eqtm_annotation
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
        gene_file = rules.get_annotation_data.output.gene_manifest,
        # Produced by fetch_cpg_island_cache in rules/dmr.smk. Referenced by
        # path rather than via rules.* because dmr.smk is included after this
        # file. GENCODE v41 dropped the CGI and CGIposition columns that v36
        # carried, so islands are called from the UCSC track instead -- which
        # is where Zhou's own CGI annotation came from.
        cpg_islands = CW.dmr_cpg_island_bed,
        eqtm_file = rules.prep_bios_eqtm_annotation.output.annotation,
        script = "scripts/annotation.R"
    params:
        o_dir = CW.out_dir,
        strat = CW.stratified,
        assoc = CW.assoc_var,
        o_type = CW.out_type,
        shore_bp = CW.cpg_island_shore_bp,
        shelf_bp = CW.cpg_island_shelf_bp
    output: 
        CW.annotated_results
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --gene-anno {input.gene_file} \
        --cpg-islands {input.cpg_islands} \
        --eQTM-anno {input.eqtm_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} \
        --out-type {params.o_type} \
        --shore-bp {params.shore_bp} \
        --shelf-bp {params.shelf_bp}
        """