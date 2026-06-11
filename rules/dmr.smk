localrules: fetch_dmr_annotation_cache

rule make_bed:
    input: 
        in_file = annotated_results,
        script = "scripts/make_bed.R"
    params:
        o_dir = OUT_DIR,
        strat = STRATIFIED,
        assoc = ASSOC
    output: 
        results_bed
    conda:
        "../envs/dmr.yaml"
    shell:
        """
        Rscript {input.script} \
        --results {input.in_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} 
        """

rule run_dmr:
    input:
        in_file = results_bed
    params:
        o_prefix = OUT_DIR +  "/dmr/" + ASSOC + "_ewas",
        min_p = CW.min_pval,
        win_sz = CW.win_size,
        region_filter = CW.region_filter
    output: 
        acf = dmr_acf,
        args = dmr_args,
        fdr = dmr_fdr,
        regions = dmr_regions,
        slk = dmr_slk
    conda:
        "../envs/dmr.yaml"
    shell:
        """
        comb-p pipeline \
		--seed {params.min_p} \
		--dist {params.win_sz} \
		-p {params.o_prefix} \
		--region-filter-p {params.region_filter} \
		{input.in_file} 
        """

rule fetch_dmr_annotation_cache:
    output:
        refgene_txt = dmr_refgene_txt,
        cpg_txt = dmr_cpg_island_txt,
        hgnc_bb = dmr_hgnc_bb,
        refgene_bed = dmr_refgene_bed,
        cpg_bed = dmr_cpg_island_bed,
        hgnc_bed = dmr_hgnc_bed,
        manifest = dmr_annotation_manifest
    params:
        genome = CW.genome_build,
        gene_table = CW.gene_table,
        ucsc_database_base = CW.ucsc_database_base,
        ucsc_gbdb_base = CW.ucsc_gbdb_base,
        cache_tag = CW.dmr_anno_cache_tag
    conda:
        "../envs/dmr_annotation.yaml"
    shell:
        r"""
        set -euo pipefail

        cache_dir=$(dirname {output.refgene_txt})
        mkdir -p "$cache_dir"

        wget -O {output.refgene_txt}.tmp \
          {params.ucsc_database_base}/{params.genome}/database/{params.gene_table}.txt.gz
        mv {output.refgene_txt}.tmp {output.refgene_txt}

        wget -O {output.cpg_txt}.tmp \
          {params.ucsc_database_base}/{params.genome}/database/cpgIslandExt.txt.gz
        mv {output.cpg_txt}.tmp {output.cpg_txt}

        wget -O {output.hgnc_bb}.tmp \
          {params.ucsc_gbdb_base}/{params.genome}/hgnc/hgnc.bb
        mv {output.hgnc_bb}.tmp {output.hgnc_bb}

        # refGene fields:
        # 1 bin, 2 transcript/name, 3 chrom, 4 strand, 5 txStart, 6 txEnd, ... 13 name2/gene symbol
        zcat {output.refgene_txt} \
          | awk 'BEGIN{{FS=OFS="\t"}} {{print $3,$5,$6,$13,$2,$4}}' \
          | sort -k1,1 -k2,2n \
          | gzip -c > {output.refgene_bed}.tmp
        mv {output.refgene_bed}.tmp {output.refgene_bed}

        # cpgIslandExt fields:
        # 1 bin, 2 chrom, 3 chromStart, 4 chromEnd, 5 name, ...
        zcat {output.cpg_txt} \
          | awk 'BEGIN{{FS=OFS="\t"}} {{print $2,$3,$4,$5}}' \
          | sort -k1,1 -k2,2n \
          | gzip -c > {output.cpg_bed}.tmp
        mv {output.cpg_bed}.tmp {output.cpg_bed}

        # UCSC HGNC BigBed fields:
        # 1 chrom, 2 chromStart, 3 chromEnd, 4 HGNC ID,
        # 10 symbol, 11 geneName, 12 locus_group, 13 locus_type
        bigBedToBed {output.hgnc_bb} stdout \
          | awk 'BEGIN{{FS=OFS="\t"}} NF >= 14 {{print $1,$2,$3,$10,$4,$11,$12,$13}}' \
          | sort -k1,1 -k2,2n \
          | gzip -c > {output.hgnc_bed}.tmp
        mv {output.hgnc_bed}.tmp {output.hgnc_bed}

        {{
          echo -e "resource\tgenome_build\tcache_tag\tsource\tcreated"
          echo -e "{params.gene_table}\t{params.genome}\t{params.cache_tag}\t{params.ucsc_database_base}/{params.genome}/database/{params.gene_table}.txt.gz\t$(date -Iseconds)"
          echo -e "cpgIslandExt\t{params.genome}\t{params.cache_tag}\t{params.ucsc_database_base}/{params.genome}/database/cpgIslandExt.txt.gz\t$(date -Iseconds)"
          echo -e "hgnc_bigbed\t{params.genome}\t{params.cache_tag}\t{params.ucsc_gbdb_base}/{params.genome}/hgnc/hgnc.bb\t$(date -Iseconds)"
        }} > {output.manifest}
        """

rule annotate_dmrs:
    input:
        dmr_regions = rules.run_dmr.output.regions,
        ewas_bed = rules.make_bed.output,
        hgnc = rules.fetch_dmr_annotation_cache.output.hgnc_bed,
        refgene = rules.fetch_dmr_annotation_cache.output.refgene_bed,
        cpgIslandExt = rules.fetch_dmr_annotation_cache.output.cpg_bed
    params:
        o_prefix = OUT_DIR +  "/dmr",
        assoc = ASSOC
    output:
        dmr_anno_final
    conda:
        "../envs/dmr.yaml"
    shell:
        """
        Rscript scripts/dmr_annotation.R \
            --dmr-regions {input.dmr_regions} \
            --ewas-bed {input.ewas_bed} \
            --hgnc {input.hgnc} \
            --refgene {input.refgene} \
            --cpgIslandExt {input.cpgIslandExt} \
            --out-dir {params.o_prefix} \
            --assoc {params.assoc}
        """

rule plot_dmrs:
    input:
        slk = rules.run_dmr.output.slk,
        dmr_regions = rules.run_dmr.output.regions
    params:
        o_prefix = OUT_DIR +  "/dmr",
        assoc = ASSOC
    output:
        dmr_manhattan
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript scripts/dmr_plot.R \
            --slk-file {input.slk} \
            --regions-file {input.dmr_regions} \
            --out-dir {params.o_prefix} \
            --assoc {params.assoc} \
            --min-probes 2 \
            --max-y -1 \
            --make-zoom no \
            --zoom-padding 2000 \
            --cluster-gap 3000 \
            --point-jitter-bp 25 \
            --zoom-midlines yes \
            --genome-build hg38 \
            --gene-label-size 3.0 \
            --make-combined no \
            --combined-formats pdf \
            --combined-width-cm 18.3 \
            --combined-region-height-cm 10 \
            --panel-label-size 16
        """