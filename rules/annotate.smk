localrules: get_annotation_data, fetch_kycg_features

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
        platform = CW.array_platform,
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

rule fetch_kycg_features:
    # KYCG feature sets from the coherent release (zhou-lab/InfiniumAnnotation),
    # which is a different repo from the gene manifest above. These are bit- or
    # state-packed and row-aligned to <platform>.ordering.tsv.gz, so the
    # ordering file is fetched alongside them.
    #
    # Published filenames are date-stamped and the dates differ per platform
    # (EPIC ships 27 files, EPICv2 25, MSA 45), so files are matched by prefix
    # and saved under a normalised name. A set that a platform does not publish
    # is simply skipped; annotation.R adds a column only for what is present.
    output:
        kycg = directory(CW.kycg_dir),
        ordering = protected(CW.probe_ordering),
        manifest = CW.kycg_manifest
    params:
        api_url = CW.kycg_api_url,
        raw_base = CW.kycg_raw_base,
        ordering_url = CW.probe_ordering_url,
        platform = CW.array_platform,
        release = CW.zhou_release,
        prefixes = " ".join(CW.KYCG_FEATURE_SETS)
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        mkdir -p {output.kycg}

        wget -O {output.ordering}.tmp {params.ordering_url}
        mv {output.ordering}.tmp {output.ordering}

        # One unauthenticated call to the GitHub contents API resolves the
        # date-stamped filenames for this platform and release.
        listing=$(mktemp)
        curl -sSL --retry 3 --max-time 120 "{params.api_url}" \
          | grep -oE '"name": "[^"]+\\.cm"' | sed 's/.*: "//; s/"//' | sort -u > "$listing"
        if [ ! -s "$listing" ]; then
          echo "ERROR: no KYCG files listed for {params.platform} at {params.release}." >&2
          echo "       Checked {params.api_url}" >&2
          rm -f "$listing"; exit 1
        fi

        # Every set is fetched, not just the ones that become columns: the
        # enrichment rules test against all of them, including the ones too
        # large to flatten into a column (TFBSrm alone has 1188 records).
        # The .cm.idx sidecar carries the record names -- without it yame
        # reports bare indices -- so it is fetched alongside each set.
        {{
          echo -e "feature_set\tpublished_file\tplatform\trelease\tannotation_column\tsource\tcreated"
          while read -r src; do
            [ -z "$src" ] && continue
            p=$(echo "$src" | sed -E 's/(\\.[0-9]{{8}})?\\.cm$//')
            curl -sSL --retry 3 --max-time 900 -o "{output.kycg}/${{p}}.cm" \
              "{params.raw_base}/${{src}}"
            curl -sSL --retry 3 --max-time 300 -o "{output.kycg}/${{p}}.cm.idx" \
              "{params.raw_base}/${{src}}.idx" || rm -f "{output.kycg}/${{p}}.cm.idx"
            col="enrichment_only"
            case " {params.prefixes} " in *" ${{p}} "*) col="yes";; esac
            echo -e "${{p}}\t${{src}}\t{params.platform}\t{params.release}\t${{col}}\t{params.raw_base}/${{src}}\t$(date -Iseconds)"
          done < "$listing"
          echo -e "probe_ordering\t{params.platform}.ordering.tsv.gz\t{params.platform}\t{params.release}\tNA\t{params.ordering_url}\t$(date -Iseconds)"
        }} > {output.manifest}

        rm -f "$listing"
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
        kycg_dir = rules.fetch_kycg_features.output.kycg,
        ordering = rules.fetch_kycg_features.output.ordering,
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
        --gene-anno {input.gene_file} \
        --cpg-islands {input.cpg_islands} \
        --eQTM-anno {input.eqtm_file} \
        --kycg-dir {input.kycg_dir} \
        --probe-order {input.ordering} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} \
        --out-type {params.o_type}
        """