rule stratify_data:
    input:
        script = "scripts/stratify.R",
        fxns = "scripts/fxns/stratify_fxns.R",
        pheno_file = PHENO,
        methyl_file = MVALS
    params:
        strat_vars = ','.join(STRAT_VARS),
        o_dir = OUT_DIR
    threads:
        N_WORKERS
    output: 
        temp(strat_pheno_files),
        temp(strat_mvals_files)
    conda:
        "../envs/ewas.yaml"    
    shell:
        """
        Rscript {input.script} \
        --pheno {input.pheno_file} \
        --methyl {input.methyl_file} \
        --stratify {params.strat_vars} \
        --out-dir {params.o_dir} \
        --threads {threads}
        """

rule run_ewas_group:
    input:
        script = "scripts/ewas.R",
        pheno_file = lambda wildcards: str(CW.group_pheno(wildcards.group)),
        methyl_file= lambda wildcards: str(CW.group_mvals(wildcards.group))
    params:
        assoc_var = ASSOC,
        stratified = STRATIFIED,
        cs = config["chunk_size"],
        pt = config["processing_type"],
        n_workers = N_WORKERS,
        o_dir = lambda wildcards: str(CW.group_dir(wildcards.group)),
        o_type = OUT_TYPE,
        o_prefix = lambda wildcards: wildcards.group
    output:
        ewas_results = str(CW.group_ewas_results("{group}"))
    log:
        "log/{group}_ewas.log"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        export R_PROGRESSR_ENABLE=TRUE 
        Rscript {input.script} \
        --pheno {input.pheno_file} \
        --methyl {input.methyl_file} \
        --assoc {params.assoc_var} \
        --stratified {params.stratified} \
        --chunk-size {params.cs} \
        --processing-type {params.pt} \
        --workers {params.n_workers} \
        --out-dir {params.o_dir} \
        --out-type {params.o_type} \
        --out-prefix {params.o_prefix} \
        > {log} 2>&1
        """

rule run_bacon_group:
    input:
        in_file = rules.run_ewas_group.output.ewas_results,
        script = "scripts/run_bacon.R"
    params:
        o_dir = lambda wildcards: str(CW.group_dir(wildcards.group)),
        o_type = OUT_TYPE,
        o_prefix = lambda wildcards: wildcards.group
    output:
        bacon_results = str(CW.group_bacon_results("{group}")),
        plot_trace = str(CW.group_bacon_plot("{group}", "traces")),
        plot_post = str(CW.group_bacon_plot("{group}", "posteriors")),
        plot_fit = str(CW.group_bacon_plot("{group}", "fit")),
        plot_qqs = str(CW.group_bacon_plot("{group}", "qqs"))
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {params.o_dir} \
        --out-prefix {params.o_prefix} \
        --out-type {params.o_type} \
        """

rule install_metal:
    output:
        metal_bin = "software/metal/build/bin/metal"
    conda:
        "../envs/metal-build.yaml"
    shell:
        """
      	cd software/metal
	    tar -xzf METAL.tar.gz --strip-components=1
	    sed -i 's/cmake_minimum_required(VERSION [0-9.]\\+)/cmake_minimum_required(VERSION 3.5...3.27)/' CMakeLists.txt
	    sed -i '/add_executable(metal/a target_include_directories(metal PRIVATE ${{ZLIB_INCLUDE_DIRS}})' metal/CMakeLists.txt

        mkdir -p build
        cd build

        # Configure with conda compilers and explicit zlib paths
        cmake .. \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_C_COMPILER=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc \
          -DCMAKE_CXX_COMPILER=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++ \
          -DZLIB_LIBRARY=$CONDA_PREFIX/lib/libz.so \
          -DZLIB_INCLUDE_DIR=$CONDA_PREFIX/include

        # Build and install
        make
        make test
        make install
        """

rule make_metal_script:
    input:
        script = "scripts/metal_cmd.sh",
        in_files = expand(rules.run_bacon_group.output.bacon_results, group=GROUPS)
    params:
        out_prefix = str(CW.metal_out_prefix)
    output:
        metal_script = str(CW.metal_command_script)
    shell:
        "bash {input.script} {output.metal_script} {params.out_prefix} {input.in_files}"
    

rule run_metal:
    input:
        metal = rules.install_metal.output.metal_bin,
        bacon_group_results = expand(rules.run_bacon_group.output.bacon_results, group = GROUPS),
        script = rules.make_metal_script.output.metal_script
    output:
        meta_analysis_results
    shell: 
        "{input.metal} {input.script}"
