

rule get_fst:
    input:
        vcf = "data/gwasrapidd/{date}/vcfs/no_dups/{efo_id}.vcf.gz",
        pop_file = config["local_pop_file_plink"]
    output:
        "data/gwasrapidd/{date}/plink/fst/{efo_id}.fst"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(str(output)),
        pref = lambda wildcards, output: os.path.splitext(str(output))[0]
    container:
        config["plink1.9"]
    shell:
        """
        mkdir -p {params.output_dir} ;
        plink \
            --vcf {input.vcf} \
            --double-id \
            --fst \
            --within {input.pop_file} \
            --out {params.pref} \
            --vcf-half-call missing
        """

#rule get_risk_alleles:
#    input:
#        "data/gwasrapidd/{date}/assocations_raw/{efo_id}.rds"
#    output:
#        "data/gwasrapidd/{date}/assocations_risk_alleles/{efo_id}.txt"
#    container:
#        "envs/r_4.1.0.yaml"
#    script:
#        "scripts/get_risk_alleles.R"

rule recode_012:
    input:
        vcf = "data/gwasrapidd/{date}/vcfs/no_dups/{efo_id}.vcf.gz"
    output:
        "data/gwasrapidd/{date}/plink/recode_012/{efo_id}.traw"
    params:
        pref = lambda wildcards, output: os.path.splitext(str(output))[0]
    container:
        config["plink1.9"]
    shell:
        """
        plink \
            --vcf {input.vcf} \
            --double-id \
            --recode A-transpose \
            --out {params.pref}
        """

rule get_p_values:
    input:
        "data/gwasrapidd/{date}/associations_raw/{efo_id}.rds"
    output:
        "data/gwasrapidd/{date}/p-values/{efo_id}.txt"
    container:
        config["r_container"]
    script:
        "scripts/get_p_values.R"


rule clump_snps:
    input:
        vcf = "data/gwasrapidd/{date}/vcfs/no_dups/{efo_id}.vcf.gz",
        p_values = "data/gwasrapidd/{date}/p-values/{efo_id}.txt"
    output:
        "data/gwasrapidd/{date}/plink/clumped/{efo_id}.clumped"
    params:
        pref = lambda wildcards, output: os.path.splitext(str(output))[0]
    container:
        config["plink1.9"]
    shell:
        """
        plink \
          --vcf {input.vcf} \
          --clump {input.p_values} \
          --clump-p1 0.00000001 \
          --clump-p2 0.00000001 \
          --clump-r2 {config[clump_r2]} \
          --clump-kb {config[clump_kb]} \
          --out {params.pref}     
        """
# Note: some variants missing from the main dataset, e.g.
# Warning: 'rs199921354' is missing from the main dataset, and is a top variant.
# 1 more top variant ID missing; see log file.

#rule recode:
#    input:
#        os.path.join(config["working_dir"], "data/gwasrapidd/{date}/vcfs/{efo_id}.vcf.gz")
#    output:
#        os.path.join(config["working_dir"], "data/gwasrapidd/{date}/plink/recode/{efo_id}.vcf.gz")
