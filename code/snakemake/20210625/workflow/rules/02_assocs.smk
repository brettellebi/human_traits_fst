
rule get_associations:
    output:
        "data/gwasrapidd/{date}/associations_raw/{efo_id}.rds"
    params:
        date = lambda wildcards: wildcards.date,
        efo_id = lambda wildcards: wildcards.efo_id,
        output_dir = lambda wildcards, output: os.path.dirname(str(output))
    container:
        config["R"]
    shell:
        """
        mkdir -p {params.output_dir} ;
        {config[rscript]} {config[get_associations_script]} {params.efo_id} {output}
        """

rule get_studies:
    input:
        "data/gwasrapidd/{date}/associations_raw/{efo_id}.rds"
    output:
        key = "data/gwasrapidd/{date}/studies_key/{efo_id}.rds",
        studies = "data/gwasrapidd/{date}/studies_raw/{efo_id}.rds"
    container:
        config["R"]
    script:
        "scripts/get_studies.R"

rule process_associations:
    input:
        "data/gwasrapidd/{date}/associations_raw/{efo_id}.rds"
    output:
        "data/gwasrapidd/{date}/associations_clean/{efo_id}.csv"
    params:
        efo_id = lambda wildcards: wildcards.efo_id
    container:
        config["R"]
    script:
        "scripts/process_associations.R"

rule consolidate_associations:
    input:
        expand("data/gwasrapidd/{date}/associations_clean/{efo_id}.csv",
            date = DATE_OF_COLLECTION,
            efo_id = EFO_IDS)
    output:
        "data/gwasrapidd/{date}/final/final.csv.gz"
    run:
        out = pd.concat([pd.read_csv(file) for file in input])
        # Convert to csv and export
        out.to_csv(output[0], index = False)

rule get_snp_ids:
    input:
        "data/gwasrapidd/{date}/assocations_raw/{efo_id}.rds"
    output:
        "data/gwasrapidd/{date}/assocations_snp_ids/{efo_id}.txt"
    container:
        config["R"]
    script:
        "scripts/get_snp_ids.R"

rule extract_gtypes:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/1kg/20150319/reheaded/{chr}.vcf.gz"),
        snps = "data/gwasrapidd/{date}/assocations_snp_ids/{efo_id}.txt"
    output:
        os.path.join(config["working_dir"], "vcfs/1kg/20150319/filtered/{date}/{efo_id}/by_chr/{chr}.vcf.gz")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --include ID=@{input.snps} \
            --output-type z \
            --output-file {output} \
            {input.vcf}
        """

rule merge_gtypes:
    input:
        expand(os.path.join(config["working_dir"], "vcfs/1kg/20150319/filtered/{{date}}/{{efo_id}}/by_chr/{chr}.vcf.gz"),
            chr = CHRS)
    output:
        vcf = "data/gwasrapidd/{date}/vcfs/{efo_id}.vcf.gz"
    container:
        config["bcftools"]
    shell:
        """
        bcftools concat \
            --output {output.vcf} \
            --output-type z \
            {input}
        """

# Remove duplicated variants to avoid problems downstream
rule get_duplicated_sites:
    input:
        vcf = "data/gwasrapidd/{date}/vcfs/{efo_id}.vcf.gz"
    output:
        dup_sites = "data/gwasrapidd/{date}/dup_sites/{efo_id}.txt",
        vcf = "data/gwasrapidd/{date}/vcfs/no_dups/{efo_id}.vcf.gz"
    container:
        config["bcftools"]
    shell:
        """
        bcftools view {input.vcf} | grep -v '^#' | cut -f 3 | sort | uniq -d > {output.dup_sites} ;
        bcftools view \
            --exclude ID=@{output.dup_sites} \
            --output-type z \
            --output-file {output.vcf} \
            {input.vcf}
        """
