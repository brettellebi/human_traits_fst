rule download_1KG_38_annotated:
    params:
        input = os.path.join(config["ftp_dir_1kg_38_annotated"], "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP.vcf{gz_ext}")
    output:
        os.path.join(config["working_dir"], "vcfs/1kg/20150319/chrs/{chr}.vcf{gz_ext}")
    shell:
        """
        wget -O {output} {params.input}
        """

# When trying to merge, get the following error:
#Caused by: htsjdk.tribble.TribbleException$InvalidHeader: Your input file has a malformed header: Unclosed quote in header line value <ID=ssID,Number=A,Type=String,Description=dbSNP ssID of the allele">
# Therefore created a new header file with:
# bcftools view /hps/nobackup/birney/users/ian/hmn_fst/vcfs/1kg/20150319/chrs/1.vcf.gz | grep "#" > /hps/nobackup/birney/users/ian/hmn_fst/vcfs/1kg/20150319/new_header.vcf
# Then manually inserted the missing quote. Use that file to re-header all VCFs before merging

rule fix_vcf_headers:
    input:
        os.path.join(config["working_dir"], "vcfs/1kg/20150319/chrs/{chr}.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/1kg/20150319/reheaded/{chr}.vcf.gz")
    singularity:
        config["bcftools"]
    shell:
        """
        bcftools reheader \
            --header {config[new_header]} \
            --output {output} \
            {input}
        """

rule get_population_file:
    input:
        FTP.remote(config["ftp_pop_file"], keep_local = True)
    output:
        config["local_pop_file"]
    run:
        pop_file = pd.read_excel(input[0], sheet_name = "Sample Info")
        pop_file = pop_file.loc[:, ['Sample', 'Population']]
        pop_file.to_csv(output[0], index = False)

rule get_population_file_plink:
    input:
        FTP.remote(config["ftp_pop_file"], keep_local = True)
    output:
        config["local_pop_file_plink"]
    run:
        pop_file = pd.read_excel(input[0], sheet_name = "Sample Info")
        pop_file = pop_file.loc[:, ['Sample', 'Population']]
        # Create second column of samples as IIDs
        pop_file['Sample_2'] = pop_file['Sample']
        # Rename columns
        pop_file = pop_file.rename(columns = {"Sample" : "FID", "Sample_2" : "IID", "Population" : "CLUSTER"})
        # Re-order columns
        pop_file = pop_file[['FID', 'IID', 'CLUSTER']]
        # Write to file
        pop_file.to_csv(output[0], sep = "\t", header = False, index = False)
