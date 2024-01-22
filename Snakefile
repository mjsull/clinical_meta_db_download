configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]
gtdb_release: config["gtdb_release"]

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

rule get_gtdb_taxfiles:
    output:
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz"
    shell:
        "wget -O {output.ar53_tax} https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_taxonomy_r214.tsv.gz && " 
        "wget -O {output.bac120_tax} https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_taxonomy_r214.tsv.gz"

rule split_gtdb_accessions:
    input:
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz"
    params:
        datasets_binary = config["datasets_binary"]
    output:
        accession_list = "data/gtdb_accessions",
        gtdb_list = expand("data/gtdb.{num}.list", num=["{:02d}".format(x) for x in range(50)])
    shell:
        "zcat {input.ar53_tax} {input.bac120_tax} | awk -F'\t' '{{print substr($1, 4); }}' > {output.accession_list} && "
        "split -n l/50 -d --additional-suffix=.list data/gtdb_accessions data/gtdb."


rule download_gtdb:
    input:
        gtdb_list =  "data/gtdb.{num}.list"
    output:
        ncbi_file = "data/ncbi_file.{num}.zip"
    shell:
        "{params.datasets_binary} download genome accession --filename {output.ncbi_file} --include genome,gff3 --inputfile {input.gtdb_list}"

rule list_ncbi_zip_files:
    input:
        ncbi_file = expand("data/ncbi_file.{num}.zip", num=["{:02d}".format(x) for x in range(50)])
    output:
        zipfile_list = "data/all_ncbi_zipfiles.fofn"
    run:
        with open(output.zipfile_list, 'w') as o:
            for num in ["{:02d}".format(x) for x in range(50)]:
                o.write("data/ncbi_file.{num}.zip\n".format(num))

rule unzip_ncbi:
    input:
        zipfile_list = "data/all_ncbi_zipfiles.fofn"
    output:
        gtdb_fasta_fofn = "data/gtdb_fasta.fofn",
        gtdb_gff_fofn = "data/gtdb_gff.fofn"
    shell:
        "cat {input.zipfile_list} | while read line; do unzip $line; done && ls ncbi_file*/data/*/*.fna > {output.gtdb_fasta_fofn} && ls ncbi_file*/data/*/*.gff3 > {output.gtdb_gff_fofn}"





rule download_virosaurus:
    output:
        virosaurus_fasta = "genomes/virosaurus.fasta"
    shell:
        "wget -o {output.virosaurus_fasta} https://viralzone.expasy.org/resources/Virosaurus/2020%5F4/virosaurus98%5Fvertebrate-20200330.fas.gz"

rule create_virus_taxfile:
    input:
        virosaurus_fasta = "genomes/virosaurus"
    output:
        virus_tax = "phylo/virus_taxonomy.tsv.gz"
    script:
        "scripts/create_virosaurus_tax.py"

rule get_eupathdb_genomes:
    params:
        "data/veupathdb_summary.txt"
    output:
        "genomes/eupathdb/fofn"





rule download_and_derep:
    input:
        ar53_tax = "phylo/ar53_taxonomy_{gtdb_release}.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy_{gtdb_release}.tsv.gz"
    output:
        genome_fofn = "data/genomes.fofn"
    script:
        "script/dlderep.py"





