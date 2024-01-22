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

rule download_gtdb:
    input:
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz"
    params:
        datasets_binary = config["datasets_binary"]
    output:
        accession_list = "data/gtdb_accessions",
        ncbi_file = "data/ncbi_dataset.zip"
    shell:
        "zcat {input.ar53_tax} {input.bac120_tax} | awk -F'\t' '{{print substr($1, 4); }}' > {output.accession_list} && "
        "{params.datasets_binary} download genome accession --filename {output.accession_list} --include genome,gff3 --inputfile {output.accession_list}"


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





