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
    output:
        accession_list = "data/gtdb_accessions",
        gtdb_list = expand("data/gtdb.{num}.list", num=["{:02d}".format(x) for x in range(50)])
    shell:
        "zcat {input.ar53_tax} {input.bac120_tax} | awk -F'\t' '{{print substr($1, 4); }}' > {output.accession_list} && "
        "split -n l/50 -d --additional-suffix=.list data/gtdb_accessions data/gtdb."


rule download_gtdb:
    input:
        gtdb_list =  "data/gtdb.{num}.list"
    params:
        datasets_binary = config["datasets_binary"]
    output:
        ncbi_file = "data/ncbi_file.{num}.zip"
    shell:
        "{params.datasets_binary} download genome accession --filename {output.ncbi_file} --include genome,gff3 --inputfile {input.gtdb_list}"


rule unzip_ncbi:
    input:
        gtdb_list =  expand("data/gtdb.{num}.list", num=["{:02d}".format(x) for x in range(50)]),
        ncbi_file = expand("data/ncbi_file.{num}.zip", num=["{:02d}".format(x) for x in range(50)])
    output:
        gtdb_fasta_fofn = "data/gtdb_fasta.list",
        gtdb_gff_fofn = "data/gtdb_gff.list"
    run:
        import glob, sys, os, subprocess
        with open(output.gtdb_fasta_fofn,'w') as fo, open(output.gtdb_gff_fofn,'w') as go:
            for i, j in zip(input.ncbi_file, input.gtdb_list):
                subprocess.Popen("unzip -o {} -d data/".format(i), shell=True).wait()
                with open(j) as f:
                    for line in f:
                        accession = line.rstrip()
                        fastas = glob.glob("data/ncbi_dataset/data/{}/*_genomic.fna".format(accession))
                        if len(fastas) == 1:
                            fasta = fastas[0]
                        else:
                            sys.stderr.write("More than one fasta file found for accession {}.".format(accession))
                            sys.exit()
                        if os.path.exists("data/ncbi_dataset/data/{}/genomic.gff".format(accession)):
                            gff = "data/ncbi_dataset/data/{}/genomic.gff".format(accession)
                        else:
                            gff = None
                        shell("gzip {}".format(fasta))
                        fo.write("{}\t{}.gz\n".format(accession, fasta))
                        if not gff is None:
                            shell("gzip {}".format(gff))
                            go.write("{}\t{}.gz\n".format(accession, gff))

rule deduplicate:
    input:
        gtdb_fasta_fofn = "data/gtdb_fasta.list",
        gtdb_gff_fofn = "data/gtdb_gff.list",
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz"
    output:
        deduplicate_list = "data/gtdb_deduplicated_fastas.list"





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
        "genomes/eupathdb.fofn"





rule download_and_derep:
    input:
        ar53_tax = "phylo/ar53_taxonomy_{gtdb_release}.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy_{gtdb_release}.tsv.gz"
    output:
        genome_fofn = "data/genomes.fofn"
    script:
        "script/dlderep.py"





