configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]
gtdb_release: config["gtdb_release"]
eupath_summary_file = workflow.source_path("data/veupathdb_summary.txt")

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
                        if len(fastas) == 0:
                            continue
                        elif len(fastas) == 1:
                            fasta = fastas[0]
                        else:
                            sys.stderr.write("More than one fasta file found for accession {}.".format(accession))
                            sys.exit()
                        if os.path.exists("data/ncbi_dataset/data/{}/genomic.gff".format(accession)):
                            gff = "data/ncbi_dataset/data/{}/genomic.gff".format(accession)
                        else:
                            gff = None
                        fo.write("{}\t{}\n".format(accession, fasta))
                        if not gff is None:
                            go.write("{}\t{}\n".format(accession, gff))


rule download_eupathdb:
    params:
         eupath_file = eupath_summary_file
    output:
         eupath_list = "data/eupath.list"
    run:
        import os
        with open(params.eupath_file) as f, open(output.eupath_list, 'w') as o:
            f.readline()
            for line in f:
                splitline = line.rstrip().split("\t")
                gff, accession, taxid, fasta = splitline[0], splitline[7], splitline[14], splitline[17]
                if not accession.startswith(("GCF_","GCA_")):
                    accession = gff.split('/')[-1][:-4]
                shell("mkdir -p data/eupath_gffs && mkdir -p data/eupath_fastas")
                try:
                    gff_out = "data/eupath_gffs/{}.gff".format(accession)
                    if not os.path.exists(gff_out):
                        shell("wget -O {} {}".format(gff_out, gff))
                except:
                    gff_out = "none"
                try:
                    fasta_out = "data/eupath_fastas/{}.fna".format(accession)
                    if not os.path.exists(fasta_out):
                        shell("wget -O {} {}".format(fasta_out, fasta))
                except:
                    fasta_out = "none"
                o.write("{}\t{}\t{}\t{}\n".format(accession, taxid, fasta_out, gff_out))


rule download_virosaurus:
    output:
        virosaurus_fasta = "genomes/virosaurus.fasta.gz"
    shell:
        "wget -O {output.virosaurus_fasta} https://viralzone.expasy.org/resources/Virosaurus/2020%5F4/virosaurus98%5Fvertebrate-20200330.fas.gz"

rule download_ncbi_tax:
    output:
        nodes = "phylo/nodes.dmp",
        names = "phylo/names.dmp"
    shell:
        "wget -O phylo/taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && "
        "tar -C phylo/ -xzvf phylo/taxdump.tar.gz"

rule create_virus_taxfile:
    input:
        virosaurus_fasta = "genomes/virosaurus.fasta.gz",
        nodes = "phylo/nodes.dmp",
        names = "phylo/names.dmp"
    params:
        dataset = "virosaurus",
        virus_dir = "data/virus_genomes"
    output:
        virus_tax = "phylo/virus_taxonomy.tsv"
    script:
        "scripts/create_tax.py"

rule create_euk_taxfile:
    input:
        euk_list = "data/eupath.list",
        nodes = "phylo/nodes.dmp",
        names = "phylo/names.dmp"
    params:
        dataset = "euk"
    output:
        euk_tax = "phylo/euk_taxonomy.tsv"
    script:
        "scripts/create_tax.py"


rule deduplicate:
    input:
        gtdb_fasta_fofn = "data/gtdb_fasta.list",
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz",
        euk_tax = "phylo/euk_taxonomy.tsv",
        euk_list= "data/eupath.list"
    threads:
        64
    output:
        deduplicate_list = "data/gtdb_deduplicated_fastas.list"
    script:
        "scripts/dedup.py"


rule create_index_fastas:
    input:
        deduplicate_list = "data/gtdb_deduplicated_fastas.list",
        gtdb_fasta_fofn = "data/gtdb_fasta.list",
        euk_list= "data/eupath.list",
        virus_tax = "phylo/virus_taxonomy.tsv"
    output:
        taxonomy_file = "final_database/taxonomy.tsv",
        index_fofn = "final_database/index_fofn"
    script:
        "scripts/index_fasta.py"

rule index_fasta:
    input:
        index_fofn = "db/final_database/index_fofn"
    output:
        index_out_fofn = "db/final_database/index_out_fofn"
    threads:
        64
    shell:
        "cat {input.index_fofn} | while read line; do echo bwa index -p ${line:0:-3} $file; done > parallel -j {threads} && "
        "cat {input.index_fofn} | while read line; do echo ${line:0:-3}; done > {output.index_out_fofn}"

