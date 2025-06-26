configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]
gtdb_release = config["gtdb_release"]
gtdb_subrelease = config["gtdb_subrelease"]
eupath_summary_file = workflow.source_path("data/veupathdb_summary.txt")


onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

rule get_gtdb_taxfiles:
    params:
        gtdb_release = gtdb_release,
        gtdb_subrelease = gtdb_subrelease
    output:
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz"
    shell:
        "wget -O {output.ar53_tax} https://data.gtdb.ecogenomic.org/releases/release{params.gtdb_release}/{params.gtdb_release}.{params.gtdb_subrelease}/ar53_taxonomy_r{params.gtdb_release}.tsv.gz --no-check-certificate && " 
        "wget -O {output.bac120_tax} https://data.gtdb.ecogenomic.org/releases/release{params.gtdb_release}/{params.gtdb_release}.{params.gtdb_subrelease}/bac120_taxonomy_r{params.gtdb_release}.tsv.gz --no-check-certificate"

rule split_gtdb_accessions:
    input:
        ar53_tax = "phylo/ar53_taxonomy.tsv.gz",
        bac120_tax = "phylo/bac120_taxonomy.tsv.gz"
    output:
        accession_list = "data/gtdb_accessions",
        gtdb_list = expand("data/gtdb.{num}.list", num=["{:03d}".format(x) for x in range(200)])
    shell:
        "zcat {input.ar53_tax} {input.bac120_tax} | awk -F'\t' '{{print substr($1, 4); }}' > {output.accession_list} && "
        "split -n l/200 -d --additional-suffix=.list data/gtdb_accessions data/gtdb."


rule download_gtdb:
    resources:
        api_calls = 1
    input:
        gtdb_list =  "data/gtdb.{num}.list"
    params:
        datasets_binary = config["datasets_binary"]
    output:
        ncbi_file = "data/ncbi_file.{num}.zip"
    shell:
        "{params.datasets_binary} download genome accession --dehydrated --filename {output.ncbi_file} --include genome,gff3 --inputfile {input.gtdb_list}"

rule unzip_ncbi_and_rehydrate:
    resources:
        api_calls = 1
    input:
        ncbi_file = "data/ncbi_file.{num}.zip"
    output:
        ncbi_dir = directory("data/ncbi_file.{num}")
    shell:
        "unzip -o {input.ncbi_file} -d {output.ncbi_dir} && "
        "datasets rehydrate --directory {output.ncbi_dir}/"

rule collate_ncbi_list:
    input:
        gtdb_list =  expand("data/gtdb.{num}.list", num=["{:03d}".format(x) for x in range(200)]),
        ncbi_file = expand("data/ncbi_file.{num}", num=["{:03d}".format(x) for x in range(200)])
    output:
        gtdb_fasta_fofn = "data/gtdb_fasta.list",
        gtdb_gff_fofn = "data/gtdb_gff.list"
    run:
        import glob, sys, os
        with open(output.gtdb_fasta_fofn,'w') as fo, open(output.gtdb_gff_fofn,'w') as go:
            for i, j in zip(input.ncbi_file, input.gtdb_list):
                with open(j) as f:
                    for line in f:
                        accession = line.rstrip()
                        fastas = glob.glob("{}/ncbi_dataset/data/{}/*_genomic.fna".format(i, accession))
                        if len(fastas) == 0:
                            continue
                        elif len(fastas) == 1:
                            fasta = fastas[0]
                        else:
                            sys.stderr.write("More than one fasta file found for accession {}.".format(accession))
                            sys.exit()
                        if os.path.exists("{}/ncbi_dataset/data/{}/genomic.gff".format(i, accession)):
                            gff = "{}/ncbi_dataset/data/{}/genomic.gff".format(i, accession)
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
                gff, accession, taxid, fasta = splitline[0], splitline[11], splitline[18], splitline[22]
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


rule download_rvdb:
    output:
        rvdb_fasta = "genomes/virus.fasta.gz"
    shell:
        "wget -O {output.rvdb_fasta} https://rvdb.dbi.udel.edu/download/C-RVDBvCurrent.fasta.gz"


rule download_rvdb_prot:
    params:
        rvdb_prot_release = config["rvdb_prot_release"]
    output:
        rvbd_faa = "genomes/virus.faa.xz"
    shell:
        "wget -O {output.rvbd_faa} https://rvdb-prot.pasteur.fr/files/U-RVDBv{params.rvdb_prot_release}-prot.fasta.xz --timeout 2400"


rule download_taxids:
    input:
        rvdb_fasta = "genomes/virus.fasta.gz"
    output:
        taxids = "genomes/virus_taxids.txt",
        virus_accession_list = "genomes/virus_accession.list",
        acc_list_fasta = "acc_list_fasta"
    run:
        import gzip
        shell("zcat {input.rvdb_fasta} | grep '>' | awk -F'|' '{{print $3}}' > {output.virus_accession_list}")
        with open(output.virus_accession_list, 'r') as f:
            acclist = [line.strip() for line in f if line.strip()]
        with open(output.taxids, 'w'):
            for num in range(0, len(acclist), 5000):
                with open(output.acc_list_fasta, 'w') as o:
                    for i in acclist[num:num+5000]:
                        o.write("{}\n".format(i))
                shell("cat {} | epost -db nuccore | esummary | xtract -pattern DocumentSummary -element Caption,TaxId >> {}".format(output.acc_list_fasta, output.taxids))

rule download_taxids_faa:
    input:
        rvdb_fasta = "genomes/virus.faa.xz"
    output:
        taxids = "genomes/virus_prot_taxids.txt",
        virus_prot_accession_list = "genomes/virus_prot_accession.list",
        acc_list_prot = "acc_list_prot"
    run:
        import gzip
        shell("xzcat {input.rvdb_fasta} | grep '>' | awk -F'|' '{{print $3}}' > {output.virus_prot_accession_list}")
        with open(output.virus_prot_accession_list, 'r') as f:
            acclist = [line.strip() for line in f if line.strip()]
        with open(output.taxids, 'w'):
            for num in range(0, len(acclist), 5000):
                with open(output.acc_list_prot, 'w') as o:
                    for i in acclist[num:num+5000]:
                        o.write("{}\n".format(i))
                shell("cat {} | epost -db protein | esummary | xtract -pattern DocumentSummary -element Caption,TaxId >> {}".format(output.acc_list_prot, output.taxids))


rule download_ncbi_tax:
    output:
        nodes = "phylo/nodes.dmp",
        names = "phylo/names.dmp"
    shell:
        "wget -O phylo/taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && "
        "tar -C phylo/ -xzvf phylo/taxdump.tar.gz"

rule create_virus_taxfile:
    input:
        rvdb_fasta = "genomes/virus.fasta.gz",
        rvbd_faa = "genomes/virus.faa.xz",
        virus_accessions = "genomes/virus_accession.list",
        virus_prot_accessions = "genomes/virus_prot_accession.list",
        virus_taxids = "genomes/virus_taxids.txt",
        virus_prot_taxids = "genomes/virus_prot_taxids.txt",
        nodes = "phylo/nodes.dmp",
        names = "phylo/names.dmp"
    params:
        dataset = "virus",
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
    params: 
        derep_script = "dereplicator.py"
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


rule dustmask_db:
    input:
        index_fofn = "final_database/index_fofn",
        db_files = "final_database/cmdd.{num}.bwa.fa.gz"
    output:
        dustmask_db = "dustmasked_db/cmdd.{num}.bwa.fa.gz"
    params:
        dustmasker_binary = config["dustmasker_binary"]
    shell:
        "zcat {input.db_files} | {params.dustmasker_binary} -in - -out - -outfmt fasta -linker 100 -hard_masking -level 16 | gzip > {output.dustmask_db}"


rule build_kraken_db:
    input:
        taxonomy_file = "final_database/taxonomy.tsv",
        index_fofn = "final_database/index_fofn",
        db_files = "dustmasked_db/cmdd.{num}.bwa.fa"
    output:
        kraken_db = "kraken_db/kraken_metadb",
        nuc_a2t = "kraken_db/taxonomy/nucl_gb.accession2taxid",
        prot_a2t = "kraken_db/taxonomy/prot.accession2taxid",
        taxdump = "kraken_db/taxonomy/taxdump.tar.gz"
    params:
        threads = config["kraken_threads"]
    shell:
        "wget -O {output.nuc_a2t}.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz && "
        "gunzip {output.nuc_a2t}.gz && "
        "wget -O {output.prot_a2t}.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz && "
        "gunzip {output.prot_a2t}.gz && "
        "wget -O {output.taxdump}.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && "
        "tar -zxf {output.taxdump}.gz && "
        "gunzip {input.db_files}.gz && "
        "kraken2-build --add-to-library {input.db_files} --db {output.kraken_db} && "
        "kraken2-build --build --threads {params.threads} --db {output.kraken_db}"
