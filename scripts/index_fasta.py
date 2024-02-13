import os

def process_file(f, o, accession):
    for line in f:
        if line.startswith(">"):
            o.write(">" + accession + "|" + line[1:])
        else:
            o.write(line)

def process_list(dedup_file, eupath_file, gtdb_file, viro_file, prefix, taxonomy_out, fasta_out, max_size=1000000000):
    fasta_location_dict = {}
    phylo_dict = {}
    with open(viro_file) as f:
        for line in f:
            accession, fasta, tax = line.rstrip().split("\t")
            fasta_location_dict[accession] = fasta
            phylo_dict[accession] = tax
    with open(gtdb_file) as f:
        for line in f:
            accession, fasta = line.rstrip().split("\t")
            fasta_location_dict[accession] = fasta
    with open(eupath_file) as f:
        for line in f:
            accession, taxid, fasta, gff = line.rstrip().split("\t")
            fasta_location_dict[accession] = fasta
    with open(dedup_file) as f:
        for line in f:
            accession, tax = line.rstrip().split("\t")
            phylo_dict[accession] = tax
    count = 0
    genome_groups = []
    sum_size = 0
    group = []
    with open(taxonomy_out, 'w') as o:
        for i, j in phylo_dict.items():
            sum_size +=  os.path.getsize(fasta_location_dict[i])
            group.append((i, fasta_location_dict[i]))
            o.write(i + "\t" + j + "\n")
            if sum_size >= max_size:
                genome_groups.append(group)
                sum_size = 0
                group = []
    with open(fasta_out, 'w') as out_fofn:
        for num, fastas in enumerate(genome_groups):
            with open("{}.{}.bwa.fa".format(prefix, num), 'w') as o:
                for accession, fasta in fastas:
                    count += 1
                    if fasta.endswith(".gz"):
                        with gzip.open(fasta, 'rt') as f:
                            process_file(f, o, accession)
                    else:
                        with open(fasta) as f:
                            process_file(f, o, accession)
            out_fofn.write("{}.{}.bwa.fa".format(prefix, num))


process_list(snakemake.input.deduplicate_list, snakemake.input.euk_list, snakemake.input.gtdb_fasta_fofn,
             snakemake.input.virus_tax, "final_database/cmdd", snakemake.output.taxonomy_file, snakemake.output.index_fofn)