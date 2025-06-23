import os, gzip, subprocess, shutil

tax_dict = {}
fasta_locations = {}
species_dict = {}
try:
    os.makedirs("temp")
except FileExistsError:
    if len(os.listdir("temp")) > 0:
        sys.exit("Please empty temp directory in working directory.")
try:
    os.makedirs("temp_dedup")
except FileExistsError:
    if len(os.listdir("temp_dedup")) > 0:
        sys.exit("Please empty temp directory in working directory.")
with gzip.open(snakemake.input.ar53_tax, 'rt') as f:
    for line in f:
        accession, tax = line.rstrip().split("\t")
        if not tax in tax_dict:
            tax_dict[tax] = []
        tax_dict[tax].append(accession[3:])
with gzip.open(snakemake.input.bac120_tax, 'rt') as f:
    for line in f:
        accession, tax = line.rstrip().split("\t")
        if not tax in tax_dict:
            tax_dict[tax] = []
        tax_dict[tax].append(accession[3:])
with open(snakemake.input.gtdb_fasta_fofn) as f:
    for line in f:
        accession, location = line.rstrip().split("\t")
        fasta_locations[accession] = location
with open(snakemake.input.euk_tax) as f:
    for line in f:
        accession, tax = line.rstrip().split("\t")
        if not tax in tax_dict:
            tax_dict[tax] = []
        tax_dict[tax].append(accession)
with open(snakemake.input.euk_list) as f:
    for line in f:
        accession, taxid, fasta, gff = line.rstrip().split("\t")
        fasta_locations[accession] = fasta

deduplist = []

for species, accessions in tax_dict.items():
    for i in accessions:
        species_dict[i] = species
    new_accessions = []
    for i in accessions:
        if i in fasta_locations:
            new_accessions.append(i)
    accessions = new_accessions
    if len(accessions) == 0:
        continue
    elif len(accessions) == 1:
        deduplist.append(accessions[0])
    elif len(accessions) > 1:
        toremove = set()
        for i in accessions:
            destination = "temp/{}.fna".format(i)
            shutil.copy2(fasta_locations[i], destination)
            toremove.add(destination)
        subprocess.Popen("{}/{} --distance 0.02 --threads {} temp temp_dedup".format(snakemake.scriptdir, snakemake.params.derep_script, snakemake.threads), shell=True).wait()
        for i in toremove:
            os.remove(i)
        for i in os.listdir("temp_dedup"):
            accession = i[:-4]
            deduplist.append(accession)
            os.remove(os.path.join("temp_dedup", i))
with open(snakemake.output.deduplicate_list, 'w') as o:
    for i in deduplist:
        o.write("{}\t{}\n".format(i, species_dict[i]))