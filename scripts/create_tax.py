def create_tax_ref(nodes, names):
    name_dict = {}
    with open(names) as f:
        for line in f:
            node, p1, name, p2, other, p3, name_type = line.split("\t")[:7]
            if name_type == "scientific name":
                name_dict[node] = name
    graph = {}
    superkingdom = set()
    species_strains = set()
    with open(nodes) as f:
        node, p1, parent, p2, level = line.split("\t")[:5]
        graph[node] = parent
        if level == "superkingdom":
            superkingdom.add(node)
        elif level == "species" or level == "strain":
            species_strains.add(level)
    for i in species_strains:
        alist = [i]
        current = i
        while not current in superkingdom:
            current = graph[current]
            alist.append(current)
        alist.reverse()
        taxstring = name_dict[alist[0]]
        for i in alist[1:]
            taxstring += ';' + name_dict[i]
        tax_ref[i] = taxstring
    return(tax_ref)

def process_virus(viral_fasta, tax_ref, fasta_outdir):
    fasta_dict = {}
    tax_id_dict = {}
    with gzip.open(viral_fasta, 'rt') as f:
        for line in f:
            if line.startswith(">"):
                name = line.split(':')[0][1:]
                for j in line.split('; '):
                    if j.startswith("taxid="):
                        taxid = j.split('=')[0]
                        taxid.replace(';', '')
                tax_id_dict[name] = taxid
                if name in fasta_dict:
                    fasta_dict[name].append("")
                else:
                    fasta_dict[name] = [""]
            elif line.startswith('--'):
                pass
            else:
                fasta_dict[name][-1] += line.rstrip()
    for i in fasta_dict:
        with open(os.path.join(fasta_outdir, i + ".fna", 'w') as o:
            o.write(">")


