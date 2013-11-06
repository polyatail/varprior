from scipy.sparse import dok_matrix
import time
import itertools
import sys
from collections import defaultdict
import networkx
import tempfile
import subprocess

ENSGENE = "20130918_ensGene.tab"
ENSEMBL_TO_GENE_NAME = "20130918_ensemblToGeneName.tab"
STRING_NETWORK_LINKS = "human-protein.links.v9.05.txt"
STRING_NETWORK_ALIAS = "human-protein.aliases.v9.05.txt"

GENOTYPE_VCF_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/supporting/AFR.2of4intersection_allele_freq.20100804.genotypes.vcf.gz"

# retrieve allele frequencies of SNPs in given region from 1000 genomes
def fetch_af_in_region(chrom, start, end):
    tabix_out = tempfile.NamedTemporaryFile(delete=False)

    tabix = subprocess.Popen(["tabix", "-f", "-h", GENOTYPE_VCF_URL,
                              "%s:%s-%s" % (chrom, start, end)],
                             stdout=tabix_out)

    while tabix.poll() == None:
        time.sleep(0.5)

    af_table = af_from_vcf(tabix_out.name)

    return af_table

# allele frequencies from any given vcf
def af_from_vcf(vcf_in):
    vcftools_out = tempfile.NamedTemporaryFile()

    vcftools = subprocess.Popen(["vcftools", "--vcf", vcf_in,
                                 "--freq", "--out", vcftools_out.name])

    while vcftools.poll() == None:
        time.sleep(0.5)

    af = {}

    with open("%s.frq" % vcftools_out.name) as fp:
        fp.next()

        for line in fp:
            line_split = line.strip().split()

            tmp_af = {}

            for allele, freq in [x.split(":") for x in line_split[4:]]:
                tmp_af[allele] = float(freq)

            af[int(line_split[1])] = tmp_af

    return af

# find genes within n connections of given gene
def find_connected_genes(node_list, max_depth, graph):
    assert isinstance(node_list, list)

    connected_nodes = []

    # depth-limited search
    def _dls(node, depth):
        if depth == 0:
            return 

        for neighbor_node in graph.neighbors(node):
            connected_nodes.append(neighbor_node)
            _dls(neighbor_node, depth - 1)
            
    for gene in gene_list:
        _dls(gene, max_depth)

    return set(connected_nodes)

# load specified cols from VCF into sparse matrix
def load_snps_to_matrix(vcf_in, cols):
    matrix = {}

    for index in range(len(cols)):
        matrix[index] = dok_matrix((25, 250000000), dtype="string")

    for line in open(vcf_in):
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                header = line.strip().split()
                kept_cols = [header.index(x) for x in cols]

            continue

        line_split = line.strip().split()

        if line_split[0] == "X":
            chrom = 23
        elif line_split[0] == "Y":
            chrom = 24
        else:
            chrom = int(line_split[0]) 

        for matrix_index, col_num in enumerate(kept_cols):
            pos = int(line_split[1])
            ref = line_split[3]
            alt = [(str(x + 1), y) for x, y in enumerate(line_split[4].split(","))]

            code = dict([(str(0), ref)] + alt)

            sample = line_split[col_num][:3]
            sample_coded = reduce(lambda x, y: x.replace(y, code[y]), code, sample)

            matrix[matrix_index][chrom, pos] = sample_coded

    return matrix

# mendelian inheritance filter on sparse matrix
def mendelian_filter(matrix, pedigree, pattern, chrom, start, end):
    """
    matrix[0,1,2] = sparsearray(0:24, 0:N, dtype="string")
    pattern in ("recessive", "compound_het", "compound_het_denovo", "denovo_dominant")
    """
 
    if pedigree == "default":
        pedigree = {"child": 0, "mother": 1, "father": 2}
 
    mother, father, child = [matrix[pedigree[x]][chrom,start:end] for x in ("mother", "father", "child")]
 
    # all must have same positions
    mother_pos = [pos for (_, pos), _ in mother.iteritems()]
    father_pos = [pos for (_, pos), _ in father.iteritems()]
    child_pos = [pos for (_, pos), _ in child.iteritems()]
 
    assert len(set(mother_pos).difference(father_pos)) == 0
    assert len(set(father_pos).difference(child_pos)) == 0
 
    if pattern == "recessive":
        for pos in mother_pos:
            mother_genotype = mother[0, pos].split("/")
            father_genotype = father[0, pos].split("/")
            child_genotype = child[0, pos].split("/")

            # homozygous in child, heterozygous in both parents
            if len(set(child_genotype)) == 1 and \
               len(set(mother_genotype)) == 2 and \
               len(set(father_genotype)) == 2:
                print chrom, pos + start, child[0, pos], mother[0, pos], father[0, pos]

    if pattern == "denovo_dominant":
        for pos in mother_pos:
            mother_genotype = mother[0, pos].split("/")
            father_genotype = father[0, pos].split("/")
            child_genotype = child[0, pos].split("/")

            # homozygous in both parents
            if len(set(mother_genotype)) == 2 and \
               len(set(father_genotype)) == 2:
                # at least one allele in child must not be in either parent
                if len(set(child_genotype).difference(mother_genotype)) > 0 and \
                   len(set(child_genotype).difference(father_genotype)) > 0:
                    print chrom, pos + start, child[0, pos], mother[0, pos], father[0, pos]

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: %s score_cutoff_float gene_name1 [gene_name2 gene_name3]" % sys.argv[0]
        sys.exit(1)

    # load ensGene (contains chrom:start-end shit)
    enst = {}
    
    for line_num, line in enumerate([x.strip().split() for x in open(ENSGENE)]):
        if line_num == 0:
            ensgene_keys = line[2:]
        else:
            enst[line[1]] = dict(zip(ensgene_keys, line[2:]))
    
    # load ensemblToGeneName (translates ENST -> short names, e.g. IL2RG)
    enst_to_gene_name = {}
    gene_name_to_enst = defaultdict(list)
    
    for line in [x.strip().split() for x in open(ENSEMBL_TO_GENE_NAME)]:
        enst_to_gene_name[line[0]] = line[1]
        gene_name_to_enst[line[1]].append(line[0])
    
    # load string alises (translates ENSP -> ENST/ENSG)
    ensp_to_enst = defaultdict(list)
    enst_to_ensp = {}
    
    for line in [x.strip().split() for x in open(STRING_NETWORK_ALIAS)]:
        if line[2].startswith("ENST"):
            ensp_to_enst[line[1]].append(line[2])
    
            assert line[2] not in enst_to_ensp
            enst_to_ensp[line[2]] = line[1]
    
    # using above data, try to map ENSP in network to gene names and back
    ensp_to_gene_name = {}
    gene_name_to_ensp = defaultdict(list)
    
    for ensp in ensp_to_enst.keys():
        gene_names = set([enst_to_gene_name[x] for x in ensp_to_enst[ensp] if x in enst_to_gene_name])
        gene_names = list(gene_names)
    
        if len(gene_names) == 0:
            gene_names = [ensp]
    
        ensp_to_gene_name[ensp] = gene_names
    
        for name in gene_names:
            gene_name_to_ensp[name].append(ensp)
    
    # generate STRING graph with gene names
    gene_graph = networkx.Graph(directed=False)
    gene_graph.add_nodes_from(gene_name_to_ensp.keys())
    
    for line in [x.strip().split() for x in open(STRING_NETWORK_LINKS)]:
        if float(line[2]) < float(sys.argv[1]):
            continue
    
        gene1 = ensp_to_gene_name[line[0][5:]]
        gene2 = ensp_to_gene_name[line[1][5:]]
    
        for node1, node2 in itertools.product(gene1, gene2):
            gene_graph.add_edge(node1, node2, weight=float(line[2]))

    # actually do stuff
    import matplotlib.pyplot as plt

    neighbors = list(itertools.chain(*[gene_graph.neighbors(x) for x in sys.argv[2:]]))
    subgraph = gene_graph.subgraph(neighbors + sys.argv[2:])
    pos = networkx.graphviz_layout(subgraph, prog="neato")
    labels = dict(zip(neighbors + sys.argv[2:], neighbors + sys.argv[2:]))

    networkx.draw_networkx_nodes(subgraph, pos, nodelist=neighbors, node_color="#fdb462")
    networkx.draw_networkx_nodes(subgraph, pos, nodelist=sys.argv[2:], node_color="#80b1d3")

    networkx.draw_networkx_edges(subgraph, pos, alpha=0.1)

    networkx.draw_networkx_labels(subgraph, pos, labels)

    plt.figure(1, figsize=(16, 16))
    plt.axis("off")
    plt.show()

    import pdb; pdb.set_trace()

