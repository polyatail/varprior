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

if len(sys.argv) <> 2:
    print "usage: %s gene_name" % sys.argv[0]
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

for ensp in ensp_graph.nodes():
    gene_names = set([enst_to_gene_name[x] for x in ensp_to_enst[ensp] if x in enst_to_gene_name])
    gene_names = list(gene_names)

    if len(gene_names) == 0:
        gene_names = [ensp]

    ensp_to_gene_name[ensp] = gene_names

    for name in gene_names:
        gene_name_to_ensp[name].append(ensp)

# generate STRING graph with gene names
gene_graph = networkx.Graph()
gene_graph.add_nodes(gene_name_to_ensp.keys())

for line in [x.strip().split() for x in open(STRING_NETWORK_LINKS)]:
    gene_graph.add_edge(ensp_to_gene_name[line[0][5:]], ensp_to_gene_name[line[1][5:]], weight=float(line[2]))

# retrieve allele frequencies of SNPs in given region from 1000 genomes
def fetch_af_in_region(chrom, start, end):
    tabix_out = tempfile.NamedTemporaryFile()

    tabix = subprocess.Popen(["tabix", "-f", "-h", GENOTYPE_VCF_URL,
                              "%s:%s-%s" % (chrom, start, end)],
                             stdout=tabix_out.name)

    return af_from_vcf(tabix_out.name)

# allele frequencies from any given vcf
def af_from_vcf(vcf_in):
    vcftools_out = tempfile.NamedTemporaryFile()

    vcftools = subprocess.Popen(["vcftools", "--vcf", vcf_in,
                                 "--freq"],
                                stdout=vcftools_out.file)

    return vcftools_out.readlines()

# find genes within n connections of given gene
def find_connected_genes(ensp_list, max_depth):
    assert isinstance(ensp_list, list)

    connected_nodes = []

    # depth-limited search
    def _dls(node, depth):
        if depth == 0:
            return 

        for neighbor_node in ensp_graph.neighbors(node):
            connected_nodes.append(neighbor_node)
            _dls(neighbor_node, depth - 1)
            
    for gene in ensp_list:
        _dls(gene, max_depth)

    return set(connected_nodes)

if __name__ == "__main__":
    import pdb; pdb.set_trace()
