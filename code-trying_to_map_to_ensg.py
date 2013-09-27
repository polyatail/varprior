import sys
from collections import defaultdict
import networkx
import tempfile
import subprocess

ENSGENE = "20130918_ensGene.tab"
ENSGTP = "20130918_ensGtp.tab"
ENSEMBL_TO_GENE_NAME = "20130918_ensemblToGeneName.tab"
GENE_ARCHIVE = "gene_archive.txt"
STRING_NETWORK_LINKS = "human-protein.links.v9.05.txt"

GENOTYPE_VCF_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/supporting/AFR.2of4intersection_allele_freq.20100804.genotypes.vcf.gz"

if len(sys.argv) <> 2:
    print "usage: %s gene_name" % sys.argv[0]
    sys.exit(1)

# load ensGene (contains chrom:start-end shit)
enst = None

# load ensGtp (translates ENSG -> ENST -> ENSP)
enst_to_ensp = {}
ensp_to_enst = {}

ensg_to_enst = defaultdict(list)

for line in [x.strip().split() for x in open(ENSGTP)]:
    ensg_to_enst[line[0]].append(line[1])

    if len(line) < 3:
        continue

    assert line[1] not in enst_to_ensp
    enst_to_ensp[line[1]] = line[2]

    assert line[2] not in ensp_to_enst
    ensp_to_enst[line[2]] = line[1]

# load ensemblToGeneName (translates ENST -> short names, e.g. IL2RG)
enst_to_gene_name = {}
gene_name_to_enst = defaultdict(list)

for line in [x.strip().split() for x in open(ENSEMBL_TO_GENE_NAME)]:
    enst_to_gene_name[line[0]] = line[1]
    gene_name_to_enst[line[1]].append(line[0])

# load gene archive, maps retired ENSP -> ENSG
retired_ensp_to_ensg = defaultdict(list)

for line in [x.strip().split() for x in open(GENE_ARCHIVE)]:
    if line[4] == "\\N":
        continue

    retired_ensp_to_ensg[line[4]].append(line[0])

for ensp in retired_ensp_to_ensg.keys():
    retired_ensp_to_ensg[ensp] = list(set(retired_ensp_to_ensg[ensp]))

    for ensg in retired_ensp_to_ensg[ensp]:
        if ensg not in ensg_to_enst:
            retired_ensp_to_ensg[ensp].remove(ensg)

    if len(retired_ensp_to_ensg[ensp]) == 0:
        del retired_ensp_to_ensg[ensp]

# load STRING network
ensp_graph = networkx.Graph()
skipped_nodes = []

for line in [x.strip().split() for x in open(STRING_NETWORK_LINKS)]:
    line[0] = line[0].replace("9606.", "")
    line[1] = line[1].replace("9606.", "")

    if line[0] in skipped_nodes or \
       line[1] in skipped_nodes:
        continue

    if line[0] not in ensp_to_enst:
        if line[0] not in retired_ensp_to_ensg:
            skipped_nodes.append(line[0])
            print "skipping unknown node %s" % line[0]
            continue

        line[0] = retired_ensp_to_ensg[line[0]][0]

    if line[1] not in ensp_to_enst:
        if line[1] not in retired_ensp_to_ensg:
            skipped_nodes.append(line[1])
            print "skipping unknown node %s" % line[1]
            continue

        line[1] = retired_ensp_to_ensg[line[1]][0]

    ensp_graph.add_edge(line[0], line[1], weight=float(line[2]))

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
    enst = gene_name_to_enst[sys.argv[1]]
    pass

    import pdb; pdb.set_trace()
