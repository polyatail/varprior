from Bio import Seq
import os
import sqlite3
import numpy
import scipy
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

GENOTYPE_VCF = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/supporting/AFR.2of4intersection_allele_freq.20100804.genotypes.vcf.gz"

# retrieve allele frequencies of SNPs in given region from 1000 genomes
def fetch_af_in_region(chrom, start, end):
    tabix_out = tempfile.NamedTemporaryFile(delete=False)

    tabix = subprocess.Popen(["tabix", "-f", "-h", GENOTYPE_VCF,
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
def load_snps_to_matrix(vcf_in, cols, pyf_genome):
    matrix = {}

    for index in range(len(cols)):
        matrix[index] = {}

        for chrom in pyf_genome.keys():
            matrix[index][chrom] = blist([0])
            matrix[index][chrom] *= len(pyf_genome[chrom])

    for line in open(vcf_in):
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                header = line.strip().split()
                kept_cols = [header.index(x) for x in cols]

            continue

        line_split = line.strip().split()
        chrom = "chr%s" % line_split[0]

        for matrix_index, col_num in enumerate(kept_cols):
            pos = int(line_split[1])
            ref = line_split[3]
            alt = [(str(x + 1), y) for x, y in enumerate(line_split[4].split(","))]

            code = dict([(str(0), ref)] + alt)

            sample = line_split[col_num][:3]
            sample_coded = reduce(lambda x, y: x.replace(y, code[y]), code, sample)

            matrix[matrix_index][chrom][pos] = sample_coded

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

    hits = []
 
    if pattern == "recessive":
        for pos in mother_pos:
            # split genotypes into list
            mother_genotype = mother[0, pos].split("/")
            father_genotype = father[0, pos].split("/")
            child_genotype = child[0, pos].split("/")

            # homozygous in child, heterozygous in both parents
            if len(set(child_genotype)) == 1 and \
               len(set(mother_genotype)) == 2 and \
               len(set(father_genotype)) == 2:
                hits.append((chrom, pos + start))

    if pattern == "denovo_dominant":
        for pos in mother_pos:
            # split genotypes into list
            mother_genotype = mother[0, pos].split("/")
            father_genotype = father[0, pos].split("/")
            child_genotype = child[0, pos].split("/")

            # one or both alleles in child must not be in either parent
            if (child_genotype[0] not in mother_genotype and \
                child_genotype[0] not in father_genotype) or \
               (child_genotype[1] not in mother_genotype and \
                child_genotype[1] not in father_genotype):
                hits.append((chrom, pos + start))

    if pattern == "compound_het_denovo":
        """
        From http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0070151
        
        1. A variant has to be in a heterozygous state in all affected individuals.
        2. A variant must not occur in a homozygous state in any of the unaffected individuals.
        3. A variant that is heterozygous in an affected child must be heterozygous in exactly one of the parents.
        This rule is a compact version of:
        
        3a. The variant must not be heterozygous in both parents.
        3b. The variant must be present in at least one of the parents.
        4. A gene must have two or more heterozygous variants in each of the affected individuals.
        5. There must be at least one variant transmitted from the paternal side and one transmitted from the maternal side.
        """

        # test all combinations of positions in the given region
        for pos1, pos2 in itertools.combinations(mother_pos, 2):
            # split genotypes into list
            mother_genotype1 = mother[0, pos1].split("/")
            father_genotype1 = father[0, pos1].split("/")
            child_genotype1 = child[0, pos1].split("/")

            mother_genotype2 = mother[0, pos2].split("/")
            father_genotype2 = father[0, pos2].split("/")
            child_genotype2 = child[0, pos2].split("/")

            # 1. child heterozygous at both loci
            if child_genotype1[0] == child_genotype1[1] or \
               child_genotype2[0] == child_genotype2[1]:
                continue

            # there are four hypotheses. each position has two alleles,
            # and any combination of the two alleles at each position
            # could be causative
            for genotype1, genotype2 in itertools.product(child_genotype1, child_genotype2):
                # 2. parents can't be homozygous for either of these variants
                if (mother_genotype1[0] == genotype1 and \
                    mother_genotype1[1] == genotype1) or \
                   (father_genotype1[0] == genotype1 and \
                    father_genotype1[1] == genotype1):
                    continue

                if (mother_genotype2[0] == genotype2 and \
                    mother_genotype2[1] == genotype2) or \
                   (father_genotype2[0] == genotype2 and \
                    father_genotype2[1] == genotype2):
                    continue

                # 3. causative variant must be heterozygous in one or neither parent
                if (genotype1 in mother_genotype1 and \
                    genotype1 in father_genotype1) or \
                   (genotype2 in mother_genotype2 and \
                    genotype2 in father_genotype2):
                    continue

                # 5. both variants can't come from the same parent
                if (genotype1 in mother_genotype1 and \
                    genotype2 in mother_genotype2) or \
                   (genotype1 in father_genotype1 and \
                    genotype2 in father_genotype2):
                    continue

                hits.append((chrom, pos1 + start, pos2 + start))

    return hits

# given gene name and set of SNPs, does protein product potentially
# yield a different protein product from that expected?
def non_synonymous(gene_name, snp_list, ensgene, hg19):
    # snp_list = ((pos, mut), (pos, mut), ...)

    # gene_name -> ensp -> enst


    # examine each transcript isoform of this gene
    for enst in all_enst:
        # figure out which chromosome gene is on (only do this once)
        try:
            chrom
        except ValueError:
            chrom = ensgene.execute("SELECT chrom FROM ensGene WHERE name = '%s'" % enst).fetchone()[0]

            hg19_chrom_orig = hg19[chrom][:]
            hg19_chrom_mut = hg19[chrom][:]

            # make given mutations
            for mut, pos in snp_list:
                hg19_chrom_mut[pos] = mut

        # fetch exons
        exonStarts, exonEnds = [map(int, x.split(",")[:-1]) for x in ensgene.execute("SELECT exonStarts, exonEnds FROM ensGene WHERE name = '%s'" % enst).fetchone()]
        cdsStart, cdsEnd = ensgene.execute("SELECT cdsStart, cdsEnd FROM ensGene WHERE name = '%s'" % enst).fetchone()

        # non-coding
        if cdsStart == cdsEnd:
            continue

        coding_exons = []

        # reduce coordinates to those that are coding
        for start, end in zip(exonStarts, exonEnds):
            # exon occurs before coding start
            if end < cdsStart:
                continue
            # edge case: exon contains coding start AND end
            elif start <= cdsStart < end and \
                 start <= cdsEnd < end:
                coding_exons.append((cdsStart, cdsEnd))
                break
            # exon contains coding start
            elif start <= cdsStart < end:
                coding_exons.append((cdsStart, end))
            # exon contains coding end
            elif start <= cdsEnd < end:
                coding_exons.append((start, cdsEnd))
                break

            # everything else should be in the middle of the coding region
            coding_exons.append((start, end))

        # any of provided snps occur here? if not, return false
        snps_in_coding = []

        for pos, mut in snp_list:
            for start, end in coding_exons:
                if start <= pos <= end:
                    snps_in_coding.append((pos, mut))
                    break

        if len(snps_in_coding) == 0:
            return False

        # pull nucleotide sequence of regions
        orig_seq = []
        mut_seq = []

        for start, end in coding_exons:
            orig_seq.append(hg19_chrom_orig[start:end])
            mut_seq.append(hg19_chrom_mut[start:end])

        # translate original and mutated proteins
        orig_protein = str(Seq("".join(orig_seq)).translate(table=1))
        mut_protein = str(Seq("".join(mut_seq)).translate(table=1))

        if orig_protein != mut_protein:
            print enst
            for pos, (orig, mut) in enumerate(zip(orig_protein, mut_protein)):
                print "\t", orig, pos, mut

# newell-ikeda poisson distributed scan statistic
def newell_ikeda(k, pois_lambda, T, w):
    return 1 - numpy.exp(-pois_lambda ** k * w ** (k - 1) * T / scipy.misc.factorial(k - 1, exact=True))

# load ensGene table into an SQLite database creating
# if it doesn't already exist, return cursor
def load_ensgene(ensgene_file, ensgene_db):
    if os.path.isfile(ensgene_db):
        conn = sqlite3.connect(ensgene_db)
        c = conn.cursor()

        try:
            idx_version = c.execute("SELECT value FROM metadata WHERE key = 'version'").fetchone()[0]

            if idx_version != "varprior_ensgene-1.0":
                raise ValueError("ensGene version (%s) incompatible with this version of Varprior" % idx_version)

            record_count = int(c.execute("SELECT value FROM metadata WHERE key = 'records'").fetchone()[0])

            if record_count == "-1":
                raise ValueError("Unfinished/partial database provided")

            records_found = int(c.execute("SELECT COUNT(*) FROM ensGene").fetchone()[0])

            if record_count <> records_found:
                raise ValueError("Corrupt index. Expected %s records, found %s" % (record_count, records_found))
        except (OperationalError, DatabaseError), error:
            raise ValueError("Problem with SQLite database: %s" % error)
    else:
        conn = sqlite3.connect(ensgene_db)
        c = conn.cursor()

        c.execute("CREATE TABLE metadata (key text, value text)")
        c.execute("CREATE TABLE ensGene (bin int, name text, chrom text, strand text, txStart int, \
                                         txEnd int, cdsStart int, cdsEnd int, exonCount int, \
                                         exonStarts text, exonEnds text, score int, name2 text, \
                                         cdsStartStat text, cdsEndStat text, exonFrames text)")

        c.execute("INSERT INTO metadata VALUES ('version', 'varprior_ensgene-1.0')")
        c.execute("INSERT INTO metadata VALUES ('records', '-1')")
    
        for line_num, line in enumerate([x.strip().split() for x in open(ensgene_file)]):
            if line_num == 0:
                ensgene_keys = line
            else:
                c.execute("INSERT INTO ensGene VALUES (%s)" % ", ".join(["'%s'" % x for x in line]))

        c.execute("UPDATE metadata SET value = '%s' WHERE key = 'records'" % line_num)

        conn.commit()

    return c

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

