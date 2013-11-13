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

class AnalyzeTrio():
    def __init__(self, sample_names, pedigree, vcf_file, varprior_db, genome_fasta,
                 ensgene_file = None, enst_to_gene_name_file = None, string_alias_file = None):
        self.sample_names = sample_names
        self.pedigree = pedigree
        self.vcf_file = vcf_file
        self.varprior_db = varprior_db
        self.genome_fasta = genome_fasta

        self.variant_matrix = {}

        self.pyf_genome = self.load_genome()
        self.c = self.load_annotation(ensgene_file, enst_to_gene_name_file, string_alias_file)
        self.load_trio_vcf()

    ##
    ## FILE/DATABASE LOADING METHODS
    ##

    def load_genome(self):
        return pyfasta.Fasta(self.genome_fasta, record_class=pyfasta.records.MutNpyFastaRecord)

    def load_annotation(self, ensgene_file, enst_to_gene_name_file, string_alias_file):
        # load ensGene table into an SQLite database creating
        # if it doesn't already exist, return cursor

        if os.path.isfile(ensgene_db):
            conn = sqlite3.connect(ensgene_db)
            c = conn.cursor()
    
            try:
                idx_version = c.execute("SELECT value FROM metadata WHERE key = 'version'").fetchone()[0]
    
                if idx_version != "varprior-1.0":
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
            if ensgene_file == None or \
               enst_to_gene_name_file == None or \
               string_alias_file == None:
                raise ValueError("Cannot create database without input files")

            conn = sqlite3.connect(ensgene_db)
            c = conn.cursor()
    
            c.execute("CREATE TABLE metadata (key text, value text)")
            c.execute("CREATE TABLE ensGene (bin int, name text, chrom text, strand text, txStart int, \
                                             txEnd int, cdsStart int, cdsEnd int, exonCount int, \
                                             exonStarts text, exonEnds text, score int, name2 text, \
                                             cdsStartStat text, cdsEndStat text, exonFrames text)")
            c.execute("CREATE TABLE ensp_to_enst (ensp text, enst text)")   
            c.execute("CREATE TABLE enst_to_gene_name (enst text, gene_name text)")   

            c.execute("INSERT INTO metadata VALUES ('version', 'varprior-1.0')")
            c.execute("INSERT INTO metadata VALUES ('records', '-1')")

            # load ensemblToGeneName (translates ENST -> short names, e.g. IL2RG)
            for line in [x.strip().split() for x in open(enst_to_gene_name_file)]:
                c.execute("INSERT INTO enst_to_gene_name VALUES ('%s', '%s')" % (line[0], line[1]))

            # load string alises (translates ENSP -> ENST/ENSG)
            for line in [x.strip().split() for x in open(string_alias_file)]:
                if line[2].startswith("ENST"):
                    c.execute("INSERT INTO ensp_to_enst VALUES ('%s', '%s')" % (line[1], line[2]))

            # load ensgene
            for line_num, line in enumerate([x.strip().split() for x in open(ensgene_file)]):
                if line_num == 0:
                    ensgene_keys = line
                else:
                    c.execute("INSERT INTO ensGene VALUES (%s)" % ", ".join(["'%s'" % x for x in line]))
    
            c.execute("UPDATE metadata SET value = '%s' WHERE key = 'records'" % line_num)
    
            conn.commit()
    
        return c

    def load_trio_vcf(self):
        """
        Loads a given VCF file into a sparse matrix

        Input: VCF (vcf_file) containing samples of interest (sample_names)
        Output: Sparse arrays (variant_matrix) of every chromosome for every
                sample, containing variants reported in the VCF file
        """

        # generate empty sparse arrays for each sample    
        for sample in self.sample_names:
            self.variant_matrix[i] = {}
    
            for chrom in self.pyf_genome.keys():
                matrix[sample][chrom] = blist([0])
                matrix[sample][chrom] *= len(self.pyf_genome[chrom])

        for line in open(self.vcf_file, "r"):
            # parse the header, find columns of interest
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header = line.strip().split()

                    try:
                        kept_cols = [(x, header.index(x)) for x in self.sample_names]
                    except ValueError, error:
                        raise ValueError("Specified sample (%s) not in VCF (%s)" % (error.split()[0], self.vcf_file))
 
                continue
    
            line_split = line.strip().split()
            chrom = "chr%s" % line_split[0]
    
            for sample, column in kept_cols: 
                pos = int(line_split[1])
                ref = line_split[3]
                alt = [(str(x + 1), y) for x, y in enumerate(line_split[4].split(","))]
    
                code = dict([(str(0), ref)] + alt)
    
                genotype = line_split[column][:3]
                genotype_coded = reduce(lambda x, y: x.replace(y, code[y]), code, genotype)
    
                self.variant_matrix[sample][chrom][pos] = genotype_coded

    ##
    ## INHERITANCE FILTER METHODS
    ##

    def mendelian_filter(self, m_filter, chrom, start, end):
        """
        matrix[0,1,2] = sparsearray(0:24, 0:N, dtype="string")
        pattern in ("recessive", "compound_het", "compound_het_denovo", "denovo_dominant")
        """

        hits = []

        for pos in range(start, end):
            mother, father, child = [self.variant_matrix[self.pedigree[x]][chrom][pos] for x in ("mother", "father", "child")]

            # if no variant reported, continue
            if mother == 0 or \
               father == 0 or \
               child == 0:
                continue

            mother, father, child = [x.split("/") for x in ("mother", "father", "child")]

            if m_filter(mother, father, child):
                hits.append((chrom, pos + start))

    @staticmethod
    def _mf_recessive(mother, father, child) 
        # homozygous in child, heterozygous in both parents
        if len(set(child)) == 1 and \
           len(set(mother)) == 2 and \
           len(set(father)) == 2:
            return True

    @staticmethod
    def _mf_denovo_dominant(mother, father, child):
        # one or both alleles in child must not be in either parent
        if (child[0] not in mother and \
            child[0] not in father) or \
           (child[1] not in mother and \
            child[1] not in father):
            return True

    @staticmethod
    def _mf_compound_het_denovo(mother1, father1, child1, mother2, father2, child2):
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

        # 1. child heterozygous at both loci
        if child1[0] == child1[1] or \
           child2[0] == child2[1]:
            return False
    
        # there are four hypotheses. each position has two alleles,
        # and any combination of the two alleles at each position
        # could be causative
        for child1_a, child2_a in itertools.product(child1, child2):
            # 2. parents can't be homozygous for either of these variants
            if (mother1[0] == child1_a and \
                mother1[1] == child1_a) or \
               (father1[0] == child1_a and \
                father1[1] == child1_a):
                continue
    
            if (mother2[0] == child2_a and \
                mother2[1] == child2_a) or \
               (father2[0] == child2_a and \
                father2[1] == child2_a):
                continue
    
            # 3. causative variant must be heterozygous in one or neither parent
            if (child1_a in mother1 and \
                child1_a in father1) or \
               (child2_a in mother2 and \
                child2_a in father2):
                continue
    
            # 5. both variants can't come from the same parent
            if (child1_a in mother1 and \
                child2_a in mother2) or \
               (child1_a in father1 and \
                child2_a in father2):
                continue

            return True

    ##
    ## GENOME/CODING POTENTIAL METHODS
    ##

    def non_synonymous(self, gene_name, snp_list):
        # given gene name and set of SNPs, does protein product potentially
        # yield a different protein product from that expected?
        # snp_list = ((pos, mut), (pos, mut), ...)
    
        # gene_name -> enst
        all_enst = self.c.execute("SELECT enst FROM enst_to_gene_name WHERE gene_name = '%s'" % gene_name).fetchall()

        if len(all_enst) == 0:
            raise ValueError("Gene (%s) not found in database" % gene_name)
    
        # examine each transcript isoform of this gene
        for enst in all_enst:
            # figure out which chromosome gene is on (only do this once)
            try:
                chrom
            except ValueError:
                chrom = self.c.execute("SELECT chrom FROM ensGene WHERE name = '%s'" % enst).fetchone()[0]
    
                chrom_orig = self.pyf_genome[chrom][:]
                chrom_mut = self.pyf_genome[chrom][:]
    
                # make given mutations
                for mut, pos in snp_list:
                    chrom_mut[pos] = mut
    
            # fetch exons
            exonStarts, exonEnds = [map(int, x.split(",")[:-1]) for x in self.c.execute("SELECT exonStarts, exonEnds FROM ensGene WHERE name = '%s'" % enst).fetchone()]
            cdsStart, cdsEnd = self.c.execute("SELECT cdsStart, cdsEnd FROM ensGene WHERE name = '%s'" % enst).fetchone()
    
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
                orig_seq.append(chrom_orig[start:end])
                mut_seq.append(chrom_mut[start:end])
    
            # translate original and mutated proteins
            orig_protein = str(Seq("".join(orig_seq)).translate(table=1))
            mut_protein = str(Seq("".join(mut_seq)).translate(table=1))
    
            if orig_protein != mut_protein:
                print enst
                for pos, (orig, mut) in enumerate(zip(orig_protein, mut_protein)):
                    print "\t", orig, pos, mut

    ##
    ## VCF ALLELE FREQUENCY METHODS
    ##

    @staticmethod
    def tabix_af_in_region(vcf_file, chrom, start, end):
        # retrieve allele frequencies of SNPs in given region
        # of tabix-indexed VCF file
        tabix_out = tempfile.NamedTemporaryFile(delete=False)
    
        tabix = subprocess.Popen(["tabix", "-f", "-h", vcf_file,
                                  "%s:%s-%s" % (chrom, start, end)],
                                 stdout=tabix_out)
    
        while tabix.poll() == None:
            time.sleep(0.5)
    
        af_table = self.af_from_vcf(tabix_out.name)
    
        return af_table

    @staticmethod
    def af_from_vcf(vcf_file):
        # allele frequencies from any given vcf
        vcftools_out = tempfile.NamedTemporaryFile()
    
        vcftools = subprocess.Popen(["vcftools", "--vcf", vcf_in,
                                     "--freq", "--out", vcftools_out.name])
    
        while vcftools.poll() == None:
            time.sleep(0.5)
    
        af_table = {}
    
        with open("%s.frq" % vcftools_out.name) as fp:
            fp.next()
    
            for line in fp:
                line_split = line.strip().split()
    
                tmp_af = {}
    
                for allele, freq in [x.split(":") for x in line_split[4:]]:
                    tmp_af[allele] = float(freq)
    
                af_table[int(line_split[1])] = tmp_af
    
        return af_table

    ##
    ## NETWORK METHODS
    ##

    def load_string_network(self):
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

    def find_connected_genes(node_list, max_depth, graph):
        # find genes within n connections of given gene
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

    def plot_network(self):
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

    ##
    ## STATISTICS METHODS
    ##

    @staticmethod
    def newell_ikeda(k, pois_lambda, T, w):
        # newell-ikeda poisson distributed scan statistic
        return 1 - numpy.exp(-pois_lambda ** k * w ** (k - 1) * T / scipy.misc.factorial(k - 1, exact=True))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: %s score_cutoff_float gene_name1 [gene_name2 gene_name3]" % sys.argv[0]
        sys.exit(1)

    import pdb; pdb.set_trace()
