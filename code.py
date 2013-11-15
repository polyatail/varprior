import copy
from blist import blist
import pyfasta
from Bio.Seq import Seq
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
    def __init__(self, sample_names, pedigree, varprior_db, genome_fasta,
                 vcf_file = None, ensgene_file = None, enst_to_gene_name_file = None, string_alias_file = None):
        self.sample_names = sample_names
        self.pedigree = pedigree
        self.stripped_pedigree = dict([(x, y.replace("-", "")) for x, y in pedigree.items()])
        self.varprior_db = varprior_db
        self.genome_fasta = genome_fasta

        sys.stderr.write("loading genome\n")
        self.pyf_genome = self.load_genome()
        sys.stderr.write("loading annotation\n")
        self.conn, self.c = self.load_annotation(ensgene_file, enst_to_gene_name_file, string_alias_file)
        sys.stderr.write("loading vcf\n")
        self.load_trio_vcf(vcf_file)

    ##
    ## FILE/DATABASE LOADING METHODS
    ##

    def load_genome(self):
        return pyfasta.Fasta(self.genome_fasta, record_class=pyfasta.records.MutNpyFastaRecord)

    def load_annotation(self, ensgene_file, enst_to_gene_name_file, string_alias_file):
        # load ensGene table into an SQLite database creating
        # if it doesn't already exist, return cursor

        if os.path.isfile(self.varprior_db):
            conn = sqlite3.connect(self.varprior_db)
            conn.row_factory = sqlite3.Row
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

            conn = sqlite3.connect(self.varprior_db)
            conn.row_factory = sqlite3.Row
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
            c.execute("INSERT INTO metadata VALUES ('variants', '-1')")

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

        return conn, c

    def load_trio_vcf(self, vcf_file):
        """
        Loads a given VCF file into a sparse matrix

        Input: VCF (vcf_file) containing samples of interest (sample_names)
        Output: Sparse arrays (variant_matrix) of every chromosome for every
                sample, containing variants reported in the VCF file
        """

        variant_count = int(self.c.execute("SELECT value FROM metadata WHERE key = 'variants'").fetchone()[0])

        if variant_count > -1:
            if vcf_file != None:
                db_vcf_file = self.c.execute("SELECT value FROM metadata WHERE key = 'vcf_file'").fetchone()[0]

                if vcf_file != db_vcf_file:
                    raise ValueError("Database variants do not match specified file (%s != %s)" % (db_vcf_file, vcf_file))

            variants_found = int(self.c.execute("SELECT COUNT(*) FROM variants").fetchone()[0])

            if variant_count <> variants_found:
                raise ValueError("Corrupt index. Expected %s records, found %s" % (variant_count, variants_found))
        else:
            if vcf_file == None:
                raise ValueError("Database empty and no vcf_file specified")

            self.c.execute("CREATE TABLE variants (chrom text, pos int, %s)" % ", ".join(["%s text" % x.replace("-", "") for x in self.sample_names]))

            records = 0

            for line in open(vcf_file, "r"):
                # parse the header, find columns of interest
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        header = line.strip().split()
    
                        try:
                            kept_cols = [(x, header.index(x)) for x in self.sample_names]
                        except ValueError, error:
                            raise ValueError("Specified sample (%s) not in VCF (%s)" % (error.split()[0], vcf_file))
     
                    continue
    
                line_split = line.strip().split()
                chrom = "chr%s" % line_split[0]
    
                sample_to_snp = {}
    
                for sample, column in kept_cols: 
                    pos = int(line_split[1])
                    ref = line_split[3]
                    alt = [(str(x + 1), y) for x, y in enumerate(line_split[4].split(","))]
        
                    code = dict([(str(0), ref)] + alt)
        
                    genotype = line_split[column][:3]
                    genotype_coded = reduce(lambda x, y: x.replace(y, code[y]), code, genotype)

                    sample_to_snp[sample] = genotype_coded

                self.c.execute("INSERT INTO variants VALUES ('%s', '%s', %s)" % (chrom, pos, ", ".join(["'%s'" % sample_to_snp[x] for x in self.sample_names])))
                records += 1

            self.c.execute("CREATE INDEX IF NOT EXISTS chrom_index ON variants(chrom);")
            self.c.execute("CREATE INDEX IF NOT EXISTS pos_index ON variants(pos);")

            self.c.execute("INSERT INTO metadata VALUES ('vcf_file', '%s')" % vcf_file)
            self.c.execute("UPDATE metadata SET value = '%s' WHERE key = 'variants'" % records)
            self.conn.commit()

    ##
    ## INHERITANCE FILTER METHODS
    ##

    def mendelian_filter(self, m_filter, chrom, start, end, k = 1):
        """
        matrix[0,1,2] = sparsearray(0:24, 0:N, dtype="string")
        pattern in ("recessive", "compound_het", "compound_het_denovo", "denovo_dominant")
        """

        pos_in_region = self.c.execute("SELECT * FROM variants WHERE chrom = '%s' AND %s <= pos AND pos <= %s" % (chrom, start, end)).fetchall()

        hits = []

        for pos in itertools.combinations(pos_in_region, k):
            mother, father, child = {}, {}, {}

            for i in range(k):
                mother[i], father[i], child[i] = [pos[i][self.stripped_pedigree[x]].split("/") for x in ("mother", "father", "child")]

                # if no variant reported, continue
                if "." in mother[i] or \
                   "." in father[i] or \
                   "." in child[i]:
                    break
            else:
                if m_filter(mother, father, child):
                    hits.append(tuple([chrom] + [start + x["pos"] for x in pos]))

        return hits

    @staticmethod
    def _mf_recessive(mother, father, child):
        # homozygous in child, heterozygous in both parents
        if len(set(child[0])) == 1 and \
           len(set(mother[0])) == 2 and \
           len(set(father[0])) == 2:
            return True

    @staticmethod
    def _mf_denovo_dominant(mother, father, child):
        # one or both alleles in child must not be in either parent
        if (child[0][0] not in mother[0] and \
            child[0][0] not in father[0]) or \
           (child[0][1] not in mother[0] and \
            child[0][1] not in father[0]):
            return True

    @staticmethod
    def _mf_compound_het_denovo(mother, father, child):
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

        mother1, father1, child1 = mother[0], father[0], child[0]
        mother2, father2, child2 = mother[1], father[1], child[1]

        # 1. child heterozygous at both loci
        if child[0][0] == child[0][1] or \
           child[1][0] == child[1][1]:
            return False
    
        # there are four hypotheses. each position has two alleles,
        # and any combination of the two alleles at each position
        # could be causative
        for child1, child2 in itertools.product(child[0], child[1]):
            # 2. parents can't be homozygous for either of these variants
            if (mother[0][0] == child1 and \
                mother[0][1] == child1) or \
               (father[0][0] == child1 and \
                father[0][1] == child1):
                continue
    
            if (mother[1][0] == child2 and \
                mother[1][1] == child2) or \
               (father[1][0] == child2 and \
                father[1][1] == child2):
                continue
    
            # 3. causative variant must be heterozygous in one or neither parent
            if (child1 in mother[0] and \
                child1 in father[0]) or \
               (child2 in mother[1] and \
                child2 in father[1]):
                continue
    
            # 5. both variants can't come from the same parent
            if (child1 in mother[0] and \
                child2 in mother[1]) or \
               (child1 in father[0] and \
                child2 in father[1]):
                continue

            return True

    ##
    ## GENOME/CODING POTENTIAL METHODS
    ##

    def non_synonymous(self, gene_name, sample):
        # given gene name and set of SNPs, does protein product potentially
        # yield a different protein product from that expected?
        # snp_list = ((pos, mut), (pos, mut), ...)

        # sample can't have dashes in it
        stripped_sample = sample.replace("-", "")
   
        # gene_name -> enst
        all_enst = [x[0] for x in self.c.execute("SELECT enst FROM enst_to_gene_name WHERE gene_name = '%s'" % gene_name).fetchall()]

        if len(all_enst) == 0:
            raise ValueError("Gene (%s) not found in database" % gene_name)
   
        # figure out which chromosome gene is on (only do this once)
        chrom, gene_start, gene_end = self.c.execute("SELECT MAX(chrom), MIN(txStart), MAX(txEnd) FROM ensGene WHERE name IN (%s)" % ", ".join(["'%s'" % x for x in all_enst])).fetchone()

        results = {}
 
        # examine each transcript isoform of this gene
        for enst in all_enst:
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
    
            # pull nucleotide sequence of regions
            orig_seq = []

            for start, end in coding_exons:
                orig_seq.append(self.pyf_genome[chrom][start:end])

            mut_seq1, mut_seq2 = [], []
    
            # make given mutations
            for genotype_index, seq_list in ((0, mut_seq1), (1, mut_seq2)):
                pos_in_region = self.c.execute("SELECT * FROM variants WHERE chrom = '%s' AND %s <= pos AND pos <= %s" % (chrom, cdsStart, cdsEnd)).fetchall()

                for row in pos_in_region:
                    self.pyf_genome[chrom][row["pos"]] = row[stripped_sample].split("/")[genotype_index]
 
                for start, end in coding_exons:
                    seq_list.append(self.pyf_genome[chrom].__getitem__(slice(start, end), True))
    
            # translate original and mutated proteins
            orig_protein = str(Seq("".join(orig_seq)).translate(table=1))
            mut_protein1 = str(Seq("".join(mut_seq1)).translate(table=1))
            mut_protein2 = str(Seq("".join(mut_seq2)).translate(table=1))

            allele1 = []

            if orig_protein != mut_protein1:
                for pos, (orig, mut) in enumerate(zip(orig_protein, mut_protein1)):
                    if orig != mut:
                        allele1.append((orig, pos, mut))

            allele2 = []

            if orig_protein != mut_protein2:
                for pos, (orig, mut) in enumerate(zip(orig_protein, mut_protein2)):
                    if orig != mut:
                        allele2.append((orig, pos, mut))

            if allele1 or allele2:
                results[enst] = (allele1, allele2)

        return results

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

a = AnalyzeTrio(("jp-scid7a", "jp-scid7b", "jp-scid7c"),
                {"mother": "jp-scid7a", "father": "jp-scid7b", "child": "jp-scid7c"},
                "varprior.db",
                "hg19.fa",
                "exomes_49.vcf", "20130918_ensGene.tab", "20130918_ensemblToGeneName.tab", "human-protein.aliases.v9.05.txt")

print a.non_synonymous("IL23R", "jp-scid7a")
#print a.mendelian_filter(a._mf_compound_het_denovo, "chr1", 1000000, 2000000, k=2)
#print a.mendelian_filter(a._mf_recessive, "chr1", 1000000, 2000000)
#print a.mendelian_filter(a._mf_denovo_dominant, "chr1", 1000000, 2000000)
