import os
import math
import matplotlib.pyplot as plt
import pyfasta
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import sqlite3
import numpy
import scipy
import scipy.stats
import time
import itertools
import sys
import networkx
from networkx.exception import NetworkXNoPath
import tempfile
import subprocess

class AnalyzeTrio():
  def __init__(self, sample_names, pedigree, varprior_db, network_pickle,
               genome_fasta, vcf_file = None, ensgene_file = None,
               enst_to_gene_name_file = None, string_alias_file = None,
               string_links_file = None):
    self.weights = {"mendelian": 1,
                    "non-synonymous": 1,
                    "network": 1,
                    "global_af": 1,
                    "local_af": 1}
    self.network_score_cutoff = 677
    self.top_genes = ["ADAMTS8", "AIRE", "AK2", "ATM", "BTK", "CD247", "CD3D",
                      "CD3G", "CD40LG", "CD8A", "CD8B", "CHD7", "CIITA",
                      "CORO1A", "CYBB", "DCLRE1C", "DKC1", "DOCK8", "FOXN1",
                      "IKBKG", "IL2RA", "IL2RG", "IL7R", "ITK", "JAK3", "LCK",
                      "LIG4", "NBN", "NHEJ1", "ORAI1", "PNPLA2", "PRKDC",
                      "PTPRC", "RAG1", "RAG2", "RFXANK", "SH2D1A", "STAT5B",
                      "STIM1", "TAP1", "TAP2", "TAPBP", "TBX1", "WAS", "XIAP",
                      "ZAP70", "ZBTB1"]
    self.global_vcf = "AFR.2of4intersection_allele_freq" \
                      ".20100804.genotypes.vcf.gz"

    self.sample_names = sample_names
    self.stripped_samples = [x.replace("-", "") for x in sample_names]

    self.varprior_db = varprior_db
    self.network_pickle = network_pickle
    self.genome_fasta = genome_fasta

    self.vcf_file = vcf_file
    self.ensgene_file = ensgene_file
    self.enst_to_gene_name_file = enst_to_gene_name_file
    self.string_alias_file = string_alias_file
    self.string_links_file = string_links_file

    self.pedigree = pedigree
    self.stripped_pedigree = dict([(x, y.replace("-", "")) for \
                                   x, y in pedigree.items()])

    sys.stderr.write("loading db\n")
    self.conn, self.c, new_db = self.load_db()

    if new_db:
      sys.stderr.write("loading annotation\n")
      self.load_annotation(ensgene_file, enst_to_gene_name_file,
                           string_alias_file)

      sys.stderr.write("loading vcf\n")
      self.load_trio_vcf(vcf_file)

    if os.path.isfile(network_pickle):
      sys.stderr.write("loading network (pickle)\n")
      self.gene_graph = networkx.read_gpickle(network_pickle)
    else:
      sys.stderr.write("loading network\n")
      self.gene_graph = self.load_string_network(string_links_file)

      networkx.write_gpickle(self.gene_graph, network_pickle)

    sys.stderr.write("loading genome\n")
    self.pyf_genome = self.load_genome()

  ##
  ## FILE/DATABASE METHODS
  ##

  def load_db(self):
    if os.path.isfile(self.varprior_db):
      conn = sqlite3.connect(self.varprior_db)
      conn.row_factory = sqlite3.Row
      c = conn.cursor()
    
      try:
        idx_version = c.execute(
          "SELECT value FROM metadata WHERE key = 'version'").fetchone()[0]

        if idx_version != "varprior-1.0":
          raise ValueError(
            "ensGene version (%s) incompatible with this version of Varprior" %
            idx_version)

        record_count = int(c.execute(
          "SELECT value FROM metadata WHERE key = 'records'").fetchone()[0])

        if record_count == "-1":
          raise ValueError("Unfinished/partial database provided")

        records_found = int(c.execute(
          "SELECT COUNT(*) FROM ensGene").fetchone()[0])

        if record_count <> records_found:
          raise ValueError(
            "Corrupt index. Expected %s records, found %s" % (record_count,
                                                              records_found))
      except (sqlite3.OperationalError, sqlite3.DatabaseError), error:
        raise ValueError("Problem with SQLite database: %s" % error)

      new_db = False
    else:
      conn = sqlite3.connect(self.varprior_db)
      conn.row_factory = sqlite3.Row
      c = conn.cursor()

      c.execute("CREATE TABLE metadata (key text, value text)")
      c.execute("CREATE TABLE ensGene (bin int, name text, chrom text, " \
                "strand text, txStart int, txEnd int, cdsStart int, " \
                "cdsEnd int, exonCount int, exonStarts text, exonEnds text, " \
                "score int, name2 text, cdsStartStat text, cdsEndStat text, " \
                "exonFrames text)")
      c.execute("CREATE TABLE ensp_to_enst (ensp text, enst text)")   
      c.execute("CREATE TABLE enst_to_gene_name (enst text, gene_name text)")
      c.execute("CREATE TABLE variants (variant_id integer primary key, chrom text, pos int, %s)" %
                ", ".join(["%s text" % x.replace("-", "") for x in self.sample_names]))
      c.execute("CREATE TABLE gene_tests (gene_name text, net_conn_score float, " \
                "net_conn_perc float, net_rel_score float, net_rel_perc float)")
      c.execute("CREATE TABLE tx_mendel (enst text, model text, varid1 int, " \
                "allele1 text, varid2 int, allele2 text)")
      c.execute("CREATE TABLE variant_tests (variant_id int, allele text, " \
                "local_af float, global_af float)")
      c.execute("CREATE TABLE variant_nonsyn (variant_id int, allele text, " \
                "enst text, mut text)")
 
      c.execute("INSERT INTO metadata VALUES ('version', 'varprior-1.0')")
      c.execute("INSERT INTO metadata VALUES ('records', '-1')")
      c.execute("INSERT INTO metadata VALUES ('variants', '-1')")

      c.execute("CREATE INDEX IF NOT EXISTS ensgene_index ON ensGene(name);")
      c.execute("CREATE INDEX IF NOT EXISTS ensgene_index2 ON ensGene(chrom);")
      c.execute("CREATE INDEX IF NOT EXISTS ensgene_index3 ON ensGene(cdsStart);")
      c.execute("CREATE INDEX IF NOT EXISTS ensgene_index4 ON ensGene(cdsEnd);")
      c.execute("CREATE INDEX IF NOT EXISTS ensgene_index5 ON ensGene(bin);")
      c.execute("CREATE INDEX IF NOT EXISTS chrom_index ON variants(chrom);")
      c.execute("CREATE INDEX IF NOT EXISTS pos_index ON variants(pos);")
      c.execute("CREATE INDEX IF NOT EXISTS ensp_index ON ensp_to_enst(ensp);")
      c.execute("CREATE INDEX IF NOT EXISTS enst_index ON enst_to_gene_name(enst);")
      c.execute("CREATE INDEX IF NOT EXISTS gene_name_index ON enst_to_gene_name(gene_name);")
      c.execute("CREATE INDEX IF NOT EXISTS gene_name_index2 ON gene_tests(gene_name);")
      c.execute("CREATE INDEX IF NOT EXISTS variant_id_index ON variants(variant_id);")
      c.execute("CREATE INDEX IF NOT EXISTS variant_id_index2 ON variant_tests(variant_id);")
      c.execute("CREATE INDEX IF NOT EXISTS variant_id_index3 ON variant_nonsyn(variant_id);")

      new_db = True

    return conn, c, new_db

  def load_genome(self):
    return pyfasta.Fasta(self.genome_fasta,
                         record_class=pyfasta.records.MutNpyFastaRecord)

  def load_annotation(self, ensgene_file, enst_to_gene_name_file,
                      string_alias_file):
    assert self.conn, self.c

    if ensgene_file == None or \
       enst_to_gene_name_file == None or \
       string_alias_file == None:
      raise ValueError("Cannot load annotation without input files")

    # load ensemblToGeneName (translates ENST -> short names, e.g. IL2RG)
    for line in [x.strip().split() for x in open(enst_to_gene_name_file)]:
      self.c.execute(
        "INSERT INTO enst_to_gene_name VALUES ('%s', '%s')" % (line[0], line[1]))

    self.c.execute(
      "INSERT INTO gene_tests (gene_name) SELECT gene_name FROM " \
      "enst_to_gene_name GROUP BY gene_name")
 
    # load string aliases (translates ENSP -> ENST/ENSG)
    for line in [x.strip().split() for x in open(string_alias_file)]:
      if line[2].startswith("ENST"):
        self.c.execute(
          "INSERT INTO ensp_to_enst VALUES ('%s', '%s')" % (line[1], line[2]))

    # load ensgene
    for line_num, line in enumerate([x.strip().split() for x in open(ensgene_file)]):
      if line_num == 0:
        ensgene_keys = line
      else:
        self.c.execute("INSERT INTO ensGene VALUES (%s)" % ", ".join(["'%s'" % x for x in line]))

    self.c.execute(
      "UPDATE metadata SET value = '%s' WHERE key = 'records'" % line_num)
    self.conn.commit()

  def load_trio_vcf(self, vcf_file):
    assert self.conn, self.c

    variant_count = int(self.c.execute(
      "SELECT value FROM metadata WHERE key = 'variants'").fetchone()[0])

    if variant_count > -1:
      if vcf_file != None:
        db_vcf_file = self.c.execute(
          "SELECT value FROM metadata WHERE key = 'vcf_file'").fetchone()[0]

        if vcf_file != db_vcf_file:
          raise ValueError(
            "Database variants do not match specified file (%s != %s)" %
            (db_vcf_file, vcf_file))

        variants_found = int(self.c.execute(
          "SELECT COUNT(*) FROM variants").fetchone()[0])

        if variant_count <> variants_found:
          raise ValueError(
            "Corrupt index. Expected %s records, found %s" % (variant_count,
                                                              variants_found))
    else:
      if vcf_file == None:
        raise ValueError("Database empty and no vcf_file specified")

      records = 0

      for line in open(vcf_file, "r"):
        # parse the header, find columns of interest
        if line.startswith("#"):
          if line.startswith("#CHROM"):
            header = line.strip().split()
    
            try:
              kept_cols = [(x, header.index(x)) for x in self.sample_names]
            except ValueError, error:
              raise ValueError(
                "Specified sample (%s) not in VCF (%s)" % (error.split()[0],
                                                           vcf_file))

          continue

        # parse rest of file
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

        self.c.execute("INSERT INTO variants VALUES (NULL, '%s', '%s', %s)" %
          (chrom, pos,
           ", ".join(["'%s'" % sample_to_snp[x] for x in self.sample_names])))
        records += 1

      self.c.execute("INSERT INTO metadata VALUES ('vcf_file', '%s')" % vcf_file)
      self.c.execute("UPDATE metadata SET value = '%s' WHERE key = 'variants'" % records)
      self.conn.commit()

  @staticmethod
  def region2bin (start, end, join = False):
    """ 
    Converts a region to a list of bins that it may belong to, including largest
    and smallest bins.

    If join is True, will return a comma-delimited list of bins, otherwise will
    return a list.
    """

    bins = [0, 1]

    bins.extend(range(1 + (start >> 26), 2 + ((end - 1) >> 26)))
    bins.extend(range(9 + (start >> 23), 10 + ((end - 1) >> 23)))
    bins.extend(range(73 + (start >> 20), 74 + ((end - 1) >> 20)))
    bins.extend(range(585 + (start >> 17), 586 + ((end - 1) >> 17)))

    if join == True:
      return ", ".join(map(str, set(bins)))
    else:
      return list(set(bins))

  def ensp_to_gene_name(self, ensp):
    assert self.conn, self.c

    gene_names = self.c.execute(
      "SELECT enst_to_gene_name.gene_name FROM enst_to_gene_name, ensp_to_enst "
      "WHERE ensp_to_enst.enst = enst_to_gene_name.enst AND "
      "ensp_to_enst.ensp = '%s'" % ensp)

    return list(set([x["gene_name"] for x in gene_names]))

  def gene_names(self):
    assert self.conn, self.c

    gene_names = self.c.execute(
      "SELECT gene_name FROM enst_to_gene_name GROUP BY gene_name").fetchall()

    for gene in [x["gene_name"] for x in gene_names]:
      yield gene

  def tx_names(self):
    assert self.conn, self.c

    tx_names = self.c.execute(
      "SELECT enst FROM enst_to_gene_name GROUP BY enst").fetchall()

    for enst in ([x["enst"] for x in tx_names]):
      yield enst

  def fetch_tx(self, enst):
    assert self.conn, self.c

    data = self.c.execute("SELECT * FROM ensGene WHERE name = '%s'" % enst).fetchone()

    return data

  def variants(self):
    assert self.conn, self.c

    variants = self.c.execute("SELECT * FROM variants").fetchall()

    for v in variants:
      yield v

  ##
  ## INHERITANCE FILTER METHODS
  ##

  def mendelian_filter(self, m_filter, chrom, start, end, k = 1):
    pos_in_region = self.c.execute(
      "SELECT * FROM variants WHERE chrom = '%s' AND %s <= pos AND pos <= %s" %
      (chrom, start, end)).fetchall()

    hits = []

    for pos in itertools.combinations(pos_in_region, k):
      mother, father, child = {}, {}, {}
      variant_ids = []

      for i in range(k):
        mother[i], father[i], child[i] = \
          [pos[i][self.stripped_pedigree[x]].split("/") for x in \
          ("mother", "father", "child")]

        variant_ids.append(pos["variant_id"])

        # if no variant reported, continue
        if "." in mother[i] or \
           "." in father[i] or \
           "." in child[i]:
          break
      else:
        result = m_filter(mother, father, child):

        if result[0]:
          hits.append(variant_ids, result[1])

    return hits

  @staticmethod
  def _mf_recessive(mother, father, child):
    # homozygous in child, heterozygous in both parents
    if len(set(child[0])) == 1 and \
       len(set(mother[0])) == 2 and \
       len(set(father[0])) == 2:
      return (True, (child[0][0]))

    return (False, None)

  @staticmethod
  def _mf_denovo_dominant(mother, father, child):
    # one or both alleles in child must not be in either parent
    if (child[0][0] not in mother[0] and \
        child[0][0] not in father[0]):
      return (True, (child[0][0]))

    if (child[0][1] not in mother[0] and \
        child[0][1] not in father[0]):
      return (True, child[0][1])

    return (False, None)

  @staticmethod
  def _mf_compound_het_denovo(mother, father, child):
    """
    From http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0070151
        
    1. A variant has to be in a heterozygous state in all affected individuals.
    2. A variant must not occur in a homozygous state in any of the unaffected
       individuals.
    3. A variant that is heterozygous in an affected child must be heterozygous
       in exactly one of the parents.
    3a. The variant must not be heterozygous in both parents.
    3b. The variant must be present in at least one of the parents.
    4. A gene must have two or more heterozygous variants in each of the
       affected individuals.
    5. There must be at least one variant transmitted from the paternal side
       and one transmitted from the maternal side.
    """

    mother1, father1, child1 = mother[0], father[0], child[0]
    mother2, father2, child2 = mother[1], father[1], child[1]

    # 1. child heterozygous at both loci
    if child[0][0] == child[0][1] or \
       child[1][0] == child[1][1]:
      return (False, None)
    
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

      return (True, (child1, child2))

    return (False, None)

  ##
  ## GENOME/CODING POTENTIAL METHODS
  ##

  def non_synonymous_gene(self, gene_name, stripped_sample = False, mut = False):
    # given gene name and set of SNPs, does protein product potentially
    # yield a different protein product from that expected?
    # snp_list = ((pos, mut), (pos, mut), ...)

    # gene_name -> enst
    all_enst = [x[0] for x in self.c.execute(
      "SELECT enst FROM enst_to_gene_name WHERE gene_name = '%s'" % \
      gene_name).fetchall()]

    if len(all_enst) == 0:
      raise ValueError("Gene (%s) not found in database" % gene_name)

    results = {}

    # examine each transcript isoform of this gene
    for enst in all_enst:
      results[enst] = self.non_synonymous_enst(enst, stripped_sample, mut)

  def non_synonymous_tx(self, enst, stripped_sample = False, mut = False):
    # fetch exons
    chrom = self.c.execute(
      "SELECT chrom FROM ensGene WHERE name = '%s'" % enst).fetchone()[0]
    exonStarts, exonEnds = [map(int, x.split(",")[:-1]) for x in self.c.execute(
      "SELECT exonStarts, exonEnds FROM ensGene WHERE name = '%s'" % enst).fetchone()]
    cdsStart, cdsEnd = self.c.execute(
      "SELECT cdsStart, cdsEnd FROM ensGene WHERE name = '%s'" % enst).fetchone()

    # non-coding
    if cdsStart == cdsEnd:
      return ([], [])

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
        break
      # exon contains coding end
      elif start <= cdsEnd < end:
        coding_exons.append((start, cdsEnd))
        break

      # everything else should be in the middle of the coding region
      coding_exons.append((start, end))

    self.pyf_genome[chrom].clearmuts()

    # pull nucleotide sequence of regions
    orig_seq = []

    for start, end in coding_exons:
      orig_seq.append(self.pyf_genome[chrom][start:end])
 
    # make given mutations
    if mut == False and stripped_sample:
      mut_seq1, mut_seq2 = [], []

      for genotype_index, seq_list in ((0, mut_seq1), (1, mut_seq2)):
        pos_in_region = self.c.execute(
          "SELECT * FROM variants WHERE chrom = '%s' AND %s <= pos AND pos <= %s" %
          (chrom, cdsStart, cdsEnd)).fetchall()
    
        for row in pos_in_region:
          seq = row[stripped_sample].split("/")[genotype_index]

          if seq == ".":
            continue

          self.pyf_genome[chrom][row["pos"]] = seq
     
          for start, end in coding_exons:
            seq_list.append(
              self.pyf_genome[chrom].__getitem__(slice(start, end), True))

        # translate original and mutated proteins
        try:
          orig_protein = str(Seq("".join(orig_seq)).translate(table=1))
          mut_protein1 = str(Seq("".join(mut_seq1)).translate(table=1))
          mut_protein2 = str(Seq("".join(mut_seq2)).translate(table=1))
        except TranslationError:
          sys.stderr.write("Transcript %s could not be translated\n" % enst)
          return ([], [])
    
        allele1 = [(orig, pos, mut) for pos, (orig, mut) in \
          enumerate(zip(orig_protein, mut_protein1)) if orig != mut]
        allele2 = [(orig, pos, mut) for pos, (orig, mut) in \
          enumerate(zip(orig_protein, mut_protein2)) if orig != mut]

        return (allele1, allele2)
    elif stripped_sample == False and mut:
      if mut["seq"] != ".":
        self.pyf_genome[chrom][mut["pos"]] = mut["seq"]

      mut_seq = []

      for start, end in coding_exons:
        mut_seq.append(
          self.pyf_genome[chrom].__getitem__(slice(start, end), True))

        try:
          orig_protein = str(Seq("".join(orig_seq)).translate(table=1))
          mut_protein = str(Seq("".join(mut_seq)).translate(table=1))
        except TranslationError:
          sys.stderr.write("Transcript %s could not be translated\n" % enst)

        diffs = [(orig, pos, mut) for pos, (orig, mut) in \
          enumerate(zip(orig_protein, mut_protein)) if orig != mut]

        return diffs
    else:
      raise ValueError("Must specify a mutation or a sample to pull variants from")

  ##
  ## VCF ALLELE FREQUENCY METHODS
  ##

  def tabix_af_in_region(self, vcf_file, chrom, start, end):
    # retrieve allele frequencies of SNPs in given region
    # of tabix-indexed VCF file
    tabix_out = tempfile.NamedTemporaryFile(delete=False)
    
    tabix = subprocess.Popen(["tabix", "-f", "-h", vcf_file,
                              "%s:%s-%s" % (chrom, start, end)],
                             stdout=tabix_out)
    
    while tabix.poll() == None:
      time.sleep(0.01)
    
    af_table = self.af_from_vcf(tabix_out.name)

    os.unlink(tabix_out.name)

    return af_table

  def af_from_vcf(self, vcf_file):
    # allele frequencies from any given vcf
    af_table = dict([(x, {}) for x in self.pyf_genome.keys()])

    vcftools_out = tempfile.NamedTemporaryFile(delete=False)
    
    vcftools = subprocess.Popen(["vcftools", "--vcf", vcf_file,
                                 "--freq", "--out", vcftools_out.name],
                                stdout=open(os.devnull))
    
    while vcftools.poll() == None:
      time.sleep(0.01)

    if os.path.isfile("%s.frq" % vcftools_out.name):
      with open("%s.frq" % vcftools_out.name) as fp:
        fp.next()
        
        for line in fp:
          line_split = line.strip().split()
        
          for allele, freq in [x.split(":") for x in line_split[4:]]:
            chrom = "chr%s" % line_split[0]

            try:
              af_table[chrom][int(line_split[1])][allele] = float(freq)
            except KeyError:
              af_table[chrom][int(line_split[1])] = {}
              af_table[chrom][int(line_split[1])][allele] = float(freq)

      os.unlink("%s.frq" % vcftools_out.name)

    os.unlink("%s.vcfidx" % vcf_file)
    os.unlink("%s.log" % vcftools_out.name)
    os.unlink(vcftools_out.name)

    for chrom in af_table.keys():
      if len(af_table[chrom]) == 0:
        del af_table[chrom]

    return af_table

  ##
  ## NETWORK METHODS
  ##

  def load_string_network(self, string_network_links):
    # NOTE: scores are 1000 - the score so we can use shortest path algorithms
    # NOTE: we don't use a multigraph and instead update to keep edge weights minimum
    assert self.conn, self.c

    # generate STRING graph with gene names
    gene_graph = networkx.Graph(directed=False)
    sys.stderr.write("    nodes\n")
    gene_graph.add_nodes_from(self.gene_names())
        
    links = [x.strip().split() for x in open(string_network_links)]

    for line_num, line in enumerate(links):
      if float(line[2]) < self.network_score_cutoff:
        continue

      gene1 = self.ensp_to_gene_name(line[0][5:])
      gene2 = self.ensp_to_gene_name(line[1][5:])
        
      for node1, node2 in itertools.product(gene1, gene2):
        weight = 1000 - float(line[2])

        try:
          if gene_graph[node1][node2]["weight"] > weight:
            gene_graph[node1][node2]["weight"] = weight
            continue
        except KeyError:
          gene_graph.add_edge(node1, node2, weight=weight)

          if line_num % 100 == 0:
            num_stars = int(math.ceil(line_num / float(len(links)) * 10))
            progress = ["*"] * num_stars + [" "] * (10 - num_stars)

            sys.stderr.write("\r    edges [%s]" % "".join(progress))

    sys.stderr.write("\n")

    return gene_graph

  def mk_subgraph(self):
    # return subgraph induced by top genes and their first-degree neighbors
    neighbors = list(itertools.chain(*[self.gene_graph.neighbors(x) for x in self.top_genes]))
    subgraph = self.gene_graph.subgraph(neighbors + self.top_genes)

    return subgraph

  @staticmethod
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
    neighbors = list(itertools.chain(*[self.gene_graph.neighbors(x) for x in self.top_genes]))
    subgraph = self.gene_graph.subgraph(neighbors + self.top_genes)
    pos = networkx.graphviz_layout(subgraph, prog="neato")
    labels = dict(zip(neighbors + sys.argv[2:], neighbors + self.top_genes))
    
    networkx.draw_networkx_nodes(subgraph, pos, nodelist=neighbors, node_color="#fdb462")
    networkx.draw_networkx_nodes(subgraph, pos, nodelist=self.top_genes, node_color="#80b1d3")
    
    networkx.draw_networkx_edges(subgraph, pos, alpha=0.1)
    
    networkx.draw_networkx_labels(subgraph, pos, labels)
    
    plt.figure(1, figsize=(16, 16))
    plt.axis("off")
    plt.show()

  def score_all_genes(self, graph):
    # connectivity score, defines how central in the network a gene is
    # relative score, defines how close to the nearest puck gene a gene is
    gene_names = self.gene_names()

    for gene in gene_names:
      all_weights = []
      no_paths = []
  
      if gene not in graph:
        connectivity = -1
        relative = -1
      else:
        # shortest path by weight to each of top genes
        for top_gene in self.top_genes:
          try:
            path = networkx.shortest_path(graph, gene, top_gene, weight="weight")
          except (NetworkXNoPath, KeyError):
            no_paths.append(top_gene)
            continue

          weights = []
 
          for node_from, node_to in [path[i:i+2] for i in range(len(path)-1)]:
            # pseudocount of 1 to prevent log(nada) == log(1)
            weights.append(graph[node_from][node_to]["weight"] + 1)

          all_weights.append((top_gene, weights))

        all_weights_only = [x[1] for x in all_weights]

        # connectivity score
        combined_weights = sum(all_weights_only, [])
        conn_raw_score = sum([math.log(x) for x in combined_weights]) 
        connectivity = conn_raw_score * (len(combined_weights) + 10 * len(no_paths)) ** 0.5

        # relative score
        nearest_gene, shortest_path = sorted(all_weights, key=lambda x: len(x[1]))[0]
        rel_raw_score = sum([math.log(x) for x in shortest_path])
        relative = rel_raw_score * len(shortest_path) ** 0.5

      if connectivity != -1:
        print """
gene: %s
conn: %s
    conn_raw: %s
    steps:    %s
    no_path:  %s
rel:  %s
    rel_raw:  %s
    steps:    %s
    topgene:  %s
    path:     %s
""" % (gene, connectivity, conn_raw_score, len(combined_weights), no_paths,
       relative, rel_raw_score, len(shortest_path), nearest_gene, shortest_path)

      self.c.execute(
        "UPDATE gene_tests SET net_conn_score = '%s', net_rel_score = '%s'" \
        " WHERE gene_name = '%s'" % (connectivity, relative, gene))

    self.conn.commit()
 
  def score_gene_percentile(self, gene_name):
    # edge case: one of top genes
    if gene_name in self.top_genes:
      return (1, 1)

    assert self.conn, self.c

    try:
      self.net_conn_hist
    except AttributeError:
      self.net_conn_hist = numpy.array([x[0] for x in \
        self.c.execute("SELECT net_conn_score FROM gene_tests").fetchall() if \
        x[0] != -1])

    try:
      self.net_rel_hist
    except AttributeError:
      self.net_rel_hist = numpy.array([x[0] for x in \
        self.c.execute("SELECT net_rel_score FROM gene_tests").fetchall() if \
        x[0] != -1])

    net_conn_score, net_rel_score = self.c.execute(
      "SELECT net_conn_score, net_rel_score FROM gene_tests WHERE " \
      "gene_name = '%s'" % gene_name).fetchone()

    # edge case: not in network
    if net_conn_score == -1 or \
       net_rel_score == -1:
      return (0, 0)
    else:
      return (1 - scipy.stats.percentileofscore(self.net_conn_hist, net_conn_score) / 100.0,
              1 - scipy.stats.percentileofscore(self.net_rel_hist, net_rel_score) / 100.0)

  ##
  ## STATISTICS METHODS
  ##

  @staticmethod
  def newell_ikeda(k, pois_lambda, T, w):
    # newell-ikeda poisson distributed scan statistic
    return 1 - numpy.exp(-pois_lambda ** k * w ** (k - 1) * T / \
                         scipy.misc.factorial(k - 1, exact=True))

  def mcda_product(self, score_dict):
    # from http://en.wikipedia.org/wiki/Weighted_product_model
    # return weighted product score, combining scores with self.weights
    total_score = 1

    for k, v in score_dict:
        total_score *= float(v) ** self.weights[k]

    return total_score

  ##
  ## RUN THE ANALYSIS
  ##

  def go_gene(self):
    # raw scores for network placement
    self.score_all_genes(self.mk_subgraph())

    # for every gene
    for gene in self.gene_names():
      # network gene placement score
      net_conn_perc, net_rel_perc = self.score_gene_percentile(gene)

      self.c.execute(
        "UPDATE gene_tests SET net_conn_perc = '%s', net_rel_perc = '%s' WHERE " \
        "gene_name = '%s'" % (net_conn_perc, net_rel_perc, gene)) 

      # for every transcript that belongs to this gene
      gene_tx = self.c.execute(
        "SELECT enst FROM enst_to_gene_name WHERE gene_name = '%s'" % gene).fetchall()

      for tx in gene_tx:
        # fetch non-synonymous variants for this tx
        tx_nonsyn = self.c.execute(
          "SELECT * FROM variant_nonsyn WHERE enst = '%s'" % tx).fetchall()

        var_to_mut = {}

        for row in tx_nonsyn:
          var_to_mut[(row["variant_id"], row["allele"])] = row["mut"]

        # for every inheritance model in this transcript
        tx_models = self.c.execute(
          "SELECT * FROM tx_mendel WHERE enst = '%s'" % tx).fetchall()

        for model in tx_models:
          m = {0: {"varid": model["varid1"],
                   "allele": model["allele1"]},
               1: {"varid": model["varid2"],
                   "allele": model["allele2"]}}

          # model -> p-value

          # allele frequency and non-synonymous
          for i in m:
            m[i]["l_af"], m[i]["g_af"] = self.c.execute(
              "SELECT local_af, global_af FROM variant_tests WHERE " \
              "variant_id = '%s' AND allele = '%s'" %
              (m[i]["varid"], m[i]["allele"])).fetchone()

            m[i]["nonsyn"] = var_to_mut[(m[i]["varid"], m[i]["allele"])]

            # af -> p-value
            # non-syn -> p-value

          # combine model params -> pvalue

        # take best-scoring model for this tx

      # take best-scoring tx for this gene

      # combine with gene network placement score -> final gene score

    self.conn.commit()

  def go_tx(self):
    # for every enst
    for enst in self.tx_names():
      print enst

      # fetch tx data
      enst_data = self.fetch_tx(enst)
      chrom, txStart, txEnd = [enst_data[x] for x in ("chrom", "txStart", "txEnd")]

      # mendelian inheritance patterns
      comphet = self.mendelian_filter(self._mf_compound_het_denovo, chrom,
                                      txStart, txEnd, k=2)
      recessive = self.mendelian_filter(self._mf_recessive, chrom,
                                        txStart, txEnd)
      dominant = self.mendelian_filter(self._mf_denovo_dominant, chrom,
                                       txStart, txEnd)

      for test, results in zip(("comphet", "recessive", "dominant"), \
                               (comphet, recessive, dominant)):
        for result in results:
          varid1 = result[0][0]
          allele1 = result[1][0]

          if len(result[0]) == 2:
            varid2 = result[0][1]
            allele2 = result[1][1]
          else:
            varid2 = "NULL"
            allele2 = "NULL"
 
          self.c.execute(
            "INSERT INTO tx_mendel VALUES ('%s', '%s', %s, '%s', %s, '%s')" %
            (enst, test, varid1, allele1, varid2, allele2))

    self.conn.commit()

  def go_variant(self):
    # load local allele frequencies
    local_af = self.af_from_vcf(self.vcf_file)

    # for every variant
    s_time = time.time()
    processed = 0

    for variant in self.variants():
      # what are the non-reference alleles?
      all_seqs = filter(lambda x: x != ".", 
        sum([variant[x].split("/") for x in self.stripped_samples], []))
      nonref_seqs = set(all_seqs).difference(
        self.pyf_genome[variant["chrom"]][variant["pos"]])

      # non-synonymous mutations
      if nonref_seqs:
        # what transcripts overlap this variant
        enst_overlap = [x["name"] for x in self.c.execute(
          "SELECT name FROM ensGene WHERE bin IN (%s) AND chrom = '%s' AND " \
          "cdsStart <= %s AND %s <= cdsEnd" % (self.region2bin(variant["pos"],
                                                               variant["pos"] + 1,
                                                               True),
                                               variant["chrom"],
                                               variant["pos"],
                                               variant["pos"])).fetchall()]

        for enst in enst_overlap:
          for seq in nonref_seqs:
            muts = self.non_synonymous_tx(enst, mut = {"pos": variant["pos"], "seq": seq})

            for mut in muts:
              self.c.execute(
                "INSERT INTO variant_nonsyn (variant_id, allele, enst, mut) " \
                "VALUES (%s, '%s', '%s', '%s')" % (variant["variant_id"],
                                                   seq,
                                                   enst,
                                                   "".join(map(str, mut))))

      # allele freqencies
      global_af = self.tabix_af_in_region(self.global_vcf,
                                          variant["chrom"].replace("chr", ""),
                                          variant["pos"],
                                          variant["pos"] + 1)

      for seq in nonref_seqs:
        # local
        try:
          local_seq_af = local_af[variant["chrom"]][variant["pos"]][seq]
        except KeyError:
          local_seq_af = -1

        # global
        try:
          global_seq_af = global_af[variant["chrom"]][variant["pos"]][seq]
        except KeyError:
          global_seq_af = -1

        if int(self.c.execute(
           "SELECT COUNT(*) FROM variant_tests WHERE variant_id = '%s' AND allele = '%s'" %
           (variant["variant_id"], seq)).fetchone()[0]) == 0:
          self.c.execute(
            "INSERT INTO variant_tests (variant_id, allele, local_af, global_af) " \
            "VALUES (%s, '%s', %s, %s)" %
            (variant["variant_id"], seq, local_seq_af, global_seq_af))
        else:
          self.c.execute(
            "UPDATE variant_tests SET local_af = '%s', global_af = '%s' WHERE " \
            "variant_id = '%s' AND allele = '%s'" % 
            (local_seq_af, global_seq_af, variant["variant_id"], seq))

      processed += 1

      if processed % 100 == 0:
        self.conn.commit()
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    sys.stderr.write("\n")
    self.conn.commit()

a = AnalyzeTrio(("jp-scid7a", "jp-scid7b", "jp-scid7c"),
                {"mother": "jp-scid7a", "father": "jp-scid7b", "child": "jp-scid7c"},
                "varprior.db",
                "varprior.gpickle",
                "data/hg19.fa",
                "data/exomes_49.vcf",
                "data/20130918_ensGene.tab", "data/20130918_ensemblToGeneName.tab",
                "data/human-protein.aliases.v9.05.txt", "data/human-protein.links.v9.05.txt")
