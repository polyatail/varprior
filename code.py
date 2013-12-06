from multiprocessing import Pool
import pysam
from functools import partial
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

NUM_PROCS = 10

##
## GENE NETWORK SCORE FUNCTION -- MUST BE OUTSIDE CLASS FOR MP
##

def score_gene(gene, graph, top_genes): 
  if gene not in graph or \
     graph.degree(gene) == 0:
    harmonic_centrality = -1
    nearest_neighbor = -1
  else:
    dist_to_top_genes = []

    # shortest path by weight to each of top genes
    for top_gene in top_genes:
      if top_gene == gene:
        continue

      try:
        path = networkx.shortest_path(graph, gene, top_gene, weight="weight")
      except (NetworkXNoPath, KeyError):
        continue

      weights = []
 
      for node_from, node_to in [path[i:i+2] for i in range(len(path)-1)]:
        weights.append(graph[node_from][node_to]["weight"])

      dist_to_top_genes.append(float(sum(weights)))

    harmonic_centrality = sum([1 / x for x in dist_to_top_genes])
    nearest_neighbor = sorted(dist_to_top_genes)[0]

    print "gene:  %s\n  c:   %s\n  n:   %s" % (gene, harmonic_centrality, nearest_neighbor)

  return (gene, harmonic_centrality, nearest_neighbor)

def load_trio():
  globals()["a"] = AnalyzeTrio(
    ("jp-scid7a", "jp-scid7b", "jp-scid7c"),
    {"mother": "jp-scid7a", "father": "jp-scid7b", "child": "jp-scid7c"},
    "varprior.db",
    "varprior.gpickle",
    "data/hg19.fa",
    "data/exomes_49.vcf",
    "data/20130918_ensGene.tab", "data/20130918_ensemblToGeneName.tab",
    "data/human-protein.aliases.v9.05.txt", "data/human-protein.links.v9.05.txt",
    "data/evs.txt")

##
## MAIN CLASS
##

class AnalyzeTrio():
  def __init__(self, sample_names, pedigree, varprior_db, network_pickle,
               genome_fasta, vcf_file = None, ensgene_file = None,
               enst_to_gene_name_file = None, string_alias_file = None,
               string_links_file = None, evs_file = None):
    self.weights = {"mendel": 1,
                    "nonsyn": 1,
                    "qv": 1,
                    "global_af": 1,
                    "local_af": 1,
                    "net_cent": 1,
                    "net_nn": 1}
    self.network_score_cutoff = 677
    self.exome_size = 30 * 10 ** 6
    self.top_genes = ["ADAMTS8", "AIRE", "AK2", "ATM", "BTK", "CD247", "CD3D",
                      "CD3G", "CD40LG", "CD8A", "CD8B", "CHD7", "CIITA",
                      "CORO1A", "CYBB", "DCLRE1C", "DKC1", "DOCK8", "FOXN1",
                      "IKBKG", "IL2RA", "IL2RG", "IL7R", "ITK", "JAK3", "LCK",
                      "LIG4", "NBN", "NHEJ1", "ORAI1", "PNPLA2", "PRKDC",
                      "PTPRC", "RAG1", "RAG2", "RFXANK", "SH2D1A", "STAT5B",
                      "STIM1", "TAP1", "TAP2", "TAPBP", "TBX1", "WAS", "XIAP",
                      "ZAP70", "ZBTB1"]
    self.global_vcf = "data/AFR.2of4intersection_allele_freq" \
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

    sys.stderr.write("loading genome\n")
    self.pyf_genome = self.load_genome()

    if new_db:
      sys.stderr.write("loading annotation\n")
      self.load_annotation(ensgene_file, enst_to_gene_name_file,
                           string_alias_file)

      sys.stderr.write("loading vcf\n")
      self.load_trio_vcf(vcf_file)

      sys.stderr.write("loading evs data\n")
      self.load_evs(evs_file)

    if os.path.isfile(network_pickle):
      sys.stderr.write("loading network (pickle)\n")
      self.gene_graph = networkx.read_gpickle(network_pickle)
    else:
      sys.stderr.write("loading network\n")
      self.gene_graph = self.load_string_network(string_links_file)

      networkx.write_gpickle(self.gene_graph, network_pickle)

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
                ", ".join(["%s_1 text, %s_2 text, %s_QV" % ((x.replace("-", ""),) * 3) for x in self.sample_names]))
      c.execute("CREATE TABLE gene_tests (gene_name text, net_cent_score float, " \
                "net_cent_perc float, net_nn_score float, net_nn_perc float)")
      c.execute("CREATE TABLE tx_mendel (enst text, model text, varid1 int, " \
                "allele1 text, varid2 int, allele2 text)")
      c.execute("CREATE TABLE variant_tests (variant_id int, allele text, " \
                "local_af float, global_af float)")
      c.execute("CREATE TABLE variant_nonsyn (variant_id int, allele text, " \
                "enst text, mut text)")
      c.execute("CREATE TABLE evs_pos (evs_pos_id integer primary key, " \
                "chrom text, pos int, phastcons float)")
      c.execute("CREATE TABLE evs_alleles (evs_pos_id int, allele text, af float)")
      c.execute("CREATE TABLE evs_muts (evs_pos_id int, gene_name text, " \
                "tx_name text, status text, mut_aa text, mut_nt text, polyphen text)")
 
      c.execute("INSERT INTO metadata VALUES ('version', 'varprior-1.0')")
      c.execute("INSERT INTO metadata VALUES ('records', '-1')")
      c.execute("INSERT INTO metadata VALUES ('variants', '-1')")

      c.execute("CREATE INDEX IF NOT EXISTS ensp_ensp_to_enst ON ensp_to_enst(ensp);")
      c.execute("CREATE INDEX IF NOT EXISTS enst_enst_to_gene_name ON enst_to_gene_name(enst);")
      c.execute("CREATE INDEX IF NOT EXISTS gene_name_enst_to_gene_name ON enst_to_gene_name(gene_name);")
      c.execute("CREATE INDEX IF NOT EXISTS gene_tests_gene_name ON gene_tests(gene_name);")
      c.execute("CREATE INDEX IF NOT EXISTS variant_tests_variant_id_allele ON variant_tests(variant_id, allele);")
      c.execute("CREATE INDEX IF NOT EXISTS variant_nonsyn_variant_id ON variant_nonsyn(variant_id);")
      c.execute("CREATE INDEX IF NOT EXISTS enst_tx_mendel ON tx_mendel(enst);")

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

    # load string aliases (translates ENSP -> ENST/ENSG)
    for line in [x.strip().split() for x in open(string_alias_file)]:
      if line[2].startswith("ENST"):
        self.c.execute(
          "INSERT INTO ensp_to_enst VALUES ('%s', '%s')" % (line[1], line[2]))

    # load ensgene
    s_time = time.time()
    processed = 0
    batch = []

    for line in open(ensgene_file, "r"):
      line_split = line.strip().split()

      if processed == 0:
        ensgene_keys = line_split
      else:
        batch.append(line_split)

      processed += 1
      if processed % 100 == 0:
        self.c.executemany("INSERT INTO ensGene VALUES (%s)" % (",".join(("?",) * len(ensgene_keys))), batch)
        self.conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.c.executemany("INSERT INTO ensGene VALUES (%s)" % (",".join(("?",) * len(ensgene_keys))), batch)
    sys.stderr.write("\nindexing\n")
    self.c.execute("UPDATE metadata SET value = '%s' WHERE key = 'records'" % (processed - 1))
    self.c.execute("CREATE INDEX IF NOT EXISTS ensgene_name ON ensGene(name);")
    self.c.execute("CREATE INDEX IF NOT EXISTS ensgene_chrom_bin ON ensGene(chrom, bin);")
    self.c.execute("CREATE INDEX IF NOT EXISTS ensgene_cdsStart ON ensGene(cdsStart);")
    self.c.execute("CREATE INDEX IF NOT EXISTS ensgene_cdsEnd ON ensGene(cdsEnd);")
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

      s_time = time.time()
      processed = 0

      batch = []
      cols = ["chrom", "pos"] + sum([["%s_1" % x, "%s_2" % x, "%s_QV" % x] for x in self.stripped_samples], [])

      for line in open(vcf_file, "r"):
        # parse the header, find columns of interest
        if line.startswith("#"):
          if line.startswith("#CHROM"):
            header = line.strip().split()
    
            try:
              kept_cols = [(x.replace("-", ""), header.index(x)) for x in self.sample_names]
            except ValueError, error:
              raise ValueError(
                "Specified sample (%s) not in VCF (%s)" % (error.split()[0],
                                                           vcf_file))

          continue

        # parse rest of file
        line_split = line.strip().split()
        chrom = "chr%s" % line_split[0]

        sample_to_data = {}
    
        for sample, column in kept_cols: 
          pos = int(line_split[1])
          ref = line_split[3]
          alt = [(str(x + 1), y) for x, y in enumerate(line_split[4].split(","))]
        
          code = dict([(str(0), ref)] + alt)

          f_names = line_split[8].split(":")
          f_data = line_split[column].split(":")

          col_data = dict(zip(f_names, f_data))
          col_data["GT"] = reduce(lambda x, y: x.replace(y, code[y]), code, col_data["GT"]).split("/", 1)

          if "." in col_data["GT"]:
            break

          sample_to_data[sample] = col_data
        else:
          data = [chrom, pos] + sum([[sample_to_data[x]["GT"][0],
            sample_to_data[x]["GT"][1], sample_to_data[x]["GQ"] if "GQ" in sample_to_data[x] \
            else "0"] for x in self.stripped_samples], [])

          batch.append(data)
          processed += 1

        if processed % 100 == 0:
          self.c.executemany("INSERT INTO variants (%s) VALUES (%s)" %
            (", ".join(["'%s'" % x for x in cols]),
             ", ".join(["?" for _ in cols])), batch)
          self.conn.commit()
          batch = []
          sys.stderr.write("\rprocessed %s @ %.02f/s" %
            (processed, processed / (time.time() - s_time)))

      self.c.executemany("INSERT INTO variants (%s) VALUES (%s)" %
        (", ".join(["'%s'" % x for x in cols]),
         ", ".join(["?" for _ in cols])), batch)
      self.c.execute("INSERT INTO metadata VALUES ('vcf_file', '%s')" % vcf_file)
      self.c.execute("UPDATE metadata SET value = '%s' WHERE key = 'variants'" % processed)
      sys.stderr.write("\nindexing\n")
      self.c.execute("CREATE INDEX IF NOT EXISTS variants_chrom_pos ON variants(chrom, pos);")
      self.c.execute("CREATE INDEX IF NOT EXISTS variants_variant_id ON variants(variant_id);")
      self.conn.commit()

  def load_evs(self, evs_file):
    assert self.conn, self.c

    evs = dict([(x, {}) for x in self.pyf_genome.keys()])

    s_time = time.time()
    processed = 0

    for line in open(evs_file):
      if line.startswith("#"):
        continue

      line_split = line.strip().split()

      chrom, pos = line_split[0].split(":")
      chrom = "chr%s" % chrom

      try:
        evs[chrom][pos]
      except KeyError:
        evs[chrom][pos] = {}

      try:
        evs[chrom][pos]["phastcons"] = float(line_split[18])
      except ValueError:
        evs[chrom][pos]["phastcons"] = -1

      try:
        evs[chrom][pos]["muts"]
      except KeyError:
        evs[chrom][pos]["muts"] = []

      if line_split[12] != "none":
        evs[chrom][pos]["muts"].append(line_split[12:17] + [line_split[21]])

      try:
        evs[chrom][pos]["alleles"]
      except KeyError:
        evs[chrom][pos]["alleles"] = {}

        af = dict([(x, int(y)) for x, y in [x.split("=") for x in line_split[6].split("/")]])
        total_count = float(sum(af.values()))

        if "R" in af:
          o2n = dict([("A%s" % (x+1), y) for x, y in enumerate([x.split(">")[1] \
            for x in line_split[3].split(";")])])
          o2n["R"] = line_split[3].split(">")[0]

          new_af = dict([(o2n[x], y) for x, y in af.items()])
          af = new_af

        for allele, count in af.items():
          evs[chrom][pos]["alleles"][allele] = count / total_count

        processed += 1
        if processed % 100 == 0:
          sys.stderr.write("\rloaded %s @ %.02f/s" %
            (processed, processed / (time.time() - s_time)))

    sys.stderr.write("\n")
    s_time = time.time()
    processed = 0

    for chrom in sorted(evs.keys()):
      for pos in sorted(evs[chrom].keys()):
        self.c.execute(
          "INSERT INTO evs_pos (chrom, pos, phastcons) VALUES " \
          "('%s', %s, %s)" % (chrom, pos, evs[chrom][pos]["phastcons"]))

        evs_pos_id = self.c.lastrowid

        for allele in evs[chrom][pos]["alleles"]:
          self.c.execute(
            "INSERT INTO evs_alleles (evs_pos_id, allele, af) VALUES " \
            "(%s, '%s', %s)" % 
            (evs_pos_id, allele, evs[chrom][pos]["alleles"][allele]))

        for mut in evs[chrom][pos]["muts"]:
          self.c.execute(
            "INSERT INTO evs_muts (evs_pos_id, gene_name, tx_name, status, " \
            "mut_aa, mut_nt, polyphen) VALUES (%s)" %
            (", ".join(["'%s'" % x for x in [evs_pos_id] + mut])))

        processed += 1
        if processed % 100 == 0:
          self.conn.commit()
          sys.stderr.write("\rinserted %s @ %.02f/s" %
            (processed, processed / (time.time() - s_time)))

    sys.stderr.write("\nindexing\n")
    self.c.execute("CREATE INDEX IF NOT EXISTS evs_pos_chrom_pos ON evs_pos(chrom, pos);")
    self.c.execute("CREATE INDEX IF NOT EXISTS evs_pos_evs_pos_id ON evs_pos(evs_pos_id);")
    self.c.execute("CREATE INDEX IF NOT EXISTS evs_alleles_evs_pos_id ON evs_alleles(evs_pos_id);")
    self.c.execute("CREATE INDEX IF NOT EXISTS evs_muts_evs_pos_id ON evs_muts(evs_pos_id);")
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
        for var, name in zip((mother, father, child), ("mother", "father", "child")):
          var[i] = [pos[i]["%s%s" % (self.stripped_pedigree[name], x)] for x in ("_1", "_2")]

        variant_ids.append(pos[i]["variant_id"])
      else:
        result = m_filter(mother, father, child)

        if result[0]:
          hits.append((variant_ids, result[1]))

    return hits

  @staticmethod
  def _mf_dummy(mother, father, child):
    return (False, None)

  @staticmethod
  def _mf_recessive(mother, father, child):
    # homozygous in child, heterozygous in both parents
    if len(set(child[0])) == 1 and \
       len(set(mother[0])) == 2 and \
       len(set(father[0])) == 2:
      return (True, (child[0][0],))

    return (False, None)

  @staticmethod
  def _mf_denovo_dominant(mother, father, child):
    # one or both alleles in child must not be in either parent
    if (child[0][0] not in mother[0] and \
        child[0][0] not in father[0]):
      return (True, (child[0][0],))

    if (child[0][1] not in mother[0] and \
        child[0][1] not in father[0]):
      return (True, (child[0][1],))

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
    if child1[0] == child1[1] or \
       child2[0] == child2[1]:
      return (False, None)
    
    # there are four hypotheses. each position has two alleles,
    # and any combination of the two alleles at each position
    # could be causative
    for c1, c2 in itertools.product(child1, child2):
      # 2. parents can't be homozygous for either of these variants
      if (mother1[0] == c1 and \
          mother1[1] == c1) or \
         (father1[0] == c1 and \
          father1[1] == c1):
        continue
    
      if (mother2[0] == c2 and \
          mother2[1] == c2) or \
         (father2[0] == c2 and \
          father2[1] == c2):
        continue
    
      # 3. causative variant must be heterozygous in one or neither parent
      if (c1 in mother1 and \
          c1 in father1) or \
         (c2 in mother2 and \
          c2 in father2):
        continue
   
      # 5. both variants can't come from the same parent
      if (c1 in mother1 and \
          c2 in mother2) or \
         (c1 in father1 and \
          c2 in father2):
        continue

      return (True, (c1, c2))

    return (False, None)

  ##
  ## GENOME/CODING POTENTIAL METHODS
  ##

  def non_synonymous_gene(self, gene_name, mut):
    # gene_name -> enst
    all_enst = [x[0] for x in self.c.execute(
      "SELECT enst FROM enst_to_gene_name WHERE gene_name = '%s'" % \
      gene_name).fetchall()]

    if len(all_enst) == 0:
      raise ValueError("Gene (%s) not found in database" % gene_name)

    results = {}

    # examine each transcript isoform of this gene
    for enst in all_enst:
      results[enst] = self.non_synonymous_tx(enst, mut)

  def non_synonymous_tx(self, enst, mut):
    # fetch exons
    chrom, strand, exonStarts, exonEnds, cdsStart, cdsEnd = self.c.execute(
      "SELECT chrom, strand, exonStarts, exonEnds, cdsStart, cdsEnd FROM " \
      "ensGene WHERE name = '%s'" % enst).fetchone()

    # non-coding
    if cdsStart == cdsEnd:
      return []

    exonStarts, exonEnds = [map(int, x.split(",")[:-1]) for x in (exonStarts, exonEnds)]

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
    self.pyf_genome[chrom][mut["pos"]] = mut["seq"]

    mut_seq = []

    for start, end in coding_exons:
      mut_seq.append(self.pyf_genome[chrom].__getitem__(slice(start, end), True))

    try:
      if strand == "-":
        orig_protein = str(Seq("".join(orig_seq)).reverse_complement().translate(table=1))
        mut_protein = str(Seq("".join(mut_seq)).reverse_complement().translate(table=1))
      else:
        orig_protein = str(Seq("".join(orig_seq)).translate(table=1))
        mut_protein = str(Seq("".join(mut_seq)).translate(table=1))

      diffs = [(o, p, m) for p, (o, m) in \
        enumerate(zip(orig_protein, mut_protein)) if o != m]
    except TranslationError:
      sys.stderr.write("Transcript %s could not be translated\n" % enst)

      diffs = []

    return diffs

  ##
  ## VCF ALLELE FREQUENCY METHODS
  ##

  def tabix_af_in_region_native(self, vcf_file, chrom, start, end):
    tbx = pysam.Tabixfile(vcf_file)
    lines = tbx.fetch(int(chrom), int(start), int(end))

    all_af = {}

    for l in lines:
      l = l.strip().split("\t")

      alleles = [l[3]] + l[4].split(",")
      af = dict([(x, 0) for x in alleles if x != "."])

      for indiv in l[9:]:
        try:
          af[alleles[int(indiv[0])]] += 1
        except ValueError:
          pass

        try:
          af[alleles[int(indiv[2])]] += 1
        except ValueError:
          pass

      total = float(sum(af.values()))
      af = dict([(x, y / total) for x, y in af.items()])
      all_af[int(l[1])] = af

    return all_af

  def tabix_af_in_region(self, vcf_file, chrom, start, end):
    # retrieve allele frequencies of SNPs in given region
    # of tabix-indexed VCF file
    tabix_out = tempfile.NamedTemporaryFile(delete=False)
    
    tabix = subprocess.Popen(["tabix", "-f", "-h", vcf_file,
                              "%s:%s-%s" % (chrom, start, end)],
                             stdout=tabix_out)
    
    while tabix.poll() == None:
      time.sleep(0.01)

    if os.stat(tabix_out.name).st_size == 0:
      af_table = {}
    else:
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
    # two scores:
    # harmonic centrality defines how well connected a gene is all puck genes
    # nearest-neighobr defines how close a gene is to its nearest puck gene

    p = Pool(NUM_PROCS)
    partial_score_gene = partial(score_gene, graph=graph, top_genes=self.top_genes)
    result = p.map(partial_score_gene, self.gene_names())
    p.close()

    # and now go through and convert them all to percentiles
    cent_hist = numpy.array([x[1] for x in result if x[1] != -1])
    nn_hist = numpy.array([x[2] for x in result if x[2] != -1])

    batch = []

    for gene, cent_score, nn_score in result:
      # edge case: gene is a top gene
      if gene in self.top_genes:
        cent_perc = 1
        nn_perc = 1
      # edge case: gene isn't in network
      elif cent_score == -1 or \
           nn_score == -1:
        cent_perc = 0
        nn_perc = 0
      else:
        cent_perc = scipy.stats.percentileofscore(cent_hist, cent_score) / 100.0
        nn_perc = 1 - scipy.stats.percentileofscore(nn_hist, nn_score) / 100.0

        print """
gene:  %s
  c:   %s
  c_p: %s
  n:   %s
  n_p: %s
""" % (gene, cent_score, cent_perc, nn_score, nn_perc)

      batch.append((gene, cent_score, cent_perc, nn_score, nn_perc))

    self.c.executemany("INSERT INTO gene_tests (gene_name, net_cent_score, " \
      "net_cent_perc, net_nn_score, net_nn_perc) VALUES (?,?,?,?,?)", batch)
    self.conn.commit()

  ##
  ## STATISTICS METHODS
  ##

  @staticmethod
  def newell_ikeda(k, pois_lambda, T, w):
    # newell-ikeda poisson distributed scan statistic
    return 1 - numpy.exp(-pois_lambda ** k * w ** (k - 1) * T / \
                         scipy.misc.factorial(k - 1, exact=True))

  def wpm(self, score_dict):
    # from http://en.wikipedia.org/wiki/Weighted_product_model
    # return weighted product score, combining scores with self.weights
    total_score = 1

    for k, v in score_dict.items():
        total_score *= float(v) ** self.weights[k]

    return total_score

  def wsm(self, score_dict):
    # from http://en.wikipedia.org/wiki/Weighted_sum_model
    # return weighted sum score, combining scores with self.weights
    total_score = []

    for k, v in score_dict.items():
        total_score.append(float(v) * self.weights[k])

    return sum(total_score)

  ##
  ## NATIVE AND EXOME VARIANT SERVER OR OTHER TESTS
  ##

  def nonsyn_native(self, variant, nonref_seqs):
    # what transcripts overlap this variant
    enst_overlap = [x["name"] for x in self.c.execute(
      "SELECT name FROM ensGene WHERE chrom = '%s' AND bin IN (%s) AND " \
      "cdsStart <= %s AND %s <= cdsEnd" % (variant["chrom"],
                                           self.region2bin(variant["pos"],
                                                           variant["pos"] + 1,
                                                           True),
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

  def nonsyn_evs(self, variant, nonref_seqs):
    pass

  def af_native(self, variant, all_seqs):
    try:
      self.local_af
    except AttributeError:
      self.local_af = self.af_from_vcf(self.vcf_file)

    global_af = self.tabix_af_in_region_native(
      self.global_vcf, variant["chrom"].replace("chr", ""),
      variant["pos"], variant["pos"] + 1)

    for seq in all_seqs:
      # local
      try:
        local_seq_af = self.local_af[variant["chrom"]][variant["pos"]][seq]
      except KeyError:
        local_seq_af = -1

      # global
      try:
        global_seq_af = global_af[variant["chrom"]][variant["pos"]][seq]
      except KeyError:
        global_seq_af = -1

      self.c.execute(
        "UPDATE variant_tests SET local_af = '%s', global_af = '%s' WHERE " \
        "variant_id = '%s' AND allele = '%s'" % 
        (local_seq_af, global_seq_af, variant["variant_id"], seq))

  def af_evs(self, variant, all_seqs):
    try:
      self.local_af
    except AttributeError:
      self.local_af = self.af_from_vcf(self.vcf_file)

    global_af_db = self.c.execute(
      "SELECT allele, af FROM evs_pos, evs_alleles WHERE evs_pos.chrom = '%s' AND " \
      "evs_pos.pos = '%s' AND evs_pos.evs_pos_id = evs_alleles.evs_pos_id" %
      (variant["chrom"], variant["pos"])).fetchall()
    global_af = {}

    for row in global_af_db:
      global_af[row["allele"]] = row["af"]

    for seq in all_seqs:
      # local
      try:
        local_seq_af = self.local_af[variant["chrom"]][variant["pos"]][seq]
      except KeyError:
        local_seq_af = -1

      # global
      try:
        global_seq_af = global_af[seq]
      except KeyError:
        global_seq_af = -1

      self.c.execute(
        "UPDATE variant_tests SET local_af = '%s', global_af = '%s' WHERE " \
        "variant_id = '%s' AND allele = '%s'" % 
        (local_seq_af, global_seq_af, variant["variant_id"], seq))

  ##
  ## RUN THE ANALYSIS
  ##

  def go_gene(self):
    # how many recessive, dominant, and comphet models are there
    mendel_counts = {}

    mendel_counts["recessive"] = self.c.execute(
      "SELECT COUNT(*) FROM tx_mendel WHERE model = 'recessive'").fetchone()[0]
    mendel_counts["dominant"] = self.c.execute(
      "SELECT COUNT(*) FROM tx_mendel WHERE model = 'dominant'").fetchone()[0]
    mendel_counts["comphet"] = self.c.execute(
      "SELECT COUNT(*) FROM tx_mendel WHERE model = 'comphet'").fetchone()[0]

    # for every gene
    for gene in self.gene_names():
      # gene network placement
      cent_perc, nn_perc = self.c.execute("SELECT net_cent_perc, net_nn_perc " \
        "FROM gene_tests WHERE gene_name = '%s'" % gene).fetchone()

      # every tx that belongs to this gene
      gene_tx = [x["enst"] for x in self.c.execute(
        "SELECT enst FROM enst_to_gene_name WHERE gene_name = '%s'" % gene).fetchall()]

      print "gene: %s\n  cent: %s\n  nn:   %s" % (gene, cent_perc, nn_perc)

      all_tx_scores = []

      for tx in gene_tx:
        print "    tx: %s" % tx

        # fetch every inheritance model in this transcript
        tx_models = self.c.execute(
          "SELECT * FROM tx_mendel WHERE enst = '%s'" % tx).fetchall()

        if len(tx_models) == 0:
          continue

        # fetch exon size for this transcript
        exonStarts, exonEnds = [map(int, x.split(",")[:-1]) for x in self.c.execute(
          "SELECT exonStarts, exonEnds FROM ensGene WHERE name = '%s'" % tx).fetchone()]

        exon_size = sum([(x - y) for x, y in zip(exonEnds, exonStarts)])

        # fetch non-synonymous variants for this tx
        tx_nonsyn = self.c.execute(
          "SELECT * FROM variant_nonsyn WHERE enst = '%s'" % tx).fetchall()

        var_to_mut = {}

        for row in tx_nonsyn:
          var_to_mut[(row["variant_id"], row["allele"])] = row["mut"]

        # iterate through all the models
        all_model_scores = []

        for model in tx_models:
          print "      model: %s,%s=%s,%s=%s" % (model["model"], model["varid1"],
            model["allele1"], model["varid2"], model["allele2"])

          model_score = {}

          m = {0: {"varid": model["varid1"],
                   "allele": model["allele1"]}}

          if model["varid2"] != "NULL":
            m[1] = {"varid": model["varid2"],
                    "allele": model["allele2"]}

          ni_T = self.exome_size / exon_size
          ni_lambda = mendel_counts[model["model"]] / ni_T

          model_score["mendel"] = self.newell_ikeda(1, ni_lambda, ni_T, 1)

          # allele frequency, QVs, phastcons, and non-synonymous
          for i in m:
            m[i]["qv"] = [float(x) for x in self.c.execute(
              "SELECT %s FROM variants WHERE variant_id = '%s'" %
              (", ".join(["%s_QV" % x for x in self.stripped_samples]),
               m[i]["varid"])).fetchone()]

            #TODO: pull phastcons for this variant from evs_pos table

            m[i]["l_af"], m[i]["g_af"] = self.c.execute(
              "SELECT local_af, global_af FROM variant_tests WHERE " \
              "variant_id = '%s' AND allele = '%s'" %
              (m[i]["varid"], m[i]["allele"])).fetchone()

            try:
              m[i]["nonsyn"] = var_to_mut[(m[i]["varid"], m[i]["allele"])]
            except KeyError:
              m[i]["nonsyn"] = False

          print "      m: %s" % m

          model_score["nonsyn"] = sum([1 for x in m.values() if x["nonsyn"]]) / float(len(m))
          model_score["local_af"] = reduce(lambda x, y: x*y, [x["l_af"] for x in m.values()])
          model_score["global_af"] = reduce(lambda x, y: x*y, [x["g_af"] for x in m.values()])
          model_score["qv"] = reduce(lambda x, y: x*y, sum([x["qv"] for x in m.values()], []))

          all_model_scores.append((self.wsm(model_score), model_score))

        best_model_score = sorted(all_model_scores, key=lambda x: x[0])[-1]

        all_tx_scores.append(best_model_score)

      if len(all_tx_scores) == 0:
        best_tx_score = {}
      else:
        best_tx_score = sorted(all_tx_scores, key=lambda x: x[0])[-1][1]

      best_tx_score["net_cent"] = cent_perc
      best_tx_score["net_nn"] = nn_perc

      gene_score = self.wsm(best_tx_score)

      print gene, gene_score, best_tx_score

      # UPDATE gene_tests SET final_score = '%s' WHERE gene_name = '%s'

    self.conn.commit()

  def go_tx(self):
    s_time = time.time()
    processed = 0

    batch = []

    for enst in self.tx_names():
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
 
          batch.append((enst, test, varid1, allele1, varid2, allele2))

      processed += 1

      if processed % 100 == 0:
        self.c.executemany("INSERT INTO tx_mendel VALUES (?,?,?,?,?,?)", batch)
        self.conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.c.executemany("INSERT INTO tx_mendel VALUES (?,?,?,?,?,?)", batch)
    self.conn.commit()
    sys.stderr.write("\n")

  def go_variant(self):
    s_time = time.time()
    processed = 0

    for variant in self.variants():
      # what are the non-reference alleles?
      all_seqs = set([variant["%s%s" % (x, y)] for x in self.stripped_samples for y in ("_1", "_2")])
      nonref_seqs = all_seqs.difference(
        self.pyf_genome[variant["chrom"]][variant["pos"]])

      # first off, make an entry for each allele
      for seq in all_seqs:
        self.c.execute(
          "INSERT INTO variant_tests (variant_id, allele) VALUES (%s, '%s')" %
          (variant["variant_id"], seq))

      # non-synonymous mutations
      if nonref_seqs:
        self.nonsyn_native(variant, nonref_seqs)
        #self.nonsyn_evs(variant, nonref_seqs)

      # allele freqencies
      #self.af_native(variant, all_seqs)
      self.af_evs(variant, all_seqs)

      processed += 1

      if processed % 100 == 0:
        self.conn.commit()
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.conn.commit()
    sys.stderr.write("\n")
