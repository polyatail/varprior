import pysam
import os
import pyfasta
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import sqlite3
import numpy
import scipy
import time
import itertools
import sys

def load_trio():
  globals()["a"] = AnalyzeTrio(
    {"mother": "jp-scid9b", "father": "jp-scid9c", "daughter": "jp-scid9a"},
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
  def __init__(self, db_file, vcf_file, genome_fasta, pedigree):
    self.version = "analyzetrio-1.0"
    self.weights = {"mendel": 1,
                    "nonsyn": 3,
                    "qv": 2,
                    "global_af": -6,
                    "local_af": -6,
                    "net_cent": 1,
                    "net_nn": 1}
    self.exome_size = 30 * 10 ** 6

    self.db_file = db_file
    self.vcf_file = vcf_file
    self.genome_fasta = genome_fasta

    self.pedigree = pedigree
    self.stripped_pedigree = dict([(x, y.replace("-", "")) for \
                                   x, y in pedigree.items()])

    sys.stderr.write("loading genome\n")
    self.load_genome()

    sys.stderr.write("loading db\n")

    if not os.path.isfile(self.db_file):
      self.db_connect()
      self.db_init()

      sys.stderr.write("loading vcf\n")
      self.load_vcf()
      self.db_metadata("update")
    else:
      self.db_connect()
      self.db_metadata("check")

  ##
  ## DATABASE METHODS
  ##

  def db_connect(self):
    self._conn = sqlite3.connect(self.db_file)
    self._conn.row_factory = sqlite3.Row
    self._c = self._conn.cursor()

  def db_init(self):
    assert self._conn, self._c

    tables = {}

    tables["metadata"] = {"f": (("k", "text"),
                                ("v", "text"))}

    tables["vcf_variants"] = {"f": (("vcf_variant_id", "integer primary key"),
                                    ("chrom", "text"),
                                    ("pos", "int"),
                              "i": ("vcf_variant_id", ("chrom", "pos"))}

    tables["vcf_var_to_var"] = {"f": (("vcf_variant_id", "int"),
                                      ("variant_id", "int")),
                                "i": ("vcf_variant_id", "variant_id")}

    tables["vcf_alleles"] = {"f": (("vcf_allele_id", "integer primary key"),
                                   ("variant_id", "int"),
                                   ("sequence", "text"),
                                   ("local_af", "float"),
                             "i": ("vcf_allele_id", ("variant_id", "sequence"))}

    tables["vcf_all_to_all"] = {"f": (("vcf_allele_id", "int"),
                                      ("allele_id", "int")),
                                "i": ("vcf_allele_id", "allele_id")}

    tables["vcf_muts"] = {"f": (("vcf_mut_id", "integer primary key"),
                                ("vcf_allele_id", "int"),
                                ("tx_id", "int"),
                                ("mut", "text")),
                          "i": ("vcf_mut_id", "vcf_allele_id", "tx_id")}

    inserts = {}

    inserts["metadata"] = {"f": ("k", "v"),
                           "r": (("version", self.version),
                                 ("vcf_file", "-1"),
                                 ("vcf_variants", "-1"),
                                 ("vcf_alleles", "-1"),
                                 ("vcf_muts", "-1"))}

    for t in tables:
      f = ["%s %s" % (x, y) for x, y in tables[t]["f"]]

      self._c.execute("CREATE TABLE %s (%s)" % (t, ", ".join(f)))

      if "i" in tables[t]:
        for i in tables[t]["i"]:
          n = "_".join(i) if isinstance(i, tuple) else i
          f = ", ".join(i) if isinstance(i, tuple) else i

          self._c.execute(
            "CREATE INDEX IF NOT EXISTS %(t)s_%(n)s ON %(t)s(%(f)s)" %
            {"t": t, "n": n, "f": f})

    for t in inserts:
      f = ", ".join(inserts[t]["f"])
      b = ", ".join(["?" for _ in inserts[t]["f"]])

      self._c.executemany("INSERT INTO %(t)s (%(f)s) VALUES (%(b)s)" %
        {"t": t, "f": f, "b": b}, inserts[t]["r"])

    self._conn.commit()

  def db_metadata(self, action):
    assert self._conn, self._c

    md = {"version": ("l", self.version),
          "vcf_file": ("l", self.vcf_file),
          "vcf_variants": ("q", "SELECT COUNT(*) FROM vcf_variants"),
          "vcf_alleles": ("q", "SELECT COUNT(*) FROM vcf_alleles"),
          "vcf_muts": ("q", "SELECT COUNT(*) FROM vcf_muts")}

    if action == "check":
      for k, q in md.items():
        md_v = self._c.execute("SELECT v FROM metadata WHERE k = '%s'" % k).fetchone()[0]

        if q[0] == "q":
          ac_v = self._c.execute(q).fetchone()[0]
        elif q[0] == "l":
          ac_v = q[1]

        if md_v != ac_v:
          raise ValueError("Invalid DB metadata (%s, %s != %s)" % (k, md_v, ac_v))
    elif action == "update":
      for k, q in md.items():
        if q[0] == "q":
          ac_v = self._c.execute(q).fetchone()[0]

          self._c.execute("UPDATE metadata SET v = '%s' WHERE k = '%s'" % (k, ac_v))

  def load_genome(self):
    self.pyf_genome = pyfasta.Fasta(self.genome_fasta,
                         record_class=pyfasta.records.MutNpyFastaRecord)

  ##
  ## FILE LOADING METHODS
  ##

  def load_trio_vcf(self, vcf_file):
    assert self.conn, self.c

    # has the table already been loaded?
    md_v = self._c.execute("SELECT v FROM metadata WHERE k = 'vcf_variants'").fetchone()[0]

    if int(md_v) != -1:
      raise ValueError("Unexpected DB metadata (vcf_variants, %s != -1)" % md_v)

    s_time = time.time()
    processed = 0
    batch = []

    for l in open(vcf_file):
      # parse the header, find columns of interest
      if l.startswith("#"):
        if l.startswith("#CHROM"):
          header = l.strip().split()

          kept_cols = []
          cols = ["variant_id", "chrom", "pos"]

          for i, n in enumerate(header):
            if n.startswith("jp"):
              cols.extend(["%s_1" % n, "%s_2" % n, "%s_QV" % n])
              kept_cols.append((i, n))

              self._c.execute("ALTER TABLE vcf_variants ADD COLUMN %s_1 text" % n)
              self._c.execute("ALTER TABLE vcf_variants ADD COLUMN %s_2 text" % n)
              self._c.execute("ALTER TABLE vcf_variants ADD COLUMN %s_QV int" % n)

        continue

      # parse rest of file
      l = l.strip().split()

      chrom = "chr%s" % l[0]
      pos = int(l[1])
      ref = l[3]
      alt = [(str(x + 1), y) for x, y in enumerate(l[4].split(","))]
      code = dict([(str(0), ref)] + alt)

      f_names = l[8].split(":")

      samples = {}

      for i, n in kept_cols:
        f_data = l[i].split(":")

        col_data = dict(zip(f_names, f_data))
        col_data["GT"] = [code[x] for x in col_data["GT"].split("/", 1)]

        if "." in col_data["GT"]:
          break

        samples[n] = col_data
      else:
        # match variant to static db
        v = self.vd.fetch_variant(chrom, pos)

        row = [-1 if v == None else v.variant_id, chrom, pos]

        for _, n in kept_cols:
          row.append(samples[n]["GT"][0],
                     samples[n]["GT"][1],
                     samples[n]["GQ"] if "GQ" in samples[n] else "0")

        batch.append(row)
        processed += 1

      if processed % 100 == 0:
        self.c.executemany("INSERT INTO vcf_variants (%s) VALUES (%s)" %
          (", ".join(["'%s'" % x for x in cols]),
           ", ".join(["?" for _ in cols])), batch)
        self.conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.c.executemany("INSERT INTO vcf_variants (%s) VALUES (%s)" %
      (", ".join(["'%s'" % x for x in cols]),
       ", ".join(["?" for _ in cols])), batch)
    self.conn.commit()
    sys.stderr.write("\rprocessed %s @ %.02f/s" %
      (processed, processed / (time.time() - s_time)))

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
  def _mf_xlinked(mother, father, child):
    # at least one reported allele in son must not be in father and be present
    # no more than once in the mother
    for c in child[0]:
      if c not in father[0] and \
        (mother[0][0] != c or \
         mother[0][1] != c):
        return (True, (c,))

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
    middle = False

    # reduce coordinates to those that are coding
    for start, end in zip(exonStarts, exonEnds):
      if middle:
        # exon contains coding end
        if start <= cdsEnd < end:
          coding_exons.append((start, cdsEnd))
          break

        # everything else should be in the middle of the coding region
        coding_exons.append((start, end))
      else:
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
          middle = True
          continue

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
      if not k.startswith("_"):
        total_score *= float(v) ** self.weights[k]

    return total_score

  def wsm(self, score_dict):
    # from http://en.wikipedia.org/wiki/Weighted_sum_model
    # return weighted sum score, combining scores with self.weights
    total_score = []

    for k, v in score_dict.items():
      if not k.startswith("_"):
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
        muts = self.non_synonymous_tx(enst, mut = {"pos": variant["pos"] - 1, "seq": seq})

        for mut in muts:
          self.c.execute(
            "INSERT INTO variant_nonsyn (variant_id, allele, enst, mut) " \
            "VALUES (%s, '%s', '%s', '%s')" % (variant["variant_id"],
                                               seq,
                                               enst,
                                               "".join(map(str, mut))))

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

    mendel_counts["xlinked"] = self.c.execute(
      "SELECT COUNT(*) FROM tx_mendel WHERE model = 'xlinked'").fetchone()[0]
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

      #print "gene: %s\n  cent: %s\n  nn:   %s" % (gene, cent_perc, nn_perc)

      all_tx_scores = []

      for tx in gene_tx:
        #print "    tx: %s" % tx

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
          #print "      model: %s,%s=%s,%s=%s" % (model["model"], model["varid1"],
          #  model["allele1"], model["varid2"], model["allele2"])

          model_score = {}

          m = {0: {"varid": model["varid1"],
                   "allele": model["allele1"]}}

          if model["varid2"]:
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

          #print "      m: %s" % m

          model_score["nonsyn"] = sum([1 for x in m.values() if x["nonsyn"]]) / float(len(m))
          model_score["local_af"] = reduce(lambda x, y: x*y, [x["l_af"] for x in m.values()])
          model_score["global_af"] = reduce(lambda x, y: x*y, [x["g_af"] for x in m.values()])
          model_score["global_af"] = 0 if model_score["global_af"] < 0 else model_score["global_af"]
          model_score["qv"] = sum(sum([x["qv"] for x in m.values()], [])) / \
            (594.0 if model["model"] == "comphet" else 297)
          model_score["_model"] = model

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
      if chrom == "chrX" and \
         self.child_gender == "male":
        comphet, recessive, dominant = [], [], []
        xlinked = self.mendelian_filter(self._mf_xlinked, chrom,
                                        txStart, txEnd)
      else:
        xlinked = []
        comphet = self.mendelian_filter(self._mf_compound_het_denovo, chrom,
                                        txStart, txEnd, k=2)
        recessive = self.mendelian_filter(self._mf_recessive, chrom,
                                          txStart, txEnd)
        dominant = self.mendelian_filter(self._mf_denovo_dominant, chrom,
                                         txStart, txEnd)

      for test, results in zip(("xlinked", "comphet", "recessive", "dominant"), \
                               (xlinked, comphet, recessive, dominant)):
        for result in results:
          varid1 = result[0][0]
          allele1 = result[1][0]

          if len(result[0]) == 2:
            varid2 = result[0][1]
            allele2 = result[1][1]
          else:
            varid2 = None
            allele2 = None
 
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
