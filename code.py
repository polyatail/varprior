import bz2file
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

##
## AUXILIARY CLASSES
##

class VCFEntry():
  def __init__(self, vcf_id):
    self.vcf_id = vcf_id

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

    tables["vcf"] = {"f": (("vcf_id", "integer primary key"),
                           ("variant_id", "int")),
                     "i": ("vcf_id", "variant_id")}

    tables["mendel"] = {"f": (("mendel_id", "integer primary key"),
                              ("tx_id", "int"),
                              ("model", "text"),
                              ("allele_id1", "int"),
                              ("allele_id2", "int")),
                        "i": ("mendel_id", "tx_id")}
    inserts = {}

    inserts["metadata"] = {"f": ("k", "v"),
                           "r": (("version", self.version),
                                 ("vcf_file", "-1"),
                                 ("vcf_count", "-1"),
                                 ("mendel_count", "-1"),

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
          "vcf_count": ("q", "SELECT COUNT(*) FROM vcf"),
          "mendel_count": ("q", "SELECT COUNT(*) FROM mendel")}

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

  def fetch_vcf(self, vcf_id = None, variant_id = None):
    assert self._conn, self._c

    if not (vcf_id ^ variant_id):
      raise ValueError("Must specify either vcf_id or variant_id")

    if variant_id:
      d = self._c.execute("SELECT * FROM vcf WHERE variant_id = %s" % variant_id).fetchone()
    elif vcf_id:
      d = self._c.execute("SELECT * FROM vcf WHERE vcf_id = %s" % vcf_id).fetchone()

    if d == None:
      return None

    v = VCFEntry(d["vcf_id"])

    v.variant_id = d["variant_id"]
    v.samples = {}

    for k, v in d.items():
      if k not in ("vcf_id", "variant_id"):
        v.samples[k] = v

    return v

  def fetch_vcf_overlapping(self, chrom, start, end):
    assert self._conn, self._c

    vd_d = self.vd.fetch_variant_overlapping(chrom, start, end)
    vcf = [self.fetch_vcf(variant_id=x) for x in vd_d]

    return vcf

  ##
  ## FILE LOADING METHODS
  ##

  def load_genome(self):
    self.pyf_genome = pyfasta.Fasta(self.genome_fasta,
                         record_class=pyfasta.records.MutNpyFastaRecord)

  def load_vcf(self, vcf_file):
    assert self.conn, self.c

    # has the table already been loaded?
    md_v = self._c.execute("SELECT v FROM metadata WHERE k = 'vcf_count'").fetchone()[0]

    if int(md_v) != -1:
      raise ValueError("Unexpected DB metadata (vcf_count, %s != -1)" % md_v)

    s_time = time.time()
    processed = 0
    batch = []

    if vcf_file.endswith(".bz2"):
      f_iter = bz2file.BZ2File(vcf_file)
    else:
      f_iter = open(vcf_file)

    for l in f_iter:
      # parse the header, find columns of interest
      if l.startswith("#"):
        if l.startswith("#CHROM"):
          header = l.strip().split()

          kept_cols = []
          cols = ["variant_id"]

          for i, n in enumerate(header):
            if n.startswith("jp"):
              cols.extend(["%s_1" % n, "%s_2" % n, "%s_QV" % n])
              kept_cols.append((i, n))

              self._c.execute("ALTER TABLE vcf ADD COLUMN %s_1 int" % n)
              self._c.execute("ALTER TABLE vcf ADD COLUMN %s_2 int" % n)
              self._c.execute("ALTER TABLE vcf ADD COLUMN %s_QV int" % n)

        continue

      # parse rest of file
      l = l.strip().split()

      chrom = "chr%s" % l[0]
      pos = int(l[1])
      ref = l[3]
      alt = [(str(x + 1), y) for x, y in enumerate(l[4].split(","))]
      code = dict([(str(0), ref)] + alt)

      f_names = l[8].split(":")

      alleles = []
      samples = {}

      for i, n in kept_cols:
        f_data = l[i].split(":")

        col_data = dict(zip(f_names, f_data))
        col_data["GT"] = [code[x] for x in col_data["GT"].split("/", 1)]

        alleles.extend(col_data["GT"])

        if "." in col_data["GT"]:
          break

        samples[n] = col_data
      else:
        # match variant to static db
        v = self.vd.fetch_variant(chrom, pos)

        if v == None:
          # in that case, create a new variant
          v = self.vd.new_variant(chrom, pos)

          # then create new alleles
          for a in set(alleles):
            self.vd.new_allele(v.variant_id, a)

          # and test whether they're non-synonymous
          self.nonsyn_native(variant_id=v.variant_id)
        else:
          # match sure every allele is represented
          v_alleles = [a.sequence for a in v]
          new_alleles = set(alleles).difference(v_alleles)

          # add those that aren't
          for a in new_alleles:
            a_id = self.vd.new_allele(v.variant_id, a)
            self.nonsyn_native(allele_id=a_id)

        # fetch the updated variant
        v = self.vd.fetch_variant(chrom, pos)

        seq_to_allele_id = {}

        for a_seq, a_obj in v.alleles.items():
          seq_to_allele_id[a_seq] = a_obj.allele_id

        row = [v.variant_id]

        for _, n in kept_cols:
          row.append(seq_to_allele_id[samples[n]["GT"][0]],
                     seq_to_allele_id[samples[n]["GT"][1]],
                     samples[n]["GQ"] if "GQ" in samples[n] else "0")

        batch.append(row)
        processed += 1

      if processed % 100 == 0:
        self.c.executemany("INSERT INTO vcf (%s) VALUES (%s)" %
          (", ".join(["'%s'" % x for x in cols]),
           ", ".join(["?" for _ in cols])), batch)
        self.conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.c.executemany("INSERT INTO vcf (%s) VALUES (%s)" %
      (", ".join(["'%s'" % x for x in cols]),
       ", ".join(["?" for _ in cols])), batch)
    self.conn.commit()
    sys.stderr.write("\rprocessed %s @ %.02f/s" %
      (processed, processed / (time.time() - s_time)))

  ##
  ## INHERITANCE FILTER METHODS
  ##

  def mendelian_filter(self, m_filter, chrom, start, end, k = 1):
    pos_in_region = self.fetch_vcf_overlapping(chrom, start, end)

    hits = []

    for pos in itertools.combinations(pos_in_region, k):
      mother, father, son, daughter = {}, {}, {}, {}

      for i in range(k):
        for var, name in zip((mother, father, son, daughter),
                             ("mother", "father", "son", "daughter")):
          try:
            var[i] = [pos[i].samples["%s%s" % (self.stripped_pedigree[name], x)] for x in ("_1", "_2")]
          except KeyError:
            var[i] = None
      else:
        assert son ^ daughter

        result = m_filter(mother, father, son, daughter)

        if result[0]:
          hits.append(result[1])

    return hits

  @staticmethod
  def _mf_dummy(mother, father, son, daughter):
    return (False, None)

  @staticmethod
  def _mf_xlinked(mother, father, son, daughter):
    # at least one reported allele in son must not be in father and be present
    # no more than once in the mother
    for c in son[0]:
      if c not in father[0] and \
        (mother[0][0] != c or \
         mother[0][1] != c):
        return (True, (c,))

    return (False, None)

  @staticmethod
  def _mf_recessive(mother, father, son, daughter):
    child = son if son else daughter

    # homozygous in child, heterozygous in both parents
    if len(set(child[0])) == 1 and \
       len(set(mother[0])) == 2 and \
       len(set(father[0])) == 2:
      return (True, (child[0][0],))

    return (False, None)

  @staticmethod
  def _mf_denovo_dominant(mother, father, son, daughter):
    child = son if son else daughter

    # one or both alleles in child must not be in either parent
    if (child[0][0] not in mother[0] and \
        child[0][0] not in father[0]):
      return (True, (child[0][0],))

    if (child[0][1] not in mother[0] and \
        child[0][1] not in father[0]):
      return (True, (child[0][1],))

    return (False, None)

  @staticmethod
  def _mf_compound_het_denovo(mother, father, son, daughter):
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

    child = son if son else daughter

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

  def non_synonymous_tx(self, tx, pos, seq):
    # fetch exons
    d = self.vd.fetch_tx(tx)

    # non-coding
    if d.cds_start == d.cds_end:
      return []

    coding_exons = []
    middle = False

    # reduce coordinates to those that are coding
    for start, end in zip(d.exon_starts, d.exon_ends):
      if middle:
        # exon contains coding end
        if start <= d.cds_end < end:
          coding_exons.append((start, d.cds_end))
          break

        # everything else should be in the middle of the coding region
        coding_exons.append((start, end))
      else:
        # exon occurs before coding start
        if end < d.cds_start
          continue
        # edge case: exon contains coding start AND end
        elif start <= d.cds_start < end and \
             start <= d.cds_end < end:
          coding_exons.append((d.cds_start, d.cds_end))
          break
        # exon contains coding start
        elif start <= cdsStart < end:
          coding_exons.append((d.cds_start, end))
          middle = True
          continue

    self.pyf_genome[chrom].clearmuts()

    # pull nucleotide sequence of regions
    orig_seq = []

    for start, end in coding_exons:
      orig_seq.append(self.pyf_genome[chrom][start:end])

    # make given mutation
    self.pyf_genome[chrom][pos] = seq

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
      sys.stderr.write("Transcript %s could not be translated\n" % tx)

      diffs = []

    return diffs

  ##
  ## VCF ALLELE FREQUENCY METHODS
  ##

  @staticmethod
  def tabix_af_in_region_native(vcf_file, chrom, start, end):
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
    total_score = []

    for k, v in score_dict.items():
      if not k.startswith("_"):
        total_score.append(float(v) * self.weights[k])

    return sum(total_score)

  ##
  ## NATIVE TESTS
  ##

  def nonsyn_native(self, variant_id = None, allele_id = None):
    if not (variant_id ^ allele_id):
      raise ValueError("Must specify variant_id or allele_id")

    if variant_id:
      v = self.vd.fetch_variant(variant_id=variant_id)
      tx = self.vd.find_tx_overlapping_var(variant_id=variant_id)

      for t in tx:
        for a in v.alleles:
          muts = self.non_synonymous_tx(t.tx_id, v.pos, a.sequence)

          for mut in muts:
            self.vd.new_mut(a.allele_id, t.tx_id, "".join(map(str, mut)))
    elif allele_id:
      a = self.vd.fetch_allele(allele_id=allele_id)
      tx = self.vd.find_tx_overlapping_var(variant_id=a.variant_id)

      for t in tx:
        muts = self.non_synonymous_tx(t.tx_id, v.pos, a.sequence)

        for mut in muts:
          self.vd.new_mut(a.allele_id, t.tx_id, "".join(map(str, mut)))
      
  ##
  ## RUN THE ANALYSIS
  ##

  def score_genes(self):
    # how many recessive, dominant, and comphet models are there
    mendel_counts = {}

    mendel_counts["xlinked"] = self.c.execute(
      "SELECT COUNT(*) FROM mendel WHERE model = 'xlinked'").fetchone()[0]
    mendel_counts["recessive"] = self.c.execute(
      "SELECT COUNT(*) FROM mendel WHERE model = 'recessive'").fetchone()[0]
    mendel_counts["dominant"] = self.c.execute(
      "SELECT COUNT(*) FROM mendel WHERE model = 'dominant'").fetchone()[0]
    mendel_counts["comphet"] = self.c.execute(
      "SELECT COUNT(*) FROM mendel WHERE model = 'comphet'").fetchone()[0]

    # for every gene
    for g in self.vd.fetch_all_genes():
      tx_scores = []

      # take the best-scoring tx for this gene
      for t in g.transcripts:
        tx_score = {}
        tx_size = sum([(y - x) for x, y in zip(t.exon_starts, t.exon_ends)])

        models = self._c.execute("SELECT * FROM mendel WHERE tx_id = %s" % t.tx_id).fetchall()

        if models == None:
          continue

        # take the best-scoring model for this tx
        model_scores = []

        for m in models:
          model_score = {}

          m_a = {0: self.vd.fetch_allele(allele_id=m["allele_id1"])}
          m_v = {0: self.vd.fetch_variant(variant_id=m_a[0].variant_id)}
          m_vcf = {0: self.fetch_vcf(variant_id=m_a[0].variant_id)}

          if m["allele_id2"]:
            m_a[1] = self.vd.fetch_allele(allele_id=m["allele_id2"])
            m_v[1] = self.vd.fetch_variant(variant_id=m_a[1].variant_id)
            m_vcf[1] = self.fetch_vcf(variant_id=m_a[1].variant_id)

          # local af
          for i in m_a:
            m_a[i].local_af = self.local_af[m_v[i].chrom][m_v[i].pos][m_a[i].sequence]

          ni_T = self.exome_size / tx_size
          ni_lambda = mendel_counts[model["model"]] / ni_T
          model_score["mendel"] = self.newell_ikeda(1, ni_lambda, ni_T, 1)

          model_score["nonsyn"] = sum([1 for a in m_a.values() if a.muts]) / float(len(m_a))
          model_score["local_af"] = reduce(lambda x, y: x*y, [a.local_af for a in m_a.values()])
          model_score["global_af"] = reduce(lambda x, y: x*y,[a.evs_af for a in m_a.values()])
          model_score["qv"] = sum([v.samples["%s_QV" % n] for v in m_vcf.values() \
            for n in self.stripped_pedigree.values()]) / (594.0 if m["model"] == "comphet" else 297.0)

          model_scores.append((self.wsm(model_score), model_score))

        best_model_score = sorted(model_scores, key=lambda x: x[0])[-1]
        tx_scores.append(best_model_score)

      if len(tx_scores) == 0:
        best_tx_score = {}
      else:
        best_tx_score = sorted(tx_scores, key=lambda x: x[0])[-1][1]

      best_tx_score["net_cent"] = cent_perc
      best_tx_score["net_nn"] = nn_perc

      gene_score = self.wsm(best_tx_score)

      print gene, gene_score, best_tx_score

  def populate_mendel(self):
    s_time = time.time()
    processed = 0
    batch = []

    for t in self.vd.fetch_all_tx():
      xlinked, comphet, recessive, dominant = [], [], [], []

      if t.chrom == "chrX":
        xlinked = self.mendelian_filter(self._mf_xlinked, t.chrom,
          t.tx_start, t.tx_end)

      comphet = self.mendelian_filter(self._mf_compound_het_denovo, t.chrom,
        t.tx_start, t.tx_end, k=2)
      recessive = self.mendelian_filter(self._mf_recessive, t.chrom, t.tx_start,
        t.tx_end)
      dominant = self.mendelian_filter(self._mf_denovo_dominant, t.chrom,
        t.tx_start, t.tx_end)

      for test, results in zip(("xlinked", "comphet", "recessive", "dominant"), \
                               (xlinked, comphet, recessive, dominant)):
        for result in results:
          allele1 = result[0]

          try:
            allele2 = result[1]
          except IndexError:
            allele2 = None
 
          batch.append((t.tx_id, test, allele1, allele2))

      processed += 1

      if processed % 100 == 0:
        self.c.executemany("INSERT INTO mendel (tx_id, model, allele_id1, " \
          "allele_id2) VALUES (?, ?, ?, ?)", batch)
        self.conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.c.executemany("INSERT INTO mendel (tx_id, model, allele_id1, " \
      "allele_id2) VALUES (?, ?, ?, ?)", batch)
    self.conn.commit()
    sys.stderr.write("\n")
