import base64
import cPickle
import tempfile
import subprocess
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

class VCFEntry(object):
  def __init__(self, vcf_id):
    self.vcf_id = vcf_id

##
## MAIN CLASS
##

class AnalyzeTrio(object):
  def __init__(self, variantdata, db_file):
    self.version = "analyzetrio-1.0"
    self.weights = {"mendel": 1,
                    "nonsyn": 1,
                    "phastcons": 1,
                    "siphy": 1,
                    "phylop": 1,
                    "polyphen_hdiv": 1,
                    "polyphen_hvar": 1,
                    "qv": 1,
                    "global_af": 1,
                    "local_af": 1,
                    "net_cent": 1,
                    "net_nn": 1}
    self.exome_size = 30 * 10 ** 6

    self.vd = variantdata
    self.db_file = db_file

    sys.stderr.write("AnalyzeTrio: connecting to db\n")

    if not os.path.isfile(self.db_file):
      self.db_connect()
    else:
      self.db_connect()
      self.db_metadata("check")

      if bool(self.pedigree):
        self.stripped_pedigree = dict([(x, y.replace("-", "")) for \
                                      x, y in self.pedigree.items()])

      sys.stderr.write("AnalyzeTrio: loading genome\n")
      self.load_genome()

  def load(self, vcf_file, genome_fasta):
    sys.stderr.write("AnalyzeTrio: initializing new db\n")
    self.db_init()

    sys.stderr.write("AnalyzeTrio: loading genome\n")
    self.genome_fasta = genome_fasta
    self.load_genome()

    sys.stderr.write("AnalyzeTrio: loading VCF\n")
    self.vcf_file = vcf_file
    self.load_vcf()

    self.db_metadata("update", update_literal=True)

  def load_pedigree(self, pedigree):
    self.pedigree = pedigree
    self.stripped_pedigree = dict([(x, y.replace("-", "")) for \
                                   x, y in self.pedigree.items()])

    self._c.execute("DELETE FROM mendel")
    self.populate_mendel()
    self.db_metadata("update", update_literal=True)

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
                           ("variant_id", "int"),
                           ("chrom", "text"),
                           ("pos", "int")),
                     "i": ("vcf_id", "variant_id", ("chrom", "pos"))}

    tables["mendel"] = {"f": (("mendel_id", "integer primary key"),
                              ("tx_id", "int"),
                              ("model", "text"),
                              ("allele_id1", "int"),
                              ("allele_id2", "int")),
                        "i": ("mendel_id", "tx_id")}

    inserts = {}

    inserts["metadata"] = {"f": ("k", "v"),
                           "r": (("version", self.version),
                                 ("pedigree", None),
                                 ("genome_fasta", "-1"),
                                 ("vcf_file", "-1"),
                                 ("vcf_count", "-1"),
                                 ("mendel_count", "-1"))}

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

  def db_metadata(self, action, update_literal=False):
    assert self._conn, self._c

    md = {"version": ("l", "version"),
          "genome_fasta": ("u", "genome_fasta"),
          "vcf_file": ("u", "vcf_file"),
          "pedigree": ("p", "pedigree"),
          "vcf_count": ("q", "SELECT MAX(_ROWID_) FROM vcf LIMIT 1"),
          "mendel_count": ("q", "SELECT MAX(_ROWID_) FROM mendel LIMIT 1")}

    if action == "check":
      for k, q in md.items():
        md_v = self._c.execute("SELECT v FROM metadata WHERE k = '%s'" % k).fetchone()[0]

        if q[0] == "u": # update
          self.__dict__[q[1]] = str(md_v)
          continue
        elif q[0] == "p": # update, pickled
          self.__dict__[q[1]] = md_v if md_v == None else cPickle.loads(base64.b64decode(md_v))
          continue
        elif q[0] == "q": # check against query
          ac_v = self._c.execute(q[1]).fetchone()[0]
        elif q[0] == "l": # check against literal
          ac_v = self.__dict__[q[1]]

        if str(md_v) == "-1" and ac_v == None:
          pass
        elif str(md_v) != str(ac_v):
          raise ValueError("Invalid DB metadata (%s, %s != %s)" % (k, md_v, ac_v))
    elif action == "update":
      for k, q in md.items():
        if q[0] == "q":
          ac_v = self._c.execute(q[1]).fetchone()[0]

          if ac_v == None:
            ac_v = -1

          self._c.execute("UPDATE metadata SET v = '%s' WHERE k = '%s'" % (ac_v, k))

        if update_literal == True and \
           q[0] in ("l", "u", "p"):
          try:
            self.__dict__[q[1]]
          except KeyError:
            continue

          if q[0] == "p":
            new_v = base64.b64encode(cPickle.dumps(self.__dict__[q[1]], protocol=cPickle.HIGHEST_PROTOCOL))
          else:
            new_v = self.__dict__[q[1]]

          self._c.execute("UPDATE metadata SET v = '%s' WHERE k = '%s'" % (new_v, k))

      self._conn.commit()

  def fetch_vcf(self, vcf_id = None, variant_id = None):
    assert self._conn, self._c

    if not (bool(vcf_id) ^ bool(variant_id)):
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

    for k in d.keys():
      if k not in ("vcf_id", "variant_id"):
        v.samples[k] = d[k]

    return v

  def fetch_vcf_overlapping(self, chrom, start, end):
    assert self._conn, self._c

    if not (chrom and start and end):
      raise ValueError("Must specify chrom, start, and end")

    d = self._c.execute("SELECT vcf_id FROM vcf WHERE chrom = '%s' " \
      "AND pos BETWEEN %s AND %s" % (chrom, start, end)).fetchall() 

    vcf = [self.fetch_vcf(vcf_id=x["vcf_id"]) for x in d]

    return vcf

  ##
  ## FILE LOADING METHODS
  ##

  def load_genome(self):
    self.pyf_genome = pyfasta.Fasta(self.genome_fasta,
                         record_class=pyfasta.records.MutNpyFastaRecord)

  def load_vcf(self):
    assert self._conn, self._c

    # has the table already been loaded?
    md_v = self._c.execute("SELECT v FROM metadata WHERE k = 'vcf_count'").fetchone()[0]

    if int(md_v) != -1:
      raise ValueError("Unexpected DB metadata (vcf_count, %s != -1)" % md_v)

    s_time = time.time()
    processed = 0
    batch = []

    if self.vcf_file.endswith(".bz2"):
      f_iter = bz2file.BZ2File(self.vcf_file)
    else:
      f_iter = open(self.vcf_file)

    for l in f_iter:
      # parse the header, find columns of interest
      if l.startswith("#"):
        if l.startswith("#CHROM"):
          header = l.strip().split()

          kept_cols = []
          cols = ["variant_id", "chrom", "pos"]

          for i, n in enumerate(header):
            if i > 31:
              n = n.replace("-", "")

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

        if "." in col_data["GT"]:
          break

        col_data["GT"] = [code[x] for x in col_data["GT"].split("/", 1)]

        alleles.extend(col_data["GT"])

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
          v_alleles = [a for a in v.alleles]
          new_alleles = set(alleles).difference(v_alleles)

          # add those that aren't
          for a in new_alleles:
            a_obj = self.vd.new_allele(v.variant_id, a)
            self.nonsyn_native(allele_id=a_obj.allele_id)

        # fetch the updated variant
        v = self.vd.fetch_variant(chrom, pos)

        seq_to_allele_id = {}

        for a_seq, a_obj in v.alleles.items():
          seq_to_allele_id[a_seq] = a_obj.allele_id

        row = [v.variant_id, chrom, pos]

        for _, n in kept_cols:
          row.extend((seq_to_allele_id[samples[n]["GT"][0]],
                      seq_to_allele_id[samples[n]["GT"][1]],
                      samples[n]["GQ"] if "GQ" in samples[n] else "0"))

        batch.append(row)
        processed += 1

      if processed % 100 == 0:
        self._c.executemany("INSERT INTO vcf (%s) VALUES (%s)" %
          (", ".join(["'%s'" % x for x in cols]),
           ", ".join(["?" for _ in cols])), batch)
        self._conn.commit()
        self.vd._conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self._c.executemany("INSERT INTO vcf (%s) VALUES (%s)" %
      (", ".join(["'%s'" % x for x in cols]),
       ", ".join(["?" for _ in cols])), batch)
    self._conn.commit()
    self.vd._conn.commit()

    self.vd.db_metadata("update")

    sys.stderr.write("\n")

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
            pass
      else:
        assert bool(son) ^ bool(daughter)

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
        if end < d.cds_start:
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

  @staticmethod
  def af_from_vcf(vcf_file):
    # allele frequencies from any given vcf
    af_table = {}

    if vcf_file.endswith("bz2"):
      vcftools_in = tempfile.NamedTemporaryFile(delete=False)

      for l in bz2file.BZ2File(vcf_file):
        vcftools_in.write(l)

      vcftools_in.close()
      vcf_file = vcftools_in.name

    vcftools_out = tempfile.NamedTemporaryFile(delete=False)
    
    vcftools = subprocess.Popen(["vcftools", "--vcf", vcf_file,
                                 "--freq", "--out", vcftools_out.name],
                                stdout=open(os.devnull))
    
    while vcftools.poll() == None:
      time.sleep(0.01)

    if os.path.isfile("%s.frq" % vcftools_out.name):
      with open("%s.frq" % vcftools_out.name) as fp:
        fp.next()
        
        for l in fp:
          l = l.strip().split()
        
          chrom = "chr%s" % l[0]
          pos = int(l[1])

          try:
            af_table[chrom]
          except KeyError:
            af_table[chrom] = {}

          try:
            af_table[chrom][pos]
          except KeyError:
            af_table[chrom][pos] = {}

          for allele, freq in [x.split(":") for x in l[4:]]:
            af_table[chrom][pos][allele] = float(freq)

      os.unlink("%s.frq" % vcftools_out.name)

    os.unlink("%s.vcfidx" % vcf_file)
    os.unlink("%s.log" % vcftools_out.name)
    os.unlink(vcftools_out.name)

    try:
      os.unlink(vcftools_in.name)
    except NameError:
      pass

    return af_table

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
    total_score = 1

    for k, v in score_dict.items():
      if k in self.weights:
        total_score *= float(v) ** self.weights[k]

    return total_score

  def wsm(self, score_dict):
    # from http://en.wikipedia.org/wiki/Weighted_sum_model
    total_score = []

    for k, v in score_dict.items():
      if k in self.weights:
        total_score.append(float(v) * self.weights[k])

    return sum(total_score)

  ##
  ## NATIVE TESTS
  ##

  def nonsyn_native(self, variant_id = None, allele_id = None):
    if not (bool(variant_id) ^ bool(allele_id)):
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
    mult = lambda x: reduce(lambda x, y: x*y, x)

    # load local allele frequencies
    self.local_af = self.af_from_vcf(self.vcf_file)

    # how many recessive, dominant, and comphet models are there
#    mendel_counts = {}
#
#    mendel_counts["xlinked"] = self._c.execute(
#      "SELECT COUNT(*) FROM mendel WHERE model = 'xlinked'").fetchone()[0]
#    mendel_counts["recessive"] = self._c.execute(
#      "SELECT COUNT(*) FROM mendel WHERE model = 'recessive'").fetchone()[0]
#    mendel_counts["dominant"] = self._c.execute(
#      "SELECT COUNT(*) FROM mendel WHERE model = 'dominant'").fetchone()[0]
#    mendel_counts["comphet"] = self._c.execute(
#      "SELECT COUNT(*) FROM mendel WHERE model = 'comphet'").fetchone()[0]

    # for every gene
    for g in self.vd.fetch_all_genes():
      tx_scores = []

      # take the best-scoring tx for this gene
      for t in g.transcripts.values():
        tx_score = {}
        tx_size = sum([(int(y) - int(x)) for x, y in zip(t.exon_starts, t.exon_ends)])

        models = self._c.execute("SELECT * FROM mendel WHERE tx_id = %s" % t.tx_id).fetchall()

        if len(models) == 0:
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

#          ni_T = self.exome_size / tx_size
#          ni_lambda = mendel_counts[m["model"]] / ni_T
#          model_score["mendel"] = self.newell_ikeda(1, ni_lambda, ni_T, 1)

          model_score["qv"] = sum([v.samples["%s_QV" % n] for v in m_vcf.values() \
            for n in self.stripped_pedigree.values()]) / (594.0 if m["model"] == "comphet" else 297.0)
          model_score["phastcons"] = mult([v.phastcons if v.phastcons and v.phastcons != -1 else 0 for v in m_v.values()])
          model_score["siphy"] = mult([v.siphy if v.siphy and v.siphy != -1 else 0 for v in m_v.values()])
          model_score["phylop"] = mult([v.phylop if v.phylop and v.phylop != -1 else 0 for v in m_v.values()])

          model_score["polyphen_hdiv"] = sum([a.polyphen_hdiv for a in m_a.values() if a.polyphen_hdiv and a.polyphen_hdiv != -1])
          model_score["polyphen_hvar"] = sum([a.polyphen_hvar for a in m_a.values() if a.polyphen_hvar and a.polyphen_hvar != -1])
          model_score["nonsyn"] = sum([1 for a in m_a.values() if a.muts]) / float(len(m_a))
          model_score["local_af"] = mult([a.local_af for a in m_a.values()])
          model_score["global_af"] = mult([a.evs_af if a.evs_af else 0 for a in m_a.values()])

          #print m_a[0].__dict__, m_v[0].__dict__
          #if model_score["phastcons"] != 0: import pdb; pdb.set_trace()
          model_scores.append((self.wsm(model_score), model_score))

        best_model_score = sorted(model_scores, key=lambda x: x[0])[-1]
        tx_scores.append(best_model_score)

      if len(tx_scores) == 0:
        best_tx_score = {}
      else:
        best_tx_score = sorted(tx_scores, key=lambda x: x[0])[-1][1]

      best_tx_score["net_cent"] = g.cent_perc
      best_tx_score["net_nn"] = g.nn_perc

      gene_score = self.wsm(best_tx_score)

      print g.name, gene_score, best_tx_score

  def populate_mendel(self):
    s_time = time.time()
    processed = 0
    batch = []

    for t in self.vd.fetch_all_tx():
      xlinked, comphet, recessive, dominant = [], [], [], []

      if t.chrom == "chrX" and \
         "son" in self.pedigree:
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
        self._c.executemany("INSERT INTO mendel (tx_id, model, allele_id1, " \
          "allele_id2) VALUES (?, ?, ?, ?)", batch)
        self._conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self._c.executemany("INSERT INTO mendel (tx_id, model, allele_id1, " \
      "allele_id2) VALUES (?, ?, ?, ?)", batch)
    self._conn.commit()
    sys.stderr.write("\n")
