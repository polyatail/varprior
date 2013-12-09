from multiprocessing import Pool
import numpy
import scipy.stats
import time
import itertools
import os
import sys
import networkx
import sqlite3

# there should be a framework for accessing all this data that doesn't involve
# executing SQL queries manually

# load annotation
# load variants from EVS
# load variants from dbNSFP
# load network
# calculate gene network placement score

class Gene:
  # contains a bunch of transcripts
  def __init__(self, gene_id, gene_name):
    self.gene_id = gene_id
    self.gene_name = gene_name

class Transcript:
  # contains a bunch of proteins (in practice, only one)
  def __init__(self, tx_id, tx_name):
    self.tx_id = tx_id
    self.tx_name = tx_name

class Protein:
  def __init__(self, protein_id, protein_name):
    self.protein_id = protein_id
    self.protein_name = protein_name

class Variant:
  # contains a bunch of alleles
  def __init__(self, variant_id):
    self.variant_id = variant_id

class Allele:
  def __init__(self, allele_id):
    self.allele_id = allele_id

class VariantData:
  def __init__(self, db_file):
    self.version = "variantdata-1.0"
    self.db_file = db_file

  def db_connect(self):
    self._conn = sqlite3.connect(self.db_file)
    self._conn.row_factory = sqlite3.Row
    self._c = self._conn.cursor()

  def db_init(self):
    assert self._conn, self._c

    tables = {}

    tables["metadata"] = {"f": (("k", "text"),
                                ("v", "text"))}

    tables["transcripts"] = {"f": (("tx_id", "integer primary key"),
                                   ("name", "text"),
                                   ("bin", "int"),
                                   ("chrom", "text"),
                                   ("tx_start", "int"),
                                   ("tx_end", "int"),
                                   ("cds_start", "int"),
                                   ("cds_end", "int"),
                                   ("exon_starts", "text"),
                                   ("exon_ends", "text")),
                             "i": ("tx_id", "name", ("bin", "chrom"),
                                  ("tx_start", "tx_end"))}

    tables["proteins"] = {"f": (("protein_id", "integer primary key"),
                                ("name", "text")),
                          "i": ("protein_id", "name")}

    tables["genes"] = {"f": (("gene_id", "integer primary key"),
                             ("name", "text"),
                             ("cent_score", "float"),
                             ("cent_perc", "float"),
                             ("nn_score", "float"),
                             ("nn_perc", "float")),
                       "i": ("gene_id", "name")}

    tables["gene_to_tx"] = {"f": (("gene_id", "int"),
                                  ("tx_id", "int")),
                            "i": ("gene_id", "tx_id")}

    tables["protein_to_tx"] = {"f": (("protein_id", "int"),
                                     ("tx_id", "int")),
                               "i": ("protein_id", "tx_id")}

    tables["tx_to_variant"] = {"f": (("tx_id", "int"),
                                     ("variant_id", "int")),
                               "i": ("tx_id", "variant_id")}

    tables["variants"] = {"f": (("variant_id", "integer primary key"),
                                ("chrom", "text"),
                                ("pos", "int"),
                                ("phastcons", "float"),
                                ("phylop", "float"),
                                ("siphy", "float")),
                          "i": ("variant_id", ("chrom", "pos"))}

    tables["alleles"] = {"f": (("allele_id", "integer primary key"),
                               ("variant_id", "int"),
                               ("sequence", "text"),
                               ("evs_af", "float"),
                               ("tgp_af", "float"),
                               ("polyphen_hdiv", "float"),
                               ("polyphen_hvar", "float"),
                               ("mut_taster", "float"),
                               ("mut_assessor", "float")),
                         "i": ("allele_id", ("variant_id", "sequence"))}

    tables["muts"] = {"f": (("mut_id", "integer primary key"),
                            ("allele_id", "int"),
                            ("tx_id", "int"),
                            ("mut", "text")),
                      "i": ("mut_id", "allele_id", "tx_id")}

    inserts = {}

    inserts["metadata"] = {"f": ("k", "v"),
                           "r": (("version", self.version),
                                 ("transcripts", "-1"),
                                 ("proteins", "-1"),
                                 ("genes", "-1"),
                                 ("variants", "-1"),
                                 ("alleles", "-1"),
                                 ("muts", "-1"))}

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
          "transcripts": ("q", "SELECT COUNT(*) FROM transcripts"),
          "proteins": ("q", "SELECT COUNT(*) FROM proteins"),
          "genes": ("q", "SELECT COUNT(*) FROM genes"),
          "variants": ("q", "SELECT COUNT(*) FROM variants"),
          "alleles": ("q", "SELECT COUNT(*) FROM alleles"),
          "muts": ("q", "SELECT COUNT(*) FROM muts")}

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

  def load_annotation(self, ensgene, enst_to_gene_name, string_alias):
    assert self._conn, self._c

    # has the table already been loaded?
    md_v = self._c.execute("SELECT v FROM metadata WHERE k = 'transcripts'").fetchone()[0]

    if int(md_v) != -1:
      raise ValueError("Unexpected DB metadata (transcripts, %s != -1)" % md_v)

    # can we access all the input files?
    for f in (ensgene, enst_to_gene_nanme, string_alias):
      if not os.path.isfile(f):
        raise ValueError("Cannot access input file %s" % f)

    # load transcripts
    s_time = time.time()
    processed = 0
    batch = []
    db_cols = ", ".join(("name", "tx_start", "tx_end", "exon_starts", "exon_ends", "cds_start", "cds_end"))
    eg_cols = ("name", "txStart", "txEnd", "exonStarts", "exonEnds", "cdsStart", "cdsEnd")
    p_cols = ", ".join(["?" for _ in db_cols])

    for l in open(ensgene_file, "r"):
      l = l.strip().split()

      if processed == 0:
        ensgene_keys = l
      else:
        vals = dict(zip(ensgene_keys, l))
        batch.append([vals[x] for x in eg_cols])

      processed += 1
      if processed % 100 == 0:
        self.c.executemany("INSERT INTO transcripts (%s) VALUES (%s)" % (db_cols, p_cols), batch)
        self.conn.commit()
        batch = []
        sys.stderr.write("\rprocessed %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    self.c.executemany("INSERT INTO transcripts (%s) VALUES (%s)" % (db_cols, p_cols), batch)
    self.conn.commit()
    sys.stderr.write("\n")

    # load ensemblToGeneName (translates ENST -> short names, e.g. IL2RG)
    genes = []
    gene_to_tx = []

    for l in open(enst_to_gene_name):
      l = l.strip().split()

      genes.append(l[1])
      gene_to_tx.append(l[::-1])

    self._c.executemany("INSERT INTO genes (name) VALUES (?)", genes)

    for gene, tx in gene_to_tx:
      gene_id = self._c.execute("SELECT gene_id FROM genes WHERE name = '%s'" % gene).fetchone()[0]
      tx_id = self._c.execute("SELECT tx_id FROM transcripts WHERE name = '%s'" % tx).fetchone()[0]

      self._c.execute("INSERT INTO gene_to_tx (gene_id, tx_id) VALUES (%s, %s)" % (gene_id, tx_id))

    # load string aliases (translates ENSP -> ENST/ENSG)
    proteins = []
    protein_to_tx = []

    for l in open(string_alias):
      l = l.strip().split()

      if l[2].startswith("ENST"):
        proteins.append(l[1])
        protein_to_tx.append(line[1:3])

    self._c.executemany("INSERT INTO proteins (name) VALUES (?)", proteins)

    for protein, tx in protein_to_tx:
      protein_id = self._c.execute("SELECT protein_id FROM proteins WHERE name = '%s'" % protein).fetchone()[0]
      tx_id = self._c.execute("SELECT tx_id FROM transcripts WHERE name = '%s'" % tx).fetchone()[0]

      self._c.execute("INSERT INTO protein_to_tx (protein_id, tx_id) VALUES (%s, %s)" % (protein_id, tx_id))

  def load_evs(self, evs_file):
    assert self._conn, self._c

    evs = {}

    s_time = time.time()
    processed = 0

    for l in open(evs_file):
      if l.startswith("#"):
        continue

      l = l.strip().split()

      chrom, pos = l[0].split(":")
      chrom = "chr%s" % chrom

      #NOTE: here we convert 1-based to 0-based position
      pos = int(pos) - 1

      try:
        evs[chrom]
      except KeyError:
        evs[chrom] = {}

      try:
        evs[chrom][pos]
      except KeyError:
        evs[chrom][pos] = {}

      try:
        evs[chrom][pos]["phastcons"] = float(l[18])
      except ValueError:
        evs[chrom][pos]["phastcons"] = -1

      try:
        evs[chrom][pos]["alleles"]
      except KeyError:
        evs[chrom][pos]["alleles"] = {}

        af = dict([(x, int(y)) for x, y in [x.split("=") for x in l[6].split("/")]])
        total_count = float(sum(af.values()))

        if "R" in af:
          o2n = dict([("A%s" % (x+1), y) for x, y in enumerate([x.split(">")[1] \
            for x in l[3].split(";")])])
          o2n["R"] = l[3].split(">")[0]

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
        variant_id = self._c.execute("SELECT variant_id FROM variants WHERE " \
          "chrom = '%s' AND pos = '%s'" % (chrom, pos)).fetchone()

        if variant_id == None:
          self._c.execute("INSERT INTO variants (chrom, pos, phastcons) VALUES " \
            "('%s', %s, %s)" % (chrom, pos, evs[chrom][pos]["phastcons"]))

          variant_id = self._c.lastrowid
        else:
          variant_id = variant_id[0]

          self._c_.execute("UPDATE variants SET phastcons = '%s' WHERE " \
            "variant_id = '%s'" % (evs[chrom][pos]["phastcons"], variant_id))

        for allele in evs[chrom][pos]["alleles"]:
          allele_id = self._c.execute("SELECT allele_id FROM alleles WHERE " \
            "variant_id = %s AND sequence = '%s'" % (variant_id, allele)).fetchone()

          if allele_id == None:
            self._c.execute("INSERT INTO alleles (variant_id, sequence, evs_af) " \
              "VALUES (%s, '%s', %s)" % (variant_id, allele, evs[chrom][pos]["alleles"][allele]))

            allele_id = self._c.lastrowid
          else:
            allele_id = allele_id[0]

            self._c.execute("UPDATE alleles SET evs_af = %s WHERE allele_id = %s" %
              (allele_id, evs[chrom][pos]["alleles"][allele]))

        processed += 1
        if processed % 100 == 0:
          self.conn.commit()
          sys.stderr.write("\rinserted %s @ %.02f/s" %
            (processed, processed / (time.time() - s_time)))

    self.conn.commit()
    sys.stderr.write("\n")

  def load_dbnsfp(self, dbnsfp_file):
    assert self._conn, self._c

    dbnsfp = {}

    s_time = time.time()
    processed = 0

    for l in open(dbnsfp_file):
      if l.startswith("#"):
        continue

      l = l.strip().split()

      #NOTE: here we convert 1-based to 0-based position
      pos = int(l[1]) - 1
      chrom = "chr%s" % l[0]

      try:
        dbnsfp[chrom]
      except KeyError:
        dbnsfp[chrom] = {}

      try:
        dbnsfp[chrom][pos]
      except KeyError:
        dbnsfp[chrom][pos] = {}

      try:
        dbnsfp[chrom][pos]["phylop"] = float(l[34])
      except ValueError:
        dbnsfp[chrom][pos]["phylop"] = -1

      try:
        dbnsfp[chrom][pos]["siphy"] = float(l[36])
      except ValueError:
        dbnsfp[chrom][pos]["siphy"] = -1

      try:
        dbnsfp[chrom][pos]["alleles"]
      except KeyError:
        dbnsfp[chrom][pos]["alleles"] = {}

      try:
        dbnsfp[chrom][pos]["alleles"][l[3]]
      except KeyError:
        dbnsfp[chrom][pos]["alleles"][l[3]] = {"tx": l[19].split(";"),
                                               "mut": "".join([l[4], l[20], l[5]])}

        try:
          dbnsfp[chrom][pos]["alleles"][l[3]]["tgp_af"] = float(l[40])
        except KeyError:
          dbnsfp[chrom][pos]["alleles"][l[3]]["tgp_af"] = -1

        try:
          dbnsfp[chrom][pos]["alleles"][l[3]]["pp_hdiv"] = float(l[22])
        except KeyError:
          dbnsfp[chrom][pos]["alleles"][l[3]]["pp_hdiv"] = -1

        try:
          dbnsfp[chrom][pos]["alleles"][l[3]]["pp_hvar"] = float(l[24])
        except KeyError:
          dbnsfp[chrom][pos]["alleles"][l[3]]["pp_hvar"] = -1

        try:
          dbnsfp[chrom][pos]["alleles"][l[3]]["mut_taster"] = float(l[28])
        except KeyError:
          dbnsfp[chrom][pos]["alleles"][l[3]]["mut_taster"] = -1

        try:
          dbnsfp[chrom][pos]["alleles"][l[3]]["mut_assessor"] = float(l[30])
        except KeyError:
          dbnsfp[chrom][pos]["alleles"][l[3]]["mut_assessor"] = -1

      processed += 1
      if processed % 100 == 0:
        sys.stderr.write("\rloaded %s @ %.02f/s" %
          (processed, processed / (time.time() - s_time)))

    sys.stderr.write("\n")
    s_time = time.time()
    processed = 0

    for chrom in sorted(dbnsfp.keys()):
      for pos in sorted(dbnsfp[chrom].keys()):
        d = dbnsfp[chrom[pos]]

        variant_id = self._c.execute("SELECT variant_id FROM variants WHERE " \
          "chrom = '%s' AND pos = '%s'" % (chrom, pos)).fetchone()

        if variant_id == None:
          self._c.execute("INSERT INTO variants (chrom, pos, phylop, siphy) VALUES " \
            "('%s', %s, %s)" % (chrom, pos, d["phylop"], d["siphy"]))

          variant_id = self._c.lastrowid
        else:
          variant_id = variant_id[0]

          self._c_.execute("UPDATE variants SET phylop = '%s', siphy = '%s' WHERE " \
            "variant_id = '%s'" % (d["phylop"], d["siphy"], variant_id))

        for allele in d["alleles"]:
          d_a = d[allele]

          allele_id = self._c.execute("SELECT allele_id FROM alleles WHERE " \
            "variant_id = %s AND sequence = '%s'" % (variant_id, allele)).fetchone()

          if allele_id == None:
            self._c.execute("INSERT INTO alleles (variant_id, sequence, tgp_af, " \
              "polyphen_hdiv, polyphen_hvar, mut_taster, mut_assessor) VALUES " \
              "(%s, '%s', %s, %s, %s, %s, %s)" % (variant_id, allele,
              d_a["tgp_af"], d_a["polyphen_hdiv"], d_a["polyphen_hvar"],
              d_a["mut_taster"], d_a["mut_assessor"]))

            allele_id = self._c.lastrowid
          else:
            allele_id = allele_id[0]

            self._c.execute("UPDATE alleles SET tgp_af = '%s', polyphen_hdiv = '%s' " \
              "polyphen_hvar = '%s', mut_taster = '%s', mut_assessor = '%s' " \
              "WHERE allele_id = %s" % (allele_id, d_a["tgp_af"], d_a["polyphen_hdiv"],
              d_a["polyphen_hvar"], d_a["mut_taster"], d_a["mut_assessor"]))

          for tx in d_a["tx"]:
            tx_id = self._c.execute("SELECT tx_id FROM transcripts WHERE name = '%s'" % tx).fetchone()[0]

            self._c.execute("INSERT INTO muts (allele_id, tx_id, mut) VALUES " \
              "(%s, %s, '%s')" % (allele_id, tx_id, d_a["mut"]))

        processed += 1
        if processed % 100 == 0:
          self.conn.commit()
          sys.stderr.write("\rinserted %s @ %.02f/s" %
            (processed, processed / (time.time() - s_time)))

    self.conn.commit()
    sys.stderr.write("\n")

  def find_var_overlapping_tx(self):
    assert self._conn, self._c

    all_tx = self._c.execute("SELECT tx_id, chrom, txStart, txEnd FROM transcripts").fetchall()

    for tx_id, chrom, tx_start, tx_end in all_tx:
      all_vars = self._c.execute("SELECT variant_id FROM variants WHERE " \
        "chrom = '%s' AND pos BETWEEN %s AND %s" % (chrom, tx_start, tx_end)).fetchall()

      batch = [(tx_id, x["variant_id"]) for x in all_vars]

      self._c.executemany("INSERT INTO tx_to_variant (tx_id, variant_id) VALUES (?, ?)", batch)

  def fetch_gene(self, gene_name = None, gene_id = None):
    assert self._conn, self._c

    if not (gene_name ^ gene_id):
      raise ValueError("Must specify either gene_name or gene_id")

    if gene_id:
      d = self._c.execute("SELECT * FROM genes WHERE gene_id = %s" % gene_id).fetchone()
    elif gene_name:
      d = self._c.execute("SELECT * FROM genes WHERE gene_name = %s" % gene_name).fetchone()

    g = Gene(d["gene_id"], d["name"])

    g.cent_score = d["cent_score"]
    g.cent_perc = d["cent_perc"]
    g.nn_score = d["nn_score"]
    g.nn_perc = d["nn_perc"]

    g.transcripts = {}

    all_tx = self._c.execute("SELECT tx_id FROM gene_to_tx WHERE gene_id = '%s'" % g.gene_id).fetchall()

    for tx_id in all_tx:
      t = self.fetch_tx(tx_id=tx_id["tx_id"])

      g.transcripts[t["name"]] = t

  def fetch_tx(self, tx_name = None, tx_id = None):
    assert self._conn, self._c

    if not (tx_name ^ tx_id):
      raise ValueError("Must specify either tx_name or tx_id")

    if tx_id:
      d = self._c.execute("SELECT * FROM transcripts WHERE tx_id = %s" % tx_id).fetchone()
    elif tx_name:
      d = self._c.execute("SELECT * FROM transcripts WHERE tx_name = %s" % tx_name).fetchone()

    t = Transcript(d["tx_id"], d["name"])

    t.bin = d["bin"]
    t.chrom = d["chrom"]
    t.tx_start = d["tx_start"]
    t.tx_end = d["tx_end"]
    t.cds_start = d["cds_start"]
    t.cds_end = d["cds_end"]
    t.exon_starts = d["exon_starts"]
    t.exon_ends = d["exon"]

    t.proteins = {}

    all_prot = self._c.execute("SELECT protein_id FROM protein_to_tx WHERE tx_id = %s" % t.tx_id).fetchall()

    for prot_id in all_prot:
      p = self.fetch_protein(protein_id=prot_id["protein_id"])

      t.proteins[p["name"]] = p

    return t

  def fetch_protein(self, protein_name = None, protein_id = None):
    assert self._conn, self._c

    if not (protein_name ^ protein_id):
      raise ValueError("Must specify either protein_name or protein_id")

    if protein_id:
      d = self._c.execute("SELECT * FROM transcripts WHERE protein_id = %s" % protein_id).fetchone()
    elif protein_name:
      d = self._c.execute("SELECT * FROM transcripts WHERE protein_name = %s" % protein_name).fetchone()

    p = Protein(d["protein_id"], d["name"])

    return p

  def fetch_variant(self, chrom = None, pos = None, variant_id = None):
    assert self._conn, self._c

    if not ((chrom and pos) ^ variant_id):
      raise ValueError("Must specify either chrom and pos or variant_id")

    if variant_id:
      d = self._c.execute("SELECT * FROM variants WHERE variant_id = %s" % variant_id).fetchone()
    elif chrom and pos:
      d = self._c.execute("SELECT * FROM variants WHERE chrom = '%s' AND pos = %s" % (chrom, pos)).fetchone()

    if d == None:
      return None

    v = Variant(d["variant_id"])

    v.chrom = d["chrom"]
    v.pos = d["pos"]
    v.phylop = d["phylop"]
    v.siphy = d["siphy"]

    v.alleles = {}

    all_alleles = self._c.execute("SELECT allele_id FROM alleles WHERE variant_id = '%s'" % v.variant_id).fetchall()

    for allele_id in all_alleles:
      a = self.fetch_allele(allele_id=allele_id["allele_id"])

      v.alleles[a["sequence"]] = a

    return v

  def fetch_allele(self, chrom = None, pos = None, seq = None, allele_id = None):
    assert self._conn, self._c

    if not ((chrom and pos and seq) ^ allele_id):
      raise ValueError("Must specify either chrom, pos and seq or allele_id")

    if allele_id:
      d = self._c.execute("SELECT * FROM alleles WHERE allele_id = %s" % 
        allele_id).fetchone()
    elif chrom and pos and seq:
      d = self._c.execute("SELECT * FROM alleles WHERE chrom = '%s' AND " \
        "pos = %s AND seq = '%s'" % (chrom, pos, seq)).fetchone()

    a = Allele(d["allele_id"])

    a.variant_id = d["variant_id"]
    a.sequence = d["sequence"]
    a.evs_af = d["evs_af"]
    a.tgp_af = d["tgp_af"]
    a.polyphen_hdiv = d["polyphen_hdiv"]
    a.polyphen_hvar = d["polyphen_hvar"]
    a.mut_taster = d["mut_taster"]
    a.mut_assessor = d["mut_assessor"]

    return a

  def gene_names(self):
    assert self._conn, self._c

    genes = self._c.execute("SELECT name FROM genes").fetchall()

    for g in genes:
      yield g["name"]

  def protein_to_gene(self, protein_name):
    assert self._conn, self._c

    genes = self._c.execute("SELECT genes.name AS g FROM genes, proteins, " \
      "protein_to_tx, gene_to_tx WHERE proteins.name = '%s' AND " \
      "proteins.protein_id = protein_to_tx.protein_id AND " \
      "protein_to_tx.tx_id = gene_to_tx.tx_id AND " \
      "genes.gene_id = gene_to_tx.gene_id GROUP BY g").fetchall()
    genes = [x["g"] for x in genes]

    return genes

class GeneNetwork:
  def __init__(self, var_data):
    self.vd = var_data

    self.network_score_cutoff = 677
    self.top_genes = ["ADAMTS8", "AIRE", "AK2", "ATM", "BTK", "CD247", "CD3D",
                      "CD3G", "CD40LG", "CD8A", "CD8B", "CHD7", "CIITA",
                      "CORO1A", "CYBB", "DCLRE1C", "DKC1", "DOCK8", "FOXN1",
                      "IKBKG", "IL2RA", "IL2RG", "IL7R", "ITK", "JAK3", "LCK",
                      "LIG4", "NBN", "NHEJ1", "ORAI1", "PNPLA2", "PRKDC",
                      "PTPRC", "RAG1", "RAG2", "RFXANK", "SH2D1A", "STAT5B",
                      "STIM1", "TAP1", "TAP2", "TAPBP", "TBX1", "WAS", "XIAP",
                      "ZAP70", "ZBTB1"]

  def load_string_pickle(self, gpickle):
    self.gene_graph = networkx.read_gpickle(gpickle)

  def load_string_network(self, string_network_links):
    # NOTE: scores are 1000 - the score so we can use shortest path algorithms
    # NOTE: no multigraph -- instead update to keep edge weights minimum

    # generate STRING graph with gene names
    gene_graph = networkx.Graph(directed=False)
    sys.stderr.write("    nodes\n")
    gene_graph.add_nodes_from(self.vd.gene_names())

    links = [x.strip().split() for x in open(string_network_links)]

    for line_num, line in enumerate(links):
      if float(line[2]) < self.network_score_cutoff:
        continue

      gene1 = self.vd.protein_to_gene(line[0][5:])
      gene2 = self.vd.protein_to_gene(line[1][5:])

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

    self.gene_graph = gene_graph

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

  @staticmethod
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

  def score_all_genes(self, graph):
    # two scores:
    # harmonic centrality defines how well connected a gene is all puck genes
    # nearest-neighobr defines how close a gene is to its nearest puck gene
    p = Pool(NUM_PROCS)
    partial_score_gene = partial(self.score_gene, graph=graph, top_genes=self.top_genes)
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

      batch.append((cent_score, cent_perc, nn_score, nn_perc, gene))

    self.vd._c.executemany("UPDATE genes SET cent_score = ?, cent_perc = ?, " \
      "nn_score = ?, nn_perc = ? WHERE name = ?", batch)
