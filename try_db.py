### * Test

import os
import hashlib
import db

rootDir = "/home/mabrunea/work/experiments/projects_running/2015-02-05_Ecoli-available-genomes/data/derived/010-fetch-from-genbank/genbank-records"
files = os.listdir(rootDir)
paths = [os.path.join(rootDir, x) for x in files]

g = db.GeneTable()
#g.parseGenBankRecord(paths)
g.loadTable("totoGene")
#[r.addGenBankRecord(x) for x in records]

e = g.extractSimplifiedPeptides(0.10)
g.writeSimplifiedPeptides("totoSimplifiedSeqs.fa")
