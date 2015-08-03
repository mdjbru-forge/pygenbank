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

### * Wishlist

script search "e.coli blabla" -o genbank-records/

script tablify genbank-records/* --genes genes.table --records records.table
script updateTable genes.table --hash md5
script mergePeptides genes.table --maxDissimilarity 0.05 --mapping mergedPep.mapping

script blastp mergedPep.mapping --ABC blastp.out.ABC --parallel 4
script clusterMcl blastp.out.ABC --inflation 1.4 -o clusters.groups
script buildAlignments genes.table --clusters clusters.groups --mapping mergedPep.mapping -o clusters/
script refineAlignments

### * Commands hierarchy

- genbank
  + search
  + tablify

- pyGene
  + hash
  + mergePeptides
  
- blastp
  + convertBlastpABC

- cluster
  + convertBlastpABC
  + mcl
  + makeAlignments
  + refineAlignments
  + extendAlignments
  
- SNP
  + call SNPs
