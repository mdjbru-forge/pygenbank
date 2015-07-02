### * Description

# Complementary module/rewriting of genbank.py

### * Setup

### ** Import

import collections
import itertools
import warnings

### * Functions

### ** seqDistance(seq1, seq2)

def seqDistance(seq1, seq2) :
    """Calculate the dissimilarity between two sequences.
    Compare pairs of characters, and count number of mismatches and pairs 
    with at least one X (also considered mismatch).
    Note: matching gaps would be considered matches here! No gaps should be 
    present in the sequences.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence (same length as the first sequence)

    Returns:
        float: Dissimilarity (0-1) between the two sequences
    """
    mismatch = 0
    assert len(seq1) == len(seq2)
    for (x, y) in zip(seq1, seq2) :
        if x != y or "x" in (x + y) or "X" in (x + y) :
            mismatch += 1
    return mismatch * 1. / len(seq1)

### ** seqConsensus(seq1, seq2)

def seqConsensus(seq1, seq2) :
    """Determine the consensus between two sequences, replacing mismatches by X

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence (same length as the first sequence)
    """
    o = ""
    for (x, y) in zip(seq1, seq2) :
        if x == y :
            o += x
        else :
            o += "X"
    return o

### ** mergeSequences(sequences)

def mergeSequences(sequences, maxDistance) :
    """Merge biological sequences of same length based on their similarity

    Args:
        sequences (iterable of str): List or set of strings
        maxDistance (float): Maximum distance allowed to merge sequences

    Returns:
        dictionary: Mapping between the original sequences and the merged 
          sequences
    """
    # https://en.wikipedia.org/wiki/Hierarchical_clustering
    # www.ijetae.com/files/Volume2Issue5/IJETAE_0512_48.pdf (Rambabu 2012, IJETAE,
    # "Clustering Orthologous Protein Sequences thru Python Based Program")
    #
    # Initialization
    sequences = list(sequences)
    clusters = set(sequences)
    assert len(clusters) > 2
    ancestors = dict()
    descendants = dict()
    for c in clusters :
        ancestors[c] = c
    for c in clusters :
        descendants[c] = c
    distances = dict()
    for i in clusters :
        for j in clusters :
            if i != j and not distances.get(frozenset([i, j]), False) :
                distances[frozenset([i, j])] = seqDistance(i, j)
    # Make clusters
    done = False
    while (not done) :
        sortedDistances = sorted(distances.items(), key = lambda x: x[1])
        if sortedDistances[0][1] > maxDistance :
            done = True
        else :
            # Merge
            toMerge = [x[0] for x in sortedDistances if x[1] == sortedDistances[0][1]]
            for pair in toMerge :
                distances.pop(pair)
                pair = list(pair)
                try  :
                    clusters.remove(descendants[pair[0]])
                except :
                    pass
                try :
                    clusters.remove(descendants[pair[1]])
                except :
                    pass
                newCluster = seqConsensus(descendants[pair[0]], descendants[pair[1]])
                clusters.add(newCluster)
                for ancestor in ancestors[pair[0]] :
                    descendants[ancestor] = newCluster
                for ancestor in ancestors[pair[1]] :
                    descendants[ancestor] = newCluster
                descendants[newCluster] = newCluster
                ancestors[newCluster] = ancestors.get(newCluster, [])
                ancestors[newCluster] += list(pair)
            # Refresh distances
            for i in clusters :
                for j in clusters :
                    if i != j and not distances.get(frozenset([i, j]), False) :
                        distances[frozenset([i, j])] = seqDistance(i, j)
            if len(clusters) < 2 :
                done = True
    # Return
    return dict(zip(sequences, [descendants[x] for x in sequences]))
    
### ** extractCodingSeqFast(CDS, seqRecord)

def extractCodingSeqFast(CDS, seqRecord) :
    """Helper function to get the CDS sequence faster than with the extract
    method of a CDS object.
    Note: This will NOT cope with complex locations. It will just extract the 
    nucleotide sequence from the start to the end positions of the CDS location
    and possibly reverse-complement it if needed. Use extractCodingSeqReliable 
    for a safer extraction.
    
    Args:
        CDS (Bio.SeqFeature.SeqFeature): CDS object
        seqRecord (Bio.SeqRecord.SeqRecord): Original record to extract the 
          nucleotide from
    """
    if len(CDS.location.parts) > 1 :
        warnings.warn("CDS with complex location, using "
                      "extractCodingSeqReliable()\n" +
                      str(CDS.location))
        return extractCodingSeqReliable(CDS, seqRecord)
    seq = seqRecord.seq[CDS.location.start:CDS.location.end]
    if CDS.location.strand == -1 :
        seq = seq.reverse_complement()
    return str(seq)

### ** extractCodingSeqReliable(CDS, seqRecord)

def extractCodingSeqReliable(CDS, seqRecord) :
    """Helper function to get the CDS sequence for a CDS object.
    Note: This will cope with complex locations, but is slow. See 
    extractCodingSeqFast for an alternative for simple locations.
    
    Args:
        CDS (Bio.SeqFeature.SeqFeature): CDS object
        seqRecord (Bio.SeqRecord.SeqRecord): Original record to extract the 
          nucleotide from
    """
    return str(CDS.extract(seqRecord).seq)

### * Named tuples

# How to set default values for a named tuple:
# http://stackoverflow.com/questions/11351032/named-tuple-and-optional-keyword-arguments

Gene = collections.namedtuple("Gene", ["recordId", "peptideSeq", "codingSeq",
                                       "peptideLength", 
                                       "translationTable", "gene", "product",
                                       "proteinId", "function", "essentiality",
                                       "peptideHash"])
Gene.__new__.__defaults__ = ("None", ) * 11

Record = collections.namedtuple("Record", ["recordId", "strainName", "database",
                                           "sequence", "organism", "description",
                                           "references"])
Record.__new__.__defaults__ = ("None", ) * 7
                                

### * Classes

### ** ObjectTable()

class ObjectTable(object) :
    """Parent class for more specific table classes"""

### *** __init__(self)

    def __init__(self) :
        self.items = []
        self.itemType = None

### *** __len__(self)

    def __len__(self) :
        return len(self.items)

### *** __getitem__(self, key)

    def __getitem__(self, key) :
        return self.items[key]

### *** __repr__(self)

    def __repr__(self) :
        return ("<ObjectTable (" + str(len(self)) + " items) for " +
                self.itemType.__doc__ + ">")

### *** _oneCol(self, colName)

    def _oneCol(self, colName) :
        """Return an iterator on the "column" with name colName

        Args:
            colName (str): One of the named attributes of the item type

        Returns:
            iterator
        """
        for x in self :
            yield x.__getattribute__(colName)

### *** col(self, *colNames)

    def col(self, *colNames) :
        """Return an iterator on the "columns" with names colNames

        Args:
            colNames (str): One of the named attributes of the item type

        Returns:
            iterator
        """
        # http://stackoverflow.com/questions/243865/how-do-i-merge-two-python-iterators
        return itertools.izip(*[self._oneCol(x) for x in colNames])

### *** loadTable(self, path)

    def loadTable(self, path) :
        """Load gene information from a tabular file. The first line
        contains the headers.

        Args:
            path (str): Path to the file
        """
        with open(path, "r") as fi :
            headers = fi.readline().strip("\n").strip("#").split("\t")
            for line in fi :
                line = line.strip("\n").split()
                data = dict(zip(headers, line))
                self.items.append(self.itemType(**data))
            
### *** writeTable(self, path)

    def writeTable(self, path) :
        """Write gene information to a tabular file

        Args:
            path (str): Path to the file
        """
        with open(path, "w") as fo :
            headers = list(self.itemType()._asdict().keys())
            fo.write("#" + "\t".join(headers) + "\n")
            for item in self.items :
                itemDict = item._asdict()
                fo.write("\t".join([str(itemDict[x]) for x in headers]) + "\n")
        
### ** RecordTable()

class RecordTable(ObjectTable) :
    """Store a table containing record information"""

### *** __init__(self)

    def __init__(self) :
        ObjectTable.__init__(self)
        self.itemType = Record

### *** addGenBankRecord(self, gbRecord)

    def addGenBankRecord(self, gbRecord) :
        """Add a GenBank record data

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object
        """
        d = dict()
        d["recordId"] = "GI:" + gbRecord.annotations["gi"]
        if gbRecord.description.endswith(", complete genome.") :
            d["strainName"] = gbRecord.description[:-18]
        else :
            d["strainName"] = "NA"
        d["database"] = "GenBank"
        d["sequence"] = str(gbRecord.seq)
        d["organism"] = gbRecord.annotations["organism"]
        d["description"] = gbRecord.description
        d["references"] = "<REFSEP>".join([str(x) for x in gbRecord.annotations["references"]]).replace("\n", "<FIELDSEP>")
        self.items.append(Record(**d))

### ** GeneTable()

class GeneTable(ObjectTable) :
    """Store a table containing bacterial gene information"""

### *** __init__(self)

    def __init__(self) :
        ObjectTable.__init__(self)
        self.itemType = Gene

### *** parseGenBankRecord(self, gbRecord)

    def parseGenBankRecord(self, gbRecord) :
        """Parse the content of a GenBank record

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object
        """
        allCDS = [x for x in gbRecord.features if x.type == "CDS"]
        for CDS in allCDS :
            gene = self.itemType(recordId = "GI:" + gbRecord.annotations["gi"],
                        peptideSeq = ";".join(CDS.qualifiers.get("translation", ["None"])),
                        peptideLength = str(len(";".join(CDS.qualifiers.get("translation", ["None"])))),
                        codingSeq = extractCodingSeqFast(CDS, gbRecord),
                        translationTable = ";".join(CDS.qualifiers.get("transl_table", ["None"])),
                        gene = ";".join(CDS.qualifiers.get("gene", ["None"])),
                        product = ";".join(CDS.qualifiers.get("product", ["None"])),
                        proteinId = ";".join(CDS.qualifiers.get("protein_id", ["None"])))
            self.items.append(gene)

### *** hashPeptides(self, hashConstructor)

    def hashPeptides(self, hashConstructor) :
        """Calculate hash value for each gene peptide sequence

        Args:
            hashConstructor (function): Hash algorithm to be used (from the 
              ``hashlib`` module)
        """
        for (i, g) in enumerate(self.items) :
            h = hashConstructor()
            h.update(g.peptideSeq)
            hStr = h.hexdigest()
            geneData = g._asdict()
            geneData["peptideHash"] = hStr
            self.items[i] = self.itemType(**geneData)
            
### *** extractUniquePeptides(self)

    def extractUniquePeptides(self) :
        """Extract the hash and sequences of unique peptides
        """
        uniquePep = []
        uniqueHash = set([])
        for g in self.items :
            assert g.peptideHash is not None
            if g.peptideHash not in uniqueHash :
                uniqueHash.add(g.peptideHash)
                uniquePep.append((g.peptideHash, g.peptideSeq))
        return uniquePep

### *** writeUniquePeptides(self, path)

    def writeUniquePeptides(self, path) :
        """Write the unique peptides to a fasta file

        Args:
            path (str): Path to the fasta file
        """
        uniquePep = self.extractUniquePeptides()
        with open(path, "w") as fo :
            for pep in uniquePep :
                fo.write(">" + pep[0] + "\n")
                fo.write(pep[1] + "\n")

### * Test

from Bio import SeqIO
import os
import hashlib

rootDir = "/home/mabrunea/work/experiments/projects_running/2015-02-05_Ecoli-available-genomes/data/derived/010-fetch-from-genbank/genbank-records"
files = os.listdir(rootDir)
paths = [os.path.join(rootDir, x) for x in files]
n = 5
records = [SeqIO.read(x, "genbank") for x in paths[0:n]]

g = GeneTable()
r = RecordTable()

[g.parseGenBankRecord(x) for x in records]
[r.addGenBankRecord(x) for x in records]

g.writeTable("totoGene")
r.writeTable("totoRecord")
