### * Description

# Complementary module/rewriting of genbank.py

### * Setup

### ** Import

import collections
import warnings

### * Functions

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
                                       "translationTable", "gene", "product",
                                       "proteinId", "function", "essentiality",
                                       "peptideHash"])
Gene.__new__.__defaults__ = (None, ) * 10

Record = collections.namedtuple("Record", ["recordId", "strainName", "database",
                                           "sequence", "organism", "description",
                                           "references"])
Record.__new__.__defaults__ = (None, ) * 7
                                

### * Classes

### ** RecordTable()

class RecordTable(object) :
    """Store a table containing record information"""

### *** __init__(self)

    def __init__(self) :
        self.records = []

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
        self.records.append(Record(**d))

### ** GeneTable()

class GeneTable(object) :
    """Store a table containing bacterial gene information"""

### *** __init__(self)

    def __init__(self) :
        self.genes = []

### *** parseRecord(self, gbRecord)

    def parseRecord(self, gbRecord) :
        """Parse the content of a GenBank record

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object
        """
        allCDS = [x for x in gbRecord.features if x.type == "CDS"]
        for CDS in allCDS :
            gene = Gene(recordId = "GI:" + gbRecord.annotations["gi"],
                        peptideSeq = ";".join(CDS.qualifiers.get("translation", ["None"])),
                        codingSeq = extractCodingSeqFast(CDS, gbRecord),
                        translationTable = ";".join(CDS.qualifiers.get("transl_table", ["None"])),
                        gene = ";".join(CDS.qualifiers.get("gene", ["None"])),
                        product = ";".join(CDS.qualifiers.get("product", ["None"])),
                        proteinId = ";".join(CDS.qualifiers.get("protein_id", ["None"])))
            self.genes.append(gene)

### *** hashPeptides(self, hashConstructor)

    def hashPeptides(self, hashConstructor) :
        """Calculate hash value for each gene peptide sequence

        Args:
            hashConstructor (function): Hash algorithm to be used (from the 
              ``hashlib`` module)
        """
        for (i, g) in enumerate(self.genes) :
            h = hashConstructor()
            h.update(g.peptideSeq)
            hStr = h.hexdigest()
            geneData = g._asdict()
            geneData["peptideHash"] = hStr
            self.genes[i] = Gene(**geneData)
            
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
                self.genes.add(Gene(**data))
            
### *** writeTable(self, path)

    def writeTable(self, path) :
        """Write gene information to a tabular file

        Args:
            path (str): Path to the file
        """
        with open(path, "w") as fo :
            headers = list(Gene()._asdict().keys())
            fo.write("#" + "\t".join(headers) + "\n")
            for gene in self.genes :
                geneDict = gene._asdict()
                fo.write("\t".join([str(geneDict[x]) for x in headers]) + "\n")

### *** extractUniquePeptides(self)

    def extractUniquePeptides(self) :
        """Extract the hash and sequences of unique peptides
        """
        uniquePep = []
        uniqueHash = set([])
        for g in self.genes :
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
