### * Description

# Complementary module/rewriting of genbank.py

### * Setup

### ** Import

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

### * Classes

### ** GeneTable()

class GeneTable(object) :
    """Store a table containing bacterial gene information"""

### *** __init__(self)

    def __init__(self) :
        self.genes = set([])

### *** parseRecord(self, gbRecord)

    def parseRecord(self, gbRecord) :
        """Parse the content of a GenBank record

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object
        """
        allCDS = [x for x in gbRecord.features if x.type == "CDS"]
        for CDS in allCDS :
            self.genes.add(GeneRow(CDS, gbRecord))

### ** GeneRow()

class GeneRow(object) :
    """Store the information about one gene"""

### *** __init__(self, CDS, seqRecord)

    def __init__(self, CDS, seqRecord) :
        """Args:
            CDS (Bio.SeqFeature.SeqFeature): Biopython feature object containing
              the information about a CDS
            seqRecord (Bio.SeqRecord.SeqRecord): GenBank Record object
        """
        self.recordId = "GI:" + seqRecord.annotations["gi"]
        self.peptideSeq = ";".join(CDS.qualifiers.get("translation", []))
        self.codingSeq = extractCodingSeqFast(CDS, seqRecord)
        self.translationTable = ";".join(CDS.qualifiers.get("transl_table", []))
        self.gene = ";".join(CDS.qualifiers.get("gene", []))
        self.product = ";".join(CDS.qualifiers.get("product", []))
        self.proteinId = ";".join(CDS.qualifiers.get("protein_id", []))
        self.function = ""
        
