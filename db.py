### * Description

# Complementary module/rewriting of genbank.py

### * Setup

### ** Import

import sys
import collections
import itertools
import warnings
import hashlib
from Bio import SeqIO

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

### ** seqDistances(seqList, distances = None)

def seqDistances(seqList, distances = None, minDistance = None) :
    """Calculate the distances between all pairs of sequences in seqList

    Args:
        seqList (list of str): List of sequences
        distances (dict): A previous output from :func:`seqDistances`, from
          which precalculated distances can be taken
        minDistance (float): Previous minimum distance valuex

    Returns:
        tuple (dict, float): Dictionary mapping frozen sets {s1, s2} to the 
          distance between s1 and s2, where {s1, s2} are all the possible sets 
          of distinct sequences; minimum distance observed, 1 if no 
          pairs of sequences.

    """
    distances_old = distances
    if distances_old is None :
        distances_old = dict()
    if minDistance is None :
        minDistance = 1
    distances = dict()
    for i in seqList :
        for j in seqList :
            if i != j :
                if not distances.get(frozenset([i, j]), False) :
                    if distances_old.get(frozenset([i, j]), False) :
                        distances[frozenset([i, j])] = distances_old[frozenset([i, j])]
                    else :
                        d = seqDistance(i, j)
                        distances[frozenset([i, j])] = d
                        if d < minDistance :
                            minDistance = d
    return (distances, minDistance)

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

### ** seqConsensusList(seqList)

def seqConsensusList(seqList) :
    """Determine the consensus between sequences, replacing mismatches by X

    Args:
        seqList (list of str): Sequences

    Returns:
        str: The consensus sequence

    """
    o = ""
    assert len(seqList) > 0
    length = len(seqList[0])
    assert all([len(x) == length for x in seqList])
    for i in range(length) :
        chars = set([x[i] for x in seqList])
        if len(chars) == 1 :
            o += seqList[0][i]
        else :
            o += "X"
    return o

### ** groupSets(setList)

def groupSets(setList) :
    """Merge sets into larger groups when elements are shared between them.

    Args:
        setList (list of set): List of sets

    Returns:
        list of set: List of sets with sets with shared elements merged
    """
    bag = set([frozenset(x) for x in setList])
    groups = set()
    while len(bag) > 0 :
        previousSize = 0
        element = set(bag.pop())
        while previousSize != len(element) :
            previousSize = len(element)
            overlapping = [x for x in bag if len(element & x) > 0]
            [bag.remove(x) for x in overlapping]
            [element.update(x) for x in overlapping]
        group = frozenset(element)
        groups.add(group)
    return list(groups)

### ** mergeSequences(sequences, maxDistance)

def mergeSequences(sequences, maxDistance) :
    """We can imagine a dendrogram based on the distance between sequences.

    Terminal nodes are the original sequences, and intermediates nodes are
    consensus sequences of their children nodes (parent to the root, children
    towards the leaves) where mismatches between children are replaced by "X".

    One intermediate node can have more than two children.

    We call "simplified sequences" a set of nodes so that every leaf has
    exactly one of these nodes among its parents. At every instant, there is a
    set of current "simplified nodes" and each leaf is assigned to one of them.

    The starting state is with "simplified sequences" being the leaves.

    The target state is with "simplified sequences" being only one or with
    distance between them greater than the maxDistance value. We then return
    the mapping between the leaves and the simplified nodes.

    We proceed iteratively:
    1. Check that there is more than one "simplified node"
    2. Calculate all the distances between the current "simplified nodes"
    3. Determine the minimum distance and check that it is not greater than the
       threshold
    4. Select all the pairs with this minimum distance
    5. Make groups of connected pairs among those
    6. For each independent group of connected pairs, merge them:
         a. Determine the simplified sequence from those simplified sequences.
            This is a newly determined node in the dendrogram.
         b. Map the leaves which were mapped to the simplified sequences of the
            group to this new node
         c. Remove the simplified sequences of the group from the current set 
            of simplified sequences and add the new simplified sequence

    Args:
        sequences (iterable of str): Original sequences
        maxDistance (float): Maximum dissimilarity between sequences to be
          merged

    Returns:
        dict: Dictionary mapping (original sequence, simplified sequence)
    
    """
    leaves = list(sequences)
    simpleNodes = list(sequences)
    mapping = dict(zip(leaves, leaves))
    (distances, minDistance) = seqDistances(simpleNodes)
    stop = (len(simpleNodes) < 2) or (minDistance > maxDistance)
    while not stop :
        bestPairs = [x for x in distances.keys() if distances[x] == minDistance]
        mergingGroups = groupSets(bestPairs)
        for group in mergingGroups :
            simpleNode = seqConsensusList(list(group))
            for leaf in leaves :
                if mapping[leaf] in group :
                    mapping[leaf] = simpleNode
            [simpleNodes.remove(x) for x in group]
            simpleNodes.append(simpleNode)
        (distances, minDistance) = seqDistances(simpleNodes, distances)
        stop = (len(simpleNodes) < 2) or (minDistance > maxDistance)
    return mapping
        
### ** mergeSequencesOld(sequences, maxDistance)

def mergeSequencesOld(sequences, maxDistance) :

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
        
### ** RecordTable() ObjectTable

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

### ** GeneTable() ObjectTable

class GeneTable(ObjectTable) :
    """Store a table containing bacterial gene information"""

### *** __init__(self)

    def __init__(self) :
        ObjectTable.__init__(self)
        self.itemType = Gene
        self.stderr = sys.stderr
        self.nParsedRecords = 0
        self.simplifiedSeqs = None

### *** parseGenBankRecord(self, gbRecord)

    def parseGenBankRecord(self, gbRecord) :
        """Parse the content of a GenBank record

        Args:
            gbRecord (Bio.SeqRecord.SeqRecord): GenBank Record object. It can 
              also be a list of GenBank Record objects, a path to a GenBank
              record or a list of paths to GenBank records.
        """
        # TODO: Set simplifiedSeqs to None when new sequences are parsed?
        if isinstance(gbRecord, list) :
            for r in gbRecord :
                self.parseGenBankRecord(r)
        else :
            self.nParsedRecords += 1
            msg = "Parsing GenBank record " + str(self.nParsedRecords)
            self.stderr.write(msg + "\n")
            if isinstance(gbRecord, str) :
                gbRecord = SeqIO.read(gbRecord, "genbank")
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

### *** extractSimplifiedPeptides(self, maxDissimilarity) :

    def extractSimplifiedPeptides(self, maxDissimilarity) :
        """From the unique peptide sequences, produce simplified sequences which 
        result from merging similar sequences together. Only sequences of same 
        length can be merged, based on their dissimilarity.

        This function returns the mapping and sets the ``simplifiedSeqs`` 
        attribute.

        Args:
            maxDissimilarity (float): Comprised between 0 and 1, maximum 
              dissimilarity for merging

        Returns:
            dict: Dictionary mapping (original peptide, simplified peptide)

        """
        originalSeqs = set([x.peptideSeq for x in self.items])
        originalSeqsByLen = dict()
        for x in originalSeqs :
            l = len(x)
            originalSeqsByLen[l] = originalSeqsByLen.get(l, [])
            originalSeqsByLen[l].append(x)
        mapping = dict()
        for seqs in originalSeqsByLen.values() :
            print("Processing length " + str(len(seqs[0])) +
                  "(" + str(len(seqs)) + " sequences)")
            mapping.update(mergeSequences(seqs, maxDissimilarity))
        self.simplifiedSeqs = mapping
        return mapping
            
### *** writeSimplifiedPeptides(self, path)

    def writeSimplifiedPeptides(self, path) :
        """Write the simplified peptides to a fasta file (i.e. consensus 
        sequences with X at polymorphic positions, obtained by merging 
        sequences of equal length if not too dissimilar). The simplified 
        sequences must have been produced by the
        :func:`extractSimplifiedPeptides` beforehand.

        This is a good input for a blastp run.

        Args:
            path (str): Path to the fasta file

        """
        if self.simplifiedSeqs is None :
            raise Exception("The method to simplify sequences must be called "
                            "first (extractSimplifiedPeptides)")
        with open(path, "w") as fo :
            for pep in set(self.simplifiedSeqs.values()) :
                h = hashlib.md5()
                h.update(pep)
                hStr = h.hexdigest()
                fo.write(">" + hStr + "\n")
                fo.write(pep + "\n")

