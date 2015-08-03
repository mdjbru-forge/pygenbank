### * Description

# Module providing functions and command line scripts to retrieve and process
# data from GenBank

SCRIPT_DESCRIPTION_SEARCH = (""
 "Perform a search on GenBank, retrieve document summaries and download full "
 "GenBank records.")

SCRIPT_DESCRIPTION_EXTRACT_CDS = (""
 "Extract CDS information from GenBank records and produce a summary table. "
 "Can also produce a hash for each CDS and output the unique sequences in "
 "fasta format. If unique sequences in fasta format are required, use the -u "
 "option.")

### * Setup

### ** Import

import os
import time
import argparse
import sys
import xml.etree.ElementTree as ET
import StringIO
import hashlib
import urllib2
import gzip

from Bio import Entrez
from Bio import SeqIO
import pyprind

### ** Default variables

_GB_RECORD_FMTDICT = {
    "id" : lambda x: x.id,
    "name" : lambda x: x.name,
    "dscrp" : lambda x: x.description,
    "acc" : lambda x: ",".join(x.annotations["accessions"]),
    "date" : lambda x: x.annotations["date"],
    "gi" : lambda x: x.annotations["gi"],
    "orgn" : lambda x: x.annotations["organism"],
    "src" : lambda x: x.annotations["source"]
}

def extractNt(c,r,h) :
    # (c, r, h) = (CDS feature, GB record, seq hash)
    seq = r.seq[c.location.start:c.location.end]
    if c.location.strand == -1 :
        seq = seq.reverse_complement()
    a = c.qualifiers.get("translation", ["NA"])[0][1:].strip("*")
    b = str(seq.translate(table = 11))[1:].strip("*")
    if a != b :
        sys.stderr.write(a + "\n")
        sys.stderr.write(b + "\n")
    return str(seq)

_GB_CDS_FMTDICT = {
    "loc" : lambda x: str(x.location),
    "prot" : lambda x: ",".join(x.qualifiers.get("translation", ["NA"])),
    "nuc": extractNt,
    "prod" : lambda x: ",".join(x.qualifiers.get("product", ["NA"])),
    "gene" : lambda x: ",".join(x.qualifiers.get("gene", ["NA"])),
    "hash" : lambda x: x
}

GET_ATTR_FUNCS = {
    # Those are lambda functions taking
    # (c, r, h) = (CDS feature, GB record, seq hash)
    # Record attributes
    "id" : lambda c,r,h: r.id,
    "name" : lambda c,r,h: r.name,
    "dscrp" : lambda c,r,h: r.description,
    "acc" : lambda c,r,h: ",".join(r.annotations["accessions"]),
    "date" : lambda c,r,h: r.annotations["date"],
    "gi" : lambda c,r,h: r.annotations["gi"],
    "orgn" : lambda c,r,h: r.annotations["organism"],
    "src" : lambda c,r,h: r.annotations["source"],
    # CDS attributes
    "loc" : lambda c,r,h: str(c.location),
    "prot" : lambda c,r,h: ",".join(c.qualifiers.get("translation", ["NA"])),
    "nuc": extractNt,
    "prod" : lambda c,r,h: ",".join(c.qualifiers.get("product", ["NA"])),
    "gene" : lambda c,r,h: ",".join(c.qualifiers.get("gene", ["NA"])),
    "hash" : lambda c,r,h: h
}

### * Functions

### ** Related to main_search

### *** search(term, retmax)

def search(term, retmax = None) :
    """Search GenBank for a given query string. Perform a Bio.Entrez.esearch on
    db="nuccore", using history.

    Args:
        term (str): Query to submit to GenBank
        retmax (int): Maximum number of returned ids

    Returns
        Bio.Entrez.Parser.DictionaryElement: The result of the search, with 
        detailed information accessible as if it was a Python dictionary. Keys 
        are::            

            "Count", "RetMax", "IdList", "TranslationStack", "QueryTranslation",
            "TranslationSet", "RetStart", "QueryKey", "WebEnv"
        
        This object can be used with other function of this module to get the 
        actual data (:func:`getDocSum` and :func:`_getFullRecord`)

    """
    handle = Entrez.esearch(db = "nuccore", term = term, retmax =retmax,
                            usehistory = "y")
    results = Entrez.read(handle)
    handle.close()
    return results

### *** getDocSum(searchResult, retmax = None)

def getDocSum(searchResult, retmax = None) :
    """Fetch the documents summaries for the entries from an Entrez.esearch.

    Args:
        searchResult (Bio.Entrez.Parser.DictionaryElement): Object containing
          the result of an Entrez.esearch, typically the output from
          :func: `search` (or at minimum a dictionary with WebEnv 
          and QueryKey entries)
        retmax (int): Maximum number of document summaries to get. If no number
          is given, uses the `RetMax` element from `searchResult`.

    Returns:
        list of dictionaries: A list of dictionaries containing the document 
          summaries

    """
    if retmax is None:
        retmax = searchResult["RetMax"]
    docSumsXML = _getDocSumXML(searchResult = searchResult,
                                     retmax = retmax)
    docSums = _parseDocSumXML(xmlContent = docSumsXML)
    return docSums

### *** getDocSumFromId(listId, retmax = None)

def getDocSumFromId(listId, retmax = None) :
    """Fetch the documents summaries from a list of GenBank identifiers.

    Args:
        listId (list of str): A list of GenBank identifiers
        retmax (int): Maximum number of document summaries to get. If None 
          (default), returns all the document summaries.

    Returns:
        list of dictionaries: A list of dictionaries containing the document 
          summaries

    """
    if retmax is None :
        retmax = len(listId)
    # Epost modified fromt the Biopython cookbook
    mySearch = Entrez.read(Entrez.epost(db = "nuccore", id = ",".join(listId)))
    docSumsXML = _getDocSumXML(searchResult = mySearch,
                                      retmax = retmax)
    docSums = _parseDocSumXML(xmlContent = docSumsXML)
    return docSums    

### *** _getDocSumXML(searchResult, retmax = None)

def _getDocSumXML(searchResult, retmax = None) :
    """Fetch the documents summaries in XML format for the entries from an
    Entrez.esearch.

    Args:
        searchResult (Bio.Entrez.Parser.DictionaryElement): Object containing
          the result of an Entrez.esearch, typically the output from
          :func:`search` (or at minimum a dictionary with WebEnv 
          and QueryKey entries)
        retmax (int): Maximum number of document summaries to get. If no number
          is given, uses the `RetMax` element from `searchResult`.

    Returns:
        str: A string containing the summaries in XML format

    """
    if retmax is None:
        retmax = searchResult["RetMax"]
    handle = Entrez.efetch(db = "nuccore", rettype = "docsum", retmode = "xml",
                           retstart = 0, retmax = retmax,
                           webenv = searchResult["WebEnv"],
                           query_key = searchResult["QueryKey"])
    data = handle.read()
    handle.close()
    return data

### *** _getFullRecord(searchResult, retmax)

def _getFullRecord(searchResult, retmax = 1) :
    """Fetch the full GenBank records for the entries from an Entrez.esearch.

    Args:
        searchResult (Bio.Entrez.Parser.DictionaryElement): Object containing
          the result of an Entrez.esearch, typically the output from
          :func:`search` (or at minimum a dictionary with WebEnv 
          and QueryKey entries).
        retmax (int): Maximum number of full records to get. Default is 1,
          which is a safe approach since individual records can sometimes 
          be very large (e.g. chromosomes).

    Returns:
       str: A string containing the full records in XML format
    
    """
    handle = Entrez.efetch(db = "nuccore", rettype = "gb", retmode = "xml",
                           retstart = 0, retmax = retmax,
                           webenv = searchResult["WebEnv"],
                           query_key = searchResult["QueryKey"])
    data = handle.read()
    handle.close()
    return data

### *** _parseDocSumXML(xmlContent)

def _parseDocSumXML(xmlContent) :
    """Parse the documents summaries from xml format into a list of
    dictionaries.

    Args:
        xmlContent (string): Document summaries in XML format (note: this is a
          string, not a file name). This is typically the output from 
          :func:`_getDocSumXML`.

    Returns:
        list of dictionaries: A list of dictionaries containing the 
          document summaries, or an empty list if no entry was found.

    """
    xmlFile = StringIO.StringIO(xmlContent)
    if "<ERROR>Empty result - nothing to do</ERROR>" in xmlFile.getvalue() :
        return []
    tree = ET.parse(xmlFile)
    root = tree.getroot()
    docsums = []
    for child in root :
        entry = dict()
        nKeys = 0
        for i in child :
            if i.tag == "Item" :
                if i.text == None :
                    i.text = "None"
                entry[i.attrib["Name"]] = i.text
                nKeys += 1
        assert nKeys == 12
        docsums.append(entry)
    return docsums

### *** writeDocSums(docsums, handle)

def writeDocSums(docsums, handle) :
    """Write the documents summaries into a tabular format to any handle with
    a `write` method.

    Args:
        docsums (list of dictionaries or None): A list of dictionaries 
          containing the document summaries, typically the output from 
          :func:`parseDocSumXML`. If it is an empty list, the 
          function will not write anything.
        handle (similar to file handle): Handle object with a `write` method
          (e.g. open file, sys.stdout, StringIO object)

    Returns:
        Nothing but writes the document summaries in a tabular format to the
          specified file.
    """
    if docsums == [] :
        return None
    headers = docsums[0].keys()
    handle.write("#" + "\t".join(headers) + "\n")
    for i in docsums :
        handle.write("\t".join([i[h] for h in headers]) + "\n")

### *** downloadRecords(idList, destDir, batchSize, delay)

def downloadRecords(idList, destDir, batchSize, delay = 30,
                           forceDownload = False,
                           downloadFullWGS = False) :
    """Download the GenBank records for a list of IDs and save them in a
    destination folder. This is the function to use to download data from
    GenBank. It applies a waiting delay between each batch download. Each
    record is saved as a file with name id + ".gb".

    Note that a record is not downloaded if a file with the expected name
    already exists, except if `forceDownload` is `True`.

    The downloading itself is performed by :func:`_downloadBatch` and
    :func:`_getRecordBatch`.
    
    some GenBank records do not contain actual sequence data but some reference
    to a WGS (whole genome shotgun sequencing) project. For those, setting
    `downloadFullWGS` to True is necessary to download another GenBank file
    with the actual sequence data. Note that if a GenBank record was first
    downloaded without this option, and actually contains a WGS reference, then
    the `forceDownload` option must be enabled (or the file must be removed)
    for the WGS file to be also downloaded in a new call of this function.

    Args:
        idList (list of str): List of GenBank id
        destDir (str): Path to the folder where the records will be saved
        forceDownload (boolean): Should records for which a destination file
          already exists be downloaded anyway? (default: `False`)
        downloadFullWGS (boolean): If True, also download the full GenBank 
          files corresponding to GenBank records with WGS trace reference
    
    Returns:
        Nothing, but saves the GenBank records in the destination folder
        specified by `destDir`.

    """
    # Create the destination directory if it doesn't exist
    if not os.path.isdir(destDir) :
        os.makedirs(destDir)
    # Filter the list for only records not already downloaded
    existingFileList = os.listdir(destDir)
    if forceDownload :
        newIdList = idList
    else :
        newIdList = [x for x in idList if not ((x + ".gb") in existingFileList)]

    # Download the batches
    bar = pyprind.ProgBar(len(newIdList) / batchSize, monitor=True,
                          title='Downloading the batches')
    for i in range(0, len(newIdList), batchSize):
        end = min(len(newIdList), i + batchSize)
        batch = newIdList[i: end]
        _downloadBatch(batch, destDir, downloadFullWGS)
        time.sleep(delay)
        bar.update()

### *** _downloadBatch(idBatch, destDir, downloadFullWGS = False)

def _downloadBatch(idBatch, destDir, downloadFullWGS = False) :
    """Download a batch of GenBank records to a destination directory. You should
    not call this function directly, but rather use
    :func:`downloadRecords` (which itself calls
    :func:`_downloadBatch`) for your downloads.

    :func:`_downloadBatch` calls :func:`_getRecordBatch` to
    download data from GenBank, and then takes care of separating individual
    records and writing them to files.

    Args:
        idBatch (list of str): List of GenBank id
        destDir (str): Path to the folder where the records will be saved
        downloadFullWGS (boolean): If True, also download the full GenBank 
          files corresponding to GenBank records with WGS trace reference

    Returns:
        Nothing, but saves the GenBank records in the destination folder
        specified by `destDir`.

    """
    # Download the records
    r = _getRecordBatch(idBatch)
    # Split the downloaded strings into records
    r = r.strip().split("\n//")
    # print(r)
    assert r[-1] == ""
    r = r[0: -1]
    assert len(r) == len(idBatch)
    # Save the records
    for i in range(len(idBatch)) :
        with open(os.path.join(destDir, idBatch[i] + ".gb"), "w") as fo :
            fo.write(r[i].strip())
            fo.write("\n//\n")
        WGSline = _recordIsWGS(r[i])
        if WGSline and downloadFullWGS :
            downloadWGS(r[i], destDir)

### *** _recordIsWGS(recordStr)

def _recordIsWGS(recordStr) :
    """Check if a GenBank record is from a whole genome shotgun project. This
    is done by searching for the "WGS " string at the beginning of a line

    Args:
        recordStr (str): Content of a GenBank record

    Returns:
        str or False: WGS line or False

    """
    lines = recordStr.split("\n")
    WGS_lines = [x for x in lines if x.startswith("WGS ")]
    if len(WGS_lines) == 1 :
        return WGS_lines[0]
    elif len(WGS_lines) == 0 :
        return False
    else :
        raise Exception("Several lines starting with \"WGS \" in a GenBank record")        

### *** _makeWGSurl(WGSline)

def _makeWGSurl(WGSline) :
    """Prepare the url to download the GenBank records corresponding to one
    WGS line. The WGS line is the output from 
    :func:`_genbankRecirdIsWGS`.

    Args:
        WGSline (str): WGS line from a GenBank record, output from 
          :func:`_genbankRecirdIsWGS`.

    Returns:
        str: Url to download the record

    """
    if not WGSline.startswith("WGS ") :
        raise Exception("Line does not start with \"WGS \"")
    accession = WGSline.split(" ")[-1]
    accRoot = accession.split("-")[0][0:6]
    url = "http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=" + accRoot + ".1.gbff.gz"
    return url

### *** _downloadWGS(WGSurl)

def _downloadWGS(WGSurl) :
    """Download a WGS GenBank file. The output is an uncompressed version of the
    file.

    Args:
        WGSurl (str): Url to download the gzip file, output from 
          :func:`_makeWGSurl`

    Returns:
        str: Uncompressed GenBank file content. Return None if there was a 
          problem during gzip decompression.

    """
    gzipContent = urllib2.urlopen(WGSurl).read()
    gzipFile = StringIO.StringIO(gzipContent)
    o = gzip.GzipFile(fileobj = gzipFile)
    output = None
    try :
        output = o.read()
    except IOError as e:
        print(e)
    o.close()
    return output

### *** downloadWGS(gbRecord, destDir)

def downloadWGS(gbRecord, destDir) :
    """Download and save the WGS GenBank file corresponding to a GenBank
    record with a WGS reference

    Args:
        gbRecord (str): Text content of a GenBank record with WGS reference
        destDir (str): Path to the directory where to save the GenBank file

    Returns:
        Nothing, but save the complete GenBank file corresponding to the WGS
          reference into the specified folder. The file name is the GI number
          plus "WGS"

    """
    WGS = _recordIsWGS(gbRecord)
    assert WGS
    url = _makeWGSurl(WGS)
    WGS_content = _downloadWGS(url)
    gb = SeqIO.read(StringIO.StringIO(gbRecord + "\n//"), "genbank")
    gi = gb.annotations["gi"]
    filePath = os.path.join(destDir, gi + "_WGS.gb")
    with open(filePath, "w") as fo :
        fo.write(WGS_content)
    return

### *** _getRecordBatch(idList)

def _getRecordBatch(idList) :

    """Retrieve the GenBank records for a list of GenBank id (GIs). This is a
    relatively low-level function that only gets data from GenBank but does not
    manage batches or write files. You should use the higher level wrapper
    :func:`downloadRecords` for your own downloads.

    Args:
        idBatch (list of str): List of GenBank id

    Returns:
        str: A string with all the data (can be large if many id are provided)

    """
    handle = Entrez.efetch(db = "nuccore", rettype = "gbwithparts",
                           retmode = "text", id = ",".join(idList))
    r = handle.read()
    handle.close()
    return r

### *** _fileLinesToList(filename)

def _fileLinesToList(filename) :
    """Simple function to get a list of stripped lines from a file.

    Args:
        filename (str): Path to the file to read

    Returns:
        list of str: A list containing all non-empty white-stripped lines from
          the file

    """
    o = []
    with open(filename, "r") as fi :
        for l in fi :
            if l.strip() != "" :
                o.append(l.strip())
    return o

### ** Related to main_extract_CDS

### *** _summarizeRecord(record, summaryFormat, hashConstructor, existingHashes)

def _summarizeRecord(record, summaryFormat, hashConstructor,
                            existingHashes = dict()):
    """Produce a tabular summary of all the CDS features present in a GenBank
    record and a dictionary containing hashes of the unique sequences (hash,
    protein sequence). If a dictionary of pre-existing hashes is given, update
    this one. Checks for collisions in the hash dictionary.

    Args:
        record (Bio.SeqRecord.SeqRecord): GenBank record (Biopython object)
        summaryFormat (list of str): List of attribute descriptors determining 
          the columns of the summary table
        hashConstructor (function): Hash algorithm to be used (from the 
          ``hashlib`` module)
        existingHashes (dict): Dictionary (k, v) = (hash, protein sequences).
          This is updated with the CDS hashes from the input `record` and 
          checked for collisions.

    Returns:
        (str, dict, dict): A tuple containing the following objects

          * string: tabular summary for all CDS features in the GenBank record,
            ready to be written to an output stream

          * dictionary (hash, protein sequences): dictionary given in input as
            `existingHashes` and updated with the hashes from `record`. This is
            the dictionary one can pass to another call to
            :func:`_summarizeRecord` in order to progressively build a
            complete dictionary of all hashes for several GenBank records.

          * dictionary (hash, protein sequences): dictionary containing only
            the new hashes not already present in `existingHashes`. This is
            useful if one wants to update an output stream with the unique
            hashes after each call to :func:`_summarizeRecord` when
            processing several records.

    """
    list_CDS = [x for x in record.features if x.type == "CDS"]
    summary = ""
    newHashes = dict()
    currentHashes = existingHashes.copy()
    currentHashes_keys = set(currentHashes.keys())
    for CDS in list_CDS :
        (protSeq, h) = _getProteinHashFromCDS(CDS, hashConstructor)
        summary += (_makeSummaryForCDS(record, CDS, h, summaryFormat) +
                    "\n")
        if h in currentHashes_keys :
            assert protSeq == currentHashes[h]
        else :
            newHashes[h] = protSeq
            currentHashes[h] = protSeq
            currentHashes_keys.add(h)
    return (summary, currentHashes, newHashes)

### *** _getProteinHashFromCDS(CDS, hashConstructor)

def _getProteinHashFromCDS(CDS, hashConstructor) :
    """Extract the protein sequence from a Bio.SeqFeature.SeqFeature CDS object and
    determine its hash value.

    Args:
        CDS (Bio.SeqFeature.SeqFeature): CDS of interest
        hashConstructor (function): Hash algorithm to be used (from the 
          ``hashlib`` module)

    Returns:
        (str, str): A tuple containing the protein sequence and the hash value. 
          If there is no translation available in the CDS qualifiers, returns
          "NA" as the protein sequence. If there are several translation 
          available, join them with commas.

    """
    protSeq = ",".join(CDS.qualifiers.get("translation", ["NA"])) 
    h = hashConstructor()
    h.update(protSeq)
    hStr = h.hexdigest()
    return (protSeq, hStr)

### *** _makeSummaryForCDS(record, CDS, hStr, summaryFormat)

def _makeSummaryForCDS(record, CDS, hStr, summaryFormat, getAttrFuncs = None) :

    """Make a summary for one CDS feature object.

    Args:
        record (Bio.SeqRecord.SeqRecord): Parent record
        CDS (Bio.SeqFeature.SeqFeature): CDS of interest
        hStr (str): Hash string for the protein sequence of the CDS
        summaryFormat (list of str): List of attribute descriptors determining 
          the columns of the summary table
        getAttrFuncs (dict of functions): Dictionary mapping the attribute
          descriptors to functions of the form: `f(CDS, record, hStr)`. If None
          (default), use the module-defined GET_ATTR_FUNCS dictionary.

    Returns:
        str: Summary string for this CDS

    """
    if getAttrFuncs is None :
        getAttrFuncs = GET_ATTR_FUNCS
    summaryElements = [getAttrFuncs[x](CDS, record, hStr) for x in summaryFormat]
    return "\t".join(summaryElements)

### ** Not related to a main function

### *** _recordInfo(record, outfmt, fmtdict)

def _recordInfo(record, outfmt, fmtdict = None) :
    """Get some information about a GenBank record (stored as a
    Bio.SeqRecord.SeqRecord object).

    Args:
        record (Bio.SeqRecord.SeqRecord): GenBank record
        outfmt (list of str): List of information keys
        fmtdict (dict): Dictionary mapping the information keys to simple 
          functions to retrieve the corresponding information. If None, use
          the default _GB_RECORD_FMTDICT

    Returns:
        List: List containing the information corresponding to each key in 
          `outfmt`

    """
    if fmtdict is None :
        fmtdict = _GB_RECORD_FMTDICT
    return [fmtdict[x](record) for x in outfmt]

### *** _CDSinfo(CDS, outfmt, fmtdict)

def _CDSinfo(CDS, outfmt, fmtdictCDS = None,
                    fmtdictRecord = None,
                    parentRecord = None, hashConstructor = None) :
    """Get some information about a GenBank CDS (stored as a
    Bio.SeqFeature.SeqFeature object of type "CDS").

    Args:
        CDS (Bio.SeqFeature.SeqFeature of type "CDS"): GenBank CDS
        outfmt (list of str): List of information keys
        fmtdictCDS (dict): Dictionary mapping the information keys to simple 
          functions to retrieve the corresponding information (for CDS 
          attributes). If None, default is _GB_CDS_FMTDICT.
        fmtdictRecord (dict): Dictionary mapping the information keys to simple 
          functions to retrieve the corresponding information (for record 
          attributes). If None, default is _GB_RECORD_FMTDICT.
        parentRecord (Bio.SeqRecord.SeqRecord): Parent record, needed to extract
          nucleotide sequence and other record-related information
        hashConstructor (function): Hash algorithm to be used (from the 
          ``hashlib`` module)

    Returns:
        dict: Dictionary containing the information corresponding to each key 
          in `outfmt`

    """
    if fmtdictCDS is None :
        fmtdictCDS = _GB_CDS_FMTDICT
    if fmtdictRecord is None :
        fmtdictRecord = _GB_RECORD_FMTDICT
    CDS.parentRecord = parentRecord
    info = dict()
    for k in [x for x in outfmt if x in fmtdictCDS.keys()] :
        info[k] = fmtdictCDS[k](CDS)
    for k in [x for x in outfmt if x in fmtdictRecord.keys()] :
        info[k] = fmtdictRecord[k](parentRecord)
    if hashConstructor is not None :
        protSeq = fmtdictCDS["prot"](CDS)
        h = hashConstructor()
        h.update(protSeq)
        info["hash"] = h.hexdigest()
    return info

### * Main scripts functions

### ** pygenbank-search

# Wishlist
#
# pygenbank-search --query "blabla" >> docsum.table
# pygenbank-search --query "blabla" --download
# pygenbank-search --idlist myId.list >> docsum.table
# pygenbank-search --idlist myId.list --download
# pygenbank-search --idlist myId.list --download --outputDir ./GenBank

### *** main_parser

def main_parser() :
    """Prepare the parser

    Returns:
        ArgumentParser: An argument parser

    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help = "")
    # Search
    sp_search = subparsers.add_parser("search",
                                      help = "Search GenBank",
                                      description = SCRIPT_DESCRIPTION_SEARCH)
    sp_search.add_argument("-c", "--count", action = "store_true",
                           help = "Just return the number of records, no fetch")
    # Required named arguments (http://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments)
    required = sp_search.add_argument_group("required named arguments")
    # --email
    required.add_argument("-e", "--email", type = str,
                        help = "User's email (required by Entrez)")
    # --listId
    required.add_argument("-l", "--listId", type = str,
                          help = "File containing one GenBank identifier per "
                          "line. Use - for reading from stdin. "
                          "Exactly one of --listId or "
                          "--query must be specified, but not both.")
    # --query
    required.add_argument("-q", "--query", type = str,
                          help = "Query string for GenBank search. "
                          "Exactly one of --listId or "
                          "--query must be specified, but not both.",
                          metavar = "SEARCH_TERM")
    # Download options
    download = sp_search.add_argument_group("download-related options")
    # --retmax
    download.add_argument("-r", "--retmax", type = int, default = 0,
                        help = "Maximum number of entries to retrieve from "
                        "GenBank, comprised between 1 and 10000. Use 0 for "
                        "unlimited number of returned entries. (default: 0)")
    # --download
    download.add_argument("-d", "--download", action = "store_true",
                        help = "Download the full GenBank records")
    # --forceDownload
    download.add_argument("-f", "--forceDownload", action = "store_true",
                        help = "Download record even if file already exists "
                        "(implies --download)")
    # --fullWGS
    download.add_argument("--fullWGS", action = "store_true",
                        help = "Also download full WGS sequence data when "
                        "WGS trace reference is present in a GenBank record "
                        "(only works if the original GenBank record is to be "
                        "downloaded too or if --forceDownload is used)")
    # --outputDir
    download.add_argument("-o", "--outputDir", type = str, default = ".",
                        help = "Destination folder for downloaded records "
                        "(default: current directory)")
    # --batchSize
    download.add_argument("-b", "--batchSize", type = int, default = 5,
                        help = "Batch size for full record retrieval "
                        "(default: 5)")
    # --delay
    download.add_argument("--delay", type = int, default = 15,
                        help = "Delay in seconds between successive batch "
                        "retrieval of the full records (default: 15)")
    # Action
    sp_search.set_defaults(action = "search")
    # Extract CDS
    sp_extractCDS = subparsers.add_parser("extractCDS",
                                          help = "Extract CDS",
                                          description = SCRIPT_DESCRIPTION_EXTRACT_CDS)
        # input files
    sp_extractCDS.add_argument("genbank_records", type = str, nargs = "+",
                               help = "Filename(s) for GenBank record(s)")
    # --hash
    sp_extractCDS.add_argument("--hash", metavar = "HASH_ALGORITHM",
                        choices = ["md5", "sha1", "sha224", "sha256", "sha384",
                                   "sha512"],
                        default = "md5",
                        help = "Hash algorithm to use for unique sequence signature "
                               "(default md5)")
    # --unique    
    sp_extractCDS.add_argument("-u", "--unique", metavar = "UNIQUE_FASTA",
                        type = str, default = None,
                        help = "Filename for unique fasta amino acid sequences")
    # --outfmt
    sp_extractCDS.add_argument("-o", "--outfmt", metavar = "FORMAT",
                        type = str, default = "dscrp,gi,loc,hash,prot",
                        help = "Comma-separated list of output column "
                        "identifiers for summaries. Allowed identifiers are: "
                        "id, name, dscrp, acc, date, gi, orgn, src (record "
                        "attributes) and loc, prot, nuc, prod, gene, hash "
                        "(CDS attributes). Default is dscrp,gi,loc,hash,prot.")
    # --count
    sp_extractCDS.add_argument("-c", "--count", action = "store_true",
                        help = "Count the number of CDS per record instead of "
                        " producing CDS summaries. Cannot be used at the same "
                        "time as --unique.")
    # --verbose
    sp_extractCDS.add_argument("-v", "--verbose", action = "store_true",
                        help = "Send verbose messages to stderr during execution")
    # Action
    sp_extractCDS.set_defaults(action = "extractCDS")
    return parser

### *** main

def main(args = None, stdout = None, stderr = None) :
    """Main entry point

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    """
    if args is None :
        parser = main_parser()
        args = parser.parse_args()
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    dispatch = dict()
    dispatch["search"] = main_search
    dispatch["extractCDS"] = main_extractCDS
    dispatch[args.action](args, stdout, stderr)

### *** main_search

def main_search(args, stdout, stderr) :
    args = _processArgsToLogic_search(args, stdout, stderr)
    listId = None
    # Genbank search
    if args.actionFlags.get("DoGenbankSearch", False) :
        mySearch = search(term = args.query, retmax = args.retmax)
        if args.count :
            stdout.write(mySearch["QueryTranslation"] + "\t" + str(mySearch["Count"]) + "\n")
            sys.exit(0)
        myDocSums = getDocSum(mySearch)
        writeDocSums(myDocSums, stdout)
        listId = [x["Gi"] for x in myDocSums]
    # Get docsums for a list of identifiers
    if args.actionFlags.get("DoGetList", False) :
        if args.count :
            stderr.write("-l and -c cannot be used at the same time\n")
            sys.exit(1)
        listId = _fileLinesToList(args.listId)
        myDocSums = getDocSumFromId(listId)
        writeDocSums(myDocSums, stdout)
    # Download records
    if args.download and not args.count :
        assert listId is not None
        downloadRecords(idList = listId, destDir = args.outputDir,
                               batchSize = args.batchSize, delay = args.delay,
                               forceDownload = args.forceDownload,
                               downloadFullWGS = args.fullWGS)

### *** main_extractCDS

def main_extractCDS(args, stdout, stderr) :
    args = _processArgsToLogic_extract_CDS(args, stdout, stderr,
                                           GB_RECORD_FMTDICT, GB_CDS_FMTDICT)
    # Go through the input files
    uniqueSeq = dict()
    i_file = 0
    for fi in args.genbank_records :
        i_file += 1
        if args.verbose :
            stderr.write(time.asctime() + " - " +
                         "Processing file " + str(i_file) + " : " +
                         os.path.basename(fi) + " - " +
                         "N unique seq : " + str(len(uniqueSeq.keys())) + "\n")
        record = SeqIO.parse(fi, "genbank")
        for r in record :
            if not args.actionFlags.get("DoCount", False) :
                (summaryString, uniqueSeq, newSeq) = (
                    _summarizeRecord(r, args.outfmt, args.hash, uniqueSeq))
                stdout.write(summaryString)
            else :
                count = len([x for x in r.features if x.type == "CDS"])
                stdout.write(r.annotations["gi"] + "\t" + str(count) + "\n")
    # Write unique sequences
    if args.actionFlags.get("DoUniqueSequences", False) :
        with open(args.unique, "w") as fo :
            for (k, v) in uniqueSeq.items() :
                fo.write(">" + k + "\n")
                fo.write(v + "\n")
        
### *** _main_search(args = None, stdout = None, stderr = None)

def _main_search(args = None, stdout = None, stderr = None) :

    """Main function, used by the command line script "-search". This function
    sends the arguments to :func:`_processArgsToLogic_search` to determine
    which actions must be performed, and then performs the actions.

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    Returns:
        None

    """
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    # Process arguments
    if args is None :
        parser = _makeParser_search()
        args = parser.parse_args()
    args = _processArgsToLogic_search(args, stdout, stderr)
    listId = None
    # Genbank search
    if args.actionFlags.get("DoGenbankSearch", False) :
        mySearch = search(term = args.query, retmax = args.retmax)
        if args.count :
            stdout.write(mySearch["QueryTranslation"] + "\t" + str(mySearch["Count"]) + "\n")
            sys.exit(0)
        myDocSums = getDocSum(mySearch)
        writeDocSums(myDocSums, stdout)
        listId = [x["Gi"] for x in myDocSums]
    # Get docsums for a list of identifiers
    if args.actionFlags.get("DoGetList", False) :
        if args.count :
            stderr.write("-l and -c cannot be used at the same time\n")
            sys.exit(1)
        listId = _fileLinesToList(args.listId)
        myDocSums = getDocSumFromId(listId)
        writeDocSums(myDocSums, stdout)
    # Download records
    if args.download and not args.count :
        assert listId is not None
        downloadRecords(idList = listId, destDir = args.outputDir,
                               batchSize = args.batchSize, delay = args.delay,
                               forceDownload = args.forceDownload,
                               downloadFullWGS = args.fullWGS)

### *** _makeParser_search()

def _makeParser_search() :
    """Build the argument parser for the main script "-search".
    
    Returns:
        argparse.ArgumentParser() object: An argument parser object ready to be
        used to parse the command line arguments
    """
    parser = argparse.ArgumentParser(
        description = SCRIPT_DESCRIPTION_SEARCH)
    parser.add_argument("-c", "--count", action = "store_true",
                        help = "Just return the number of records, no fetch")
    # Required named arguments (http://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments)
    required = parser.add_argument_group("required named arguments")
    # --email
    required.add_argument("-e", "--email", type = str,
                        help = "User's email (required by Entrez)")
    # --listId
    required.add_argument("-l", "--listId", type = str,
                          help = "File containing one GenBank identifier per "
                          "line. Use - for reading from stdin. "
                          "Exactly one of --listId or "
                          "--query must be specified, but not both.")
    # --query
    required.add_argument("-q", "--query", type = str,
                          help = "Query string for GenBank search. "
                          "Exactly one of --listId or "
                          "--query must be specified, but not both.",
                          metavar = "SEARCH_TERM")
    # Download options
    download = parser.add_argument_group("download-related options")
    # --retmax
    download.add_argument("-r", "--retmax", type = int, default = 0,
                        help = "Maximum number of entries to retrieve from "
                        "GenBank, comprised between 1 and 10000. Use 0 for "
                        "unlimited number of returned entries. (default: 0)")
    # --download
    download.add_argument("-d", "--download", action = "store_true",
                        help = "Download the full GenBank records")
    # --forceDownload
    download.add_argument("-f", "--forceDownload", action = "store_true",
                        help = "Download record even if file already exists "
                        "(implies --download)")
    # --fullWGS
    download.add_argument("--fullWGS", action = "store_true",
                        help = "Also download full WGS sequence data when "
                        "WGS trace reference is present in a GenBank record "
                        "(only works if the original GenBank record is to be "
                        "downloaded too or if --forceDownload is used)")
    # --outputDir
    download.add_argument("-o", "--outputDir", type = str, default = ".",
                        help = "Destination folder for downloaded records "
                        "(default: current directory)")
    # --batchSize
    download.add_argument("-b", "--batchSize", type = int, default = 5,
                        help = "Batch size for full record retrieval "
                        "(default: 5)")
    # --delay
    download.add_argument("--delay", type = int, default = 15,
                        help = "Delay in seconds between successive batch "
                        "retrieval of the full records (default: 15)")
    return parser

### *** _processArgsToLogic_search(args, stdout, stderr)

def _processArgsToLogic_search(args, stdout, stderr) :
    """Process the command line arguments and determine the action logic for the
    :func:`_main_search` function.

    Args:
        args (namespace): Argument namespace
        stdout (file): stdout stream
        stderr (file): stderr stream

    Returns:
        namespace: The argument namespace given input in `args` with added 
          flags for actions

    """

    if args.forceDownload :
        args.download = True
    # Initiliaze action flags
    args.actionFlags = dict()
    # --query and --listId
    if (args.query is not None) and (args.listId is not None) :
        stdout.write("--query and --listId options cannot be specified "
                     "simultaneously\n"
                     "Use --help for details on usage.")
        sys.exit()
    # --query and no --listId
    elif (args.query is not None) and (args.listId is None) :
        _checkRetmax(args.retmax, stderr)
        _checkEmailOption(args, stderr)
        args.actionFlags["DoGenbankSearch"] = True
    # no --query and --listId
    elif (args.query is None) and (args.listId is not None) :
        _checkEmailOption(args, stderr)
        args.actionFlags["DoGetList"] = True
    # no --query and no --listId
    else :
        assert (args.query is None) and (args.listId is None)
        stderr.write("Please specify either --listId or --query\n"
                     "Use --help for details on usage.\n")
        sys.exit()
    return args

### *** _checkEmailOption(args, stderr = None)

def _checkEmailOption(args, stderr = None) :
    """Check that an email option was provided and setup Entrez email, produce 
    a message and exit if not.

    Args:
        args (namespace): Output from `parser.parse_args()`
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    Returns:
        Nothing, but setup Entrez.email or exit the program with a message to 
          stderr if no email was provided

    """
    if stderr is None :
        stderr = sys.stderr
    if args.email is None :
        stderr.write("To make use of NCBI's E-utilities, NCBI requires you to specify\n"
              "your email address with each request. In case of excessive\n"
              "usage of the E-utilities, NCBI will attempt to contact a user\n"
              "at the email address provided before blocking access to the\n"
              "E-utilities.")
        stderr.write("\nPlease provide your email address using --email")
        sys.exit()
    Entrez.email = args.email
    return

### *** _checkRetmax(retmax, stderr)

def _checkRetmax(retmax, stderr) :
    if retmax > 10000 :
        stderr.write("Error: --retmax cannot be greater than 10000.\n\n")
        stderr.write("Use --retmax 0 for unlimited number of returned results.\n\n"
                     "You can check how many entries your search term\n"
                     "returns on the GenBank website\n"
                     "(http://www.ncbi.nlm.nih.gov/genbank/)\n")
        sys.exit()

### ** pygenbank-extract-CDS

# Wishlist
#
# pygenbank-extract-CDS --outfmt orgn,gi,cds,nuc,pos,hash --hash sha512 --summary mySummaries *.gb
# pygenbank-extract-CDS *.gb > mySummaries # default outfmt and --summary
# pygenbank-extract-CDS --unique myUniqueCDS --hash md5sum *.gb > mySummaries # produce fasta
# pygenbank-extract-CDS --unique myUniqueCDS --hash md5sum --outfmt org,gi,cds,nuc,pos,hash *.gb > mySummaries
# pygenbank-extract-CDS --unique myUniqueCDS --hash md5sum --summary mySummaries --outfmt org,gi,cds,nuc,pos,hash *.gb

### *** _main_extract_CDS(args = None, stdout = None, stderr = None,
#                         gb_record_fmtdict, gb_cds_fmtdict)

def _main_extract_CDS(args = None, stdout = None, stderr = None,
                      gb_record_fmtdict = None,
                      gb_cds_fmtdict = None) :
    """Main function, used by the command line script "-extract-CDS". This function
    sends the arguments to :func:`_processArgsToLogic_extract_CDS` to determine
    which actions must be performed, and then performs the actions.

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)
        gb_record_fmtdict (dict): Dictionary mapping outfmt specifiers to 
          functions to extract the corresponding information from GenBank 
          records. If None, default is _GB_RECORD_FMTDICT.
        gb_cds_fmtdict (dict): Dictionary mapping outfmt specifiers to 
          functions to extract the corresponding information from GenBank
          CDS. If None, default is _GB_CDS_FMTDICT.

    Returns:
        None

    """
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    if gb_record_fmtdict is None :
        gb_record_fmtdict = _GB_RECORD_FMTDICT
    if gb_cds_fmtdict is None :
        gb_cds_fmtdict = _GB_CDS_FMTDICT
    # Process arguments
    if args is None :
        parser = _makeParser_extract_CDS()
        args = parser.parse_args()
    args = _processArgsToLogic_extract_CDS(args, stdout, stderr,
                                           gb_record_fmtdict, gb_cds_fmtdict)
    # Go through the input files
    uniqueSeq = dict()
    i_file = 0
    for fi in args.genbank_records :
        i_file += 1
        if args.verbose :
            stderr.write(time.asctime() + " - " +
                         "Processing file " + str(i_file) + " : " +
                         os.path.basename(fi) + " - " +
                         "N unique seq : " + str(len(uniqueSeq.keys())) + "\n")
        record = SeqIO.parse(fi, "genbank")
        for r in record :
            if not args.actionFlags.get("DoCount", False) :
                (summaryString, uniqueSeq, newSeq) = (
                    _summarizeRecord(r, args.outfmt, args.hash, uniqueSeq))
                stdout.write(summaryString)
            else :
                count = len([x for x in r.features if x.type == "CDS"])
                stdout.write(r.annotations["gi"] + "\t" + str(count) + "\n")
    # Write unique sequences
    if args.actionFlags.get("DoUniqueSequences", False) :
        with open(args.unique, "w") as fo :
            for (k, v) in uniqueSeq.items() :
                fo.write(">" + k + "\n")
                fo.write(v + "\n")
    
### *** _makeParser_extract_CDS()

def _makeParser_extract_CDS() :
    """Build the argument parser for the main script "-extract-CDS".
    
    returns:
        argparse.ArgumentParser() object: An argument parser object ready to be
        used to parse the command line arguments
    """
    parser = argparse.ArgumentParser(
        description = SCRIPT_DESCRIPTION_EXTRACT_CDS)
    # input files
    parser.add_argument("genbank_records", type = str, nargs = "+",
                        help = "Filename(s) for GenBank record(s)")
    # --hash
    parser.add_argument("--hash", metavar = "HASH_ALGORITHM",
                        choices = ["md5", "sha1", "sha224", "sha256", "sha384",
                                   "sha512"],
                        default = "md5",
                        help = "Hash algorithm to use for unique sequence signature "
                               "(default md5)")
    # --unique    
    parser.add_argument("-u", "--unique", metavar = "UNIQUE_FASTA",
                        type = str, default = None,
                        help = "Filename for unique fasta amino acid sequences")
    # --outfmt
    parser.add_argument("-o", "--outfmt", metavar = "FORMAT",
                        type = str, default = "dscrp,gi,loc,hash,prot",
                        help = "Comma-separated list of output column "
                        "identifiers for summaries. Allowed identifiers are: "
                        "id, name, dscrp, acc, date, gi, orgn, src (record "
                        "attributes) and loc, prot, nuc, prod, gene, hash "
                        "(CDS attributes). Default is dscrp,gi,loc,hash,prot.")
    # --count
    parser.add_argument("-c", "--count", action = "store_true",
                        help = "Count the number of CDS per record instead of "
                        " producing CDS summaries. Cannot be used at the same "
                        "time as --unique.")
    # --verbose
    parser.add_argument("-v", "--verbose", action = "store_true",
                        help = "Send verbose messages to stderr during execution")
    return parser

### *** _processArgsToLogic_extract_CDS(args, stdout, stderr)

def _processArgsToLogic_extract_CDS(args, stdout, stderr,
                                    gb_record_fmtdict,
                                    gb_cds_fmtdict) :
    """Process the command line arguments and determine the action logic for the
    :func:`_main_extract_CDS` function.

    Args:
        args (namespace): Argument namespace
        stdout (file): stdout stream
        stderr (file): stderr stream
        gb_record_fmtdict (dict): Dictionary mapping outfmt specifiers to 
          functions to extract the corresponding information from GenBank 
          records
        gb_cds_fmtdict (dict): Dictionary mapping outfmt specifiers to 
          functions to extract the corresponding information from GenBank
          CDS

    Returns:
        namespace: The argument namespace given input in `args` with added 
          flags for actions

    """
    # Initialize action flags
    args.actionFlags = dict()
    # Check for --count
    if args.count :
        if (args.unique is not None) :
            stderr.write("--count and --unique cannot be used simultanesouly.\n"
                         "Use --help for details on usage.\n")
            sys.exit()
        else :
            args.actionFlags["DoCount"] = True
            return args
    # Produce unique sequences
    if args.unique is not None :
        args.actionFlags["DoUniqueSequences"] = True
    # Process outfmt
    args.outfmt = _processOutfmtArg(args.outfmt, stderr, gb_record_fmtdict,
                                    gb_cds_fmtdict)
    # Process hash
    hashFunctions = { "md5" : hashlib.md5,
                      "sha1" : hashlib.sha1,
                      "sha224" : hashlib.sha224,
                      "sha256" : hashlib.sha256,
                      "sha384" : hashlib.sha384,
                      "sha512" : hashlib.sha512}
    args.hash = hashFunctions[args.hash]
    # Return
    return args

### *** _processOutfmtArg(outfmt, gb_record_fmtdict, gb_cds_fmtdict)

def _processOutfmtArg(outfmt, stderr, gb_record_fmtdict, gb_cds_fmtdict) :
    """Check that all outfmt specifiers are allowed and return the splitted outfmt
    specifiers.

    Args:
        outfmt (str): String from the command line argument ``--outfmt``
        stderr (file): stderr stream
        gb_record_fmtdict (dict): Dictionary mapping outfmt specifiers to 
          functions to extract the corresponding information from GenBank 
          records
        gb_cds_fmtdict (dict): Dictionary mapping outfmt specifiers to 
          functions to extract the corresponding information from GenBank
          CDS

    Returns:
        list: List of outfmt specifiers ready to pass to 
          :func:`_CDSinfo`

    """
    outfmt_keys = outfmt.split(",")
    records_keys = set(gb_record_fmtdict.keys())
    cds_keys = set(gb_cds_fmtdict.keys())
    assert records_keys & cds_keys == set()
    if not all([x in records_keys | cds_keys for x in outfmt_keys]) :
        wrong_keys = [x for x in outfmt_keys if x not in records_keys | cds_keys]
        stderr.write("Bad outfmt specifier. You provided:\n")
        stderr.write(str(sorted(outfmt_keys)) + "\n")
        stderr.write("Wrong specifier(s):\n")
        stderr.write(str(sorted(wrong_keys)) + "\n")
        stderr.write("Allowed values are:\n")
        stderr.write(str(sorted(list(records_keys | cds_keys))) + "\n")
    return outfmt_keys

### ** __main__

if (__name__ == "__main__") :
    funcSuffix = sys.argv[1].replace("-", "_")
    print(sys.argv[2:])
    args = eval("_makeParser_" + funcSuffix + "().parse_args(sys.argv[2:])")
    print(args)
    exec "_main_" + funcSuffix + "(args = args)"
