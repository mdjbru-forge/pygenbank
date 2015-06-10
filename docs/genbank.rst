genbank module
==============

Description
-----------

``genbank`` provides some light wrappers around Biopython functions to download
records from the GenBank server. The module can be used from within Python, or
through the command line tool ``pygenbank``.

Tutorial
--------

This is a simple tutorial to learn how to use the :mod:`genbank` module.

Setup the environment
*********************

::

   import genbank as gb

   # Setup your email address for Entrez
   gb.Entrez.email = "yourname@youraddress"


Search GenBank and retrieve record ids
**************************************

Performing a GenBank search is as simple as::

   # Perform a GenBank search
   mySearch = gb.search(term = "hemocyanin")

The returned value is a ``Bio.Entrez.Parser.DictionaryElement`` which contains
information about the returned results::

  mySearch.keys() # Available information
  mySearch["Count"] # Number of entries found
  mySearch["QueryTranslation"] # How the query was understood by GenBank
  mySearch["IdList"] # List of GenBank identifiers returned
   
Any query string that you would be using on the GenBank web page can be used
as the term::

  mySearch = gb.search(term = "hemocyanin AND lito* [ORGN]")
  mySearch["Count"]
  mySearch = gb.search(term = "citrate synthase AND mus m* [ORGN]

To get more details about the returned entries, you can fetch the record
summaries using the previous search result::

  summaries = gb.getDocSum(mySearch)
  
The search results can be used to get summaries of the results and apply some
simple filtering on the record id before proceeding to the actual record
downloading::
   
   # Get the summaries from the results
   summaries = gb.getDocSum(mySearch)

   # Extract the id of interest ("Gi" field)
   myId = [x["Gi"] for x in summaries if int(x["Length"]) < 10000]
   len(myId)
   
Download the GenBank records
****************************

::

   # Download the GenBank records
   genbank.downloadRecords(idList = myId, destDir = ".", batchSize = 20)

Functions
---------
   
.. graphviz::

  digraph G {
  rankdir=LR;
  subgraph cluster_1 {
  node[shape=box,style=filled,fillcolor="#dfaf8f"];
  _processOutfmtArg;
  _fileLinesToList;
  _makeSummaryForCDS;
  _checkEmailOption;
  _getRecordBatch;
  _recordIsWGS;
  _downloadBatch;
  _processArgsToLogic_extract_CDS;
  _parseDocSumXML;
  _getProteinHashFromCDS;
  _getDocSumXML;
  _makeParser_search;
  _checkRetmax;
  _summarizeRecord;
  _downloadWGS;
  _makeWGSurl;
  _makeParser_extract_CDS;
  _processArgsToLogic_search;
  node[shape=box,style=filled,fillcolor="#7cb8bb"];
  getDocSumFromId;
  downloadRecords;
  getDocSum;
  search;
  downloadWGS;
  writeDocSums;
  node[shape=box,style=filled,fillcolor="#9fc59f"];
  _main_search;
  _main_extract_CDS;
  getDocSum -> _getDocSumXML;
  getDocSum -> _parseDocSumXML;
  getDocSumFromId -> _getDocSumXML;
  getDocSumFromId -> _parseDocSumXML;
  downloadRecords -> _downloadBatch;
  _processArgsToLogic_extract_CDS -> _processOutfmtArg;
  downloadWGS -> _downloadWGS;
  downloadWGS -> _makeWGSurl;
  downloadWGS -> _recordIsWGS;
  _downloadBatch -> downloadWGS;
  _downloadBatch -> _getRecordBatch;
  _downloadBatch -> _recordIsWGS;
  _main_search -> search;
  _main_search -> getDocSum;
  _main_search -> _makeParser_search;
  _main_search -> downloadRecords;
  _main_search -> getDocSumFromId;
  _main_search -> _fileLinesToList;
  _main_search -> _processArgsToLogic_search;
  _main_search -> writeDocSums;
  _processArgsToLogic_search -> _checkRetmax;
  _processArgsToLogic_search -> _checkEmailOption;
  _main_extract_CDS -> _processArgsToLogic_extract_CDS;
  _main_extract_CDS -> _makeParser_extract_CDS;
  _main_extract_CDS -> _summarizeRecord;
  _summarizeRecord -> _getProteinHashFromCDS;
  _summarizeRecord -> _makeSummaryForCDS;
  }
  }

.. automodule:: genbank
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
