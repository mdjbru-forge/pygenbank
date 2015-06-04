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

   import genbank

   # Setup your email address for Entrez
   genbank.Entrez.email = "yourname@youraddress"


Search GenBank and retrieve record id
*************************************

Performing a GenBank search is as simple as::

   # Perform a GenBank search
   mySearch = genbank.genbankSearch(term = "hemocyanin", retmax = 100)

The search results can be used to get summaries of the results and apply some
simple filtering on the record id before proceeding to the actual record
downloading::
   
   # Get the summaries from the results
   summaries = genbank.genbankGetDocSum(mySearch)

   # Extract the id of interest ("Gi" field)
   myId = [x["Gi"] for x in summaries if int(x["Length"]) < 10000]
   len(myId)
   
Download the GenBank records
****************************

::

   # Download the GenBank records
   genbank.genbankDownloadRecords(idList = myId, destDir = ".", batchSize = 20)

Functions
---------
   
.. graphviz::

  digraph G {
  rankdir=LR;
  subgraph cluster1{
  node[shape=box,style=filled,fillcolor="#dfaf8f"];
  _processOutfmtArg;
  _fileLinesToList;
  _checkEmailOption;
  _processArgsToLogic_search;
  _genbankParseDocSumXML;
  _processArgsToLogic_extract_CDS;
  _genbankSummarizeRecord;
  _genbankDownloadBatch;
  _genbankGetDocSumXML;
  _makeParser_search;
  _genbankMakeSummaryForCDS;
  _checkRetmax;
  _genbankGetRecordBatch;
  _makeParser_extract_CDS;
  _genbankGetProteinHashFromCDS;
  node[shape=box,style=filled,fillcolor="#7cb8bb"];
  genbankDownloadRecords;
  genbankGetDocSumFromId;
  genbankWriteDocSums;
  genbankSearch;
  genbankGetDocSum;
  node[shape=box,style=filled,fillcolor="#9fc59f"];
  _main_search;
  _main_extract_CDS;
  genbankDownloadRecords -> _genbankDownloadBatch;
  genbankGetDocSumFromId -> _genbankParseDocSumXML;
  genbankGetDocSumFromId -> _genbankGetDocSumXML;
  _processArgsToLogic_search -> _checkRetmax;
  _processArgsToLogic_search -> _checkEmailOption;
  _genbankSummarizeRecord -> _genbankGetProteinHashFromCDS;
  _genbankSummarizeRecord -> _genbankMakeSummaryForCDS;
  _genbankDownloadBatch -> _genbankGetRecordBatch;
  _processArgsToLogic_extract_CDS -> _processOutfmtArg;
  _main_search -> _makeParser_search;
  _main_search -> genbankGetDocSumFromId;
  _main_search -> genbankSearch;
  _main_search -> _fileLinesToList;
  _main_search -> genbankGetDocSum;
  _main_search -> genbankDownloadRecords;
  _main_search -> genbankWriteDocSums;
  _main_search -> _processArgsToLogic_search;
  genbankGetDocSum -> _genbankParseDocSumXML;
  genbankGetDocSum -> _genbankGetDocSumXML;
  _main_extract_CDS -> _processArgsToLogic_extract_CDS;
  _main_extract_CDS -> _makeParser_extract_CDS;
  _main_extract_CDS -> _genbankSummarizeRecord;
  }
  }

.. automodule:: genbank
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
