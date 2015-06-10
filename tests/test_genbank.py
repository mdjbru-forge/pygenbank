### * Description

# Test script for genbank module

### ** Requirements

# sudo pip install nose
# sudo pip install coverage

### ** Usage

# nosetests ./
# nosetests ./ --with-coverage --cover-html --cover-package=genbank

### * Setup

### ** Import

import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(".."))
import StringIO

import genbank as mod

### ** Parameters

# http://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
CUR_DIR = os.path.dirname(os.path.realpath(__file__))

### * Run

### ** Test _makeParser_search

class TestMakeParser_Search(unittest.TestCase) :

### *** setUp and tearDown

    def setUp(self) :
        # We redirect stderr so that there is no error message displayed when
        # we test for parser error message.
        self.stderr = sys.stderr
        sys.stderr = StringIO.StringIO()
        # Parser
        self.parser = mod._makeParser_search()

    def tearDown(self) :
        sys.stderr.close()
        sys.stderr = self.stderr

### *** Test --email

    def test_email_none(self) :
        args = []
        result = self.parser.parse_args(args)
        self.assertIsNone(result.email)

    def test_email_long(self) :
        args = ["--email", "toto@toto.com"]
        result = self.parser.parse_args(args)
        expected = "toto@toto.com"
        self.assertEqual(result.email, expected)

    def test_email_short(self) :
        args = ["-e", "toto@toto.com"]
        result = self.parser.parse_args(args)
        expected = "toto@toto.com"
        self.assertEqual(result.email, expected)

### *** Test --listId

    def test_listId_none(self) :
        args = []
        result = self.parser.parse_args(args)
        self.assertIsNone(result.listId)

    def test_listId_long(self) :
        args = ["--listId", "list.txt"]
        result = self.parser.parse_args(args)
        expected = "list.txt"
        self.assertEqual(result.listId, expected)
        
    def test_listId_short(self) :
        args = ["-l", "list.txt"]
        result = self.parser.parse_args(args)
        expected = "list.txt"
        self.assertEqual(result.listId, expected)

### *** Test --query

    def test_query_none(self) :
        args = []
        result = self.parser.parse_args(args)
        self.assertIsNone(result.query)

    def test_query_long(self) :
        args = ["--query", "hemocyanin AND 100:1000 [SLEN]"]
        result = self.parser.parse_args(args)
        expected = "hemocyanin AND 100:1000 [SLEN]"
        self.assertEqual(result.query, expected)
        
    def test_query_short(self) :
        args = ["-q", "hemocyanin AND 100:1000 [SLEN]"]
        result = self.parser.parse_args(args)
        expected = "hemocyanin AND 100:1000 [SLEN]"
        self.assertEqual(result.query, expected)

### *** Test --retmax

    def test_retmax_none(self) :
        args = []
        result = self.parser.parse_args(args)
        expected = 0
        self.assertEqual(result.retmax, expected)

    def test_retmax_long(self) :
        args = ["--retmax", "10"]
        result = self.parser.parse_args(args)
        expected = 10
        self.assertEqual(result.retmax, expected)
        
    def test_retmax_short(self) :
        args = ["-r", "10"]
        result = self.parser.parse_args(args)
        expected = 10
        self.assertEqual(result.retmax, expected)

    def test_retmax_string(self) :
        args = ["-r", "toto"]
        with self.assertRaises(SystemExit) :
            result = self.parser.parse_args(args)

### *** Test --download

    def test_download_none(self) :
        args = []
        result = self.parser.parse_args(args)
        expected = False
        self.assertEqual(result.download, expected)

    def test_download_long(self) :
        args = ["--download"]
        result = self.parser.parse_args(args)
        expected = True
        self.assertEqual(result.download, expected)
        
    def test_download_short(self) :
        args = ["-d"]
        result = self.parser.parse_args(args)
        expected = True
        self.assertEqual(result.download, expected)

    def test_download_error(self) :
        args = ["-d", "tata"]
        with self.assertRaises(SystemExit) :
            result = self.parser.parse_args(args)
        
### *** Test --forceDownload

    def test_forceDownload_none(self) :
        args = []
        result = self.parser.parse_args(args)
        expected = False
        self.assertEqual(result.forceDownload, expected)

    def test_forceDownload_long(self) :
        args = ["--forceDownload"]
        result = self.parser.parse_args(args)
        expected = True
        self.assertEqual(result.forceDownload, expected)
        
    def test_forceDownload_short(self) :
        args = ["-f"]
        result = self.parser.parse_args(args)
        expected = True
        self.assertEqual(result.forceDownload, expected)

    def test_forceDownload_error(self) :
        args = ["-f", "tata"]
        with self.assertRaises(SystemExit) :
            result = self.parser.parse_args(args)

### *** Test --outputDir

    def test_outputDir_none(self) :
        args = []
        result = self.parser.parse_args(args)
        expected = "."
        self.assertEqual(result.outputDir, expected)

    def test_outputDir_long(self) :
        args = ["--outputDir", "records"]
        result = self.parser.parse_args(args)
        expected = "records"
        self.assertEqual(result.outputDir, expected)
        
    def test_outputDir_short(self) :
        args = ["-o", "records"]
        result = self.parser.parse_args(args)
        expected = "records"
        self.assertEqual(result.outputDir, expected)

### *** Test --batchSize

    def test_batchSize_none(self) :
        args = []
        result = self.parser.parse_args(args)
        expected = 5
        self.assertEqual(result.batchSize, expected)

    def test_batchSize_long(self) :
        args = ["--batchSize", "10"]
        result = self.parser.parse_args(args)
        expected = 10
        self.assertEqual(result.batchSize, expected)
        
    def test_batchSize_short(self) :
        args = ["-b", "10"]
        result = self.parser.parse_args(args)
        expected = 10
        self.assertEqual(result.batchSize, expected)

    def test_batchSize_string(self) :
        args = ["-b", "toto"]
        with self.assertRaises(SystemExit) :
            result = self.parser.parse_args(args)
            
### *** Test --delay

    def test_delay_none(self) :
        args = []
        result = self.parser.parse_args(args)
        expected = 15
        self.assertEqual(result.delay, expected)

    def test_delay_000(self) :
        args = ["--delay", "10"]
        result = self.parser.parse_args(args)
        expected = 10
        self.assertEqual(result.delay, expected)
        
    def test_delay_string(self) :
        args = ["--delay", "toto"]
        with self.assertRaises(SystemExit) :
            result = self.parser.parse_args(args)

### ** Test _checkEmailOption

class TestCheckEmailOption(unittest.TestCase) :

### *** setUp and tearDown

    def setUp(self) :
        self.parser = mod._makeParser_search()
        self.stderr = StringIO.StringIO()
        self.oldEntrezEmail = mod.Entrez.email
        self.oldSysStderr = sys.stderr
        sys.stderr = StringIO.StringIO()
        
    def tearDown(self) :
        mod.Entrez.email = self.oldEntrezEmail
        sys.stderr = self.oldSysStderr
        
### *** Test

    def test_checkEmailOption_None(self) :
        args = self.parser.parse_args([])
        with self.assertRaises(SystemExit) :
            result = mod._checkEmailOption(args, self.stderr)

    def test_checkEmailOptions_000(self) :
        args = self.parser.parse_args(["-e", "darwin@evolution.fi"])
        mod._checkEmailOption(args, self.stderr)
        self.assertEqual(mod.Entrez.email, "darwin@evolution.fi")

    def test_checkEmailOption_initial_Entrez_email(self) :
        self.assertIsNone(mod.Entrez.email)

    def test_checkEmailOption_stderr_None_000(self) :
        args = self.parser.parse_args([])
        with self.assertRaises(SystemExit) :
            result = mod._checkEmailOption(args)

    def test_checkEmailOption_stderr_None_001(self) :
        args = self.parser.parse_args([])
        try :
            result = mod._checkEmailOption(args)
        except SystemExit :
            sys.stderr.seek(0)
            expected = ("To make use of NCBI's E-utilities, NCBI requires you to "
            "specify\nyour email address with each request. In case of excessive\n"
            "usage of the E-utilities, NCBI will attempt to contact a user\n"
            "at the email address provided before blocking access to the\n"
            "E-utilities.\nPlease provide your email address using --email")
            self.assertEqual(expected, sys.stderr.read())

### ** Test search

class TestSearch(unittest.TestCase) :

# Note: This test is a bit artificial - there is not much we can test without
# actually connecting to the GenBank server. Here we will just check that the
# correct function arguments are passed to Entrez.esearch.

### *** setUp and tearDown

    def setUp(self) :
        self.oldEntrezEsearch = mod.Entrez.esearch
        self.oldEntrezRead = mod.Entrez.read
        # New Entrez.esearch
        def esearch(db, term, retmax, usehistory) :
            o = StringIO.StringIO()
            o.write("\t".join([db, term, str(retmax), usehistory]))
            return o
        mod.Entrez.esearch = esearch
        # New Entrez.read
        def read(handle) :
            return handle.getvalue()
        mod.Entrez.read = read            

    def tearDown(self) :
        mod.Entrez.esearch = self.oldEntrezEsearch
        mod.Entrez.read = self.oldEntrezRead
        
### *** Test

    def test_search_db(self) :
        result = mod.search("toto", 20)
        self.assertEqual(result.split("\t")[0], "nuccore")
       
    def test_search_term(self) :
        result = mod.search("toto", 20)
        self.assertEqual(result.split("\t")[1], "toto")

    def test_search_retmax(self) :
        result = mod.search("toto", 20)
        self.assertEqual(result.split("\t")[2], "20")

    def test_search_usehistory(self) :
        result = mod.search("toto", 20)
        self.assertEqual(result.split("\t")[3], "y")

### ** Test _getRecordBatch

class TestGenbankGetRecordBatch(unittest.TestCase) :

### *** setUp and tearDown

    def setUp(self) :
        self.oldEntrezEfetch = mod.Entrez.efetch
        # New Entrez.efetch
        def efetch(db, rettype, retmode, id) :
            o = StringIO.StringIO()
            o.write("\t".join([db, rettype, retmode, id]))
            o.seek(0)
            return o
        mod.Entrez.efetch = efetch

    def tearDown(self) :
        mod.Entrez.efetch = self.oldEntrezEfetch

### *** Test

    def test_getRecordBatch_db(self) :
        inputList = ["tata", "tete", "titi", "toto"]
        result = mod._getRecordBatch(inputList)
        expected = "nuccore"
        self.assertEqual(result.split("\t")[0], expected)
                                           
    def test_getRecordBatch_rettype(self) :
        inputList = ["tata", "tete", "titi", "toto"]
        result = mod._getRecordBatch(inputList)
        expected = "gbwithparts"
        self.assertEqual(result.split("\t")[1], expected)

    def test_getRecordBatch_retmode(self) :
        inputList = ["tata", "tete", "titi", "toto"]
        result = mod._getRecordBatch(inputList)
        expected = "text"
        self.assertEqual(result.split("\t")[2], expected)

    def test_getRecordBatch_id(self) :
        inputList = ["tata", "tete", "titi", "toto"]
        result = mod._getRecordBatch(inputList)
        expected = "tata,tete,titi,toto"
        self.assertEqual(result.split("\t")[3], expected)

    def test_getRecordBatch_empty_list(self) :
        inputList = []
        result = mod._getRecordBatch(inputList)
        expected = ""
        self.assertEqual(result.split("\t")[3], expected)
                   
### ** Test _checkRetmax

class TestCheckRetmax(unittest.TestCase) :

### *** Test

    def test_checkRetmax_zero(self) :
        stderr = StringIO.StringIO()
        retmax = 0
        result = mod._checkRetmax(retmax, stderr)
        self.assertIsNone(result)
        
    def test_checkRetmax_000(self) :
        stderr = StringIO.StringIO()
        retmax = 1000
        result = mod._checkRetmax(retmax, stderr)
        self.assertIsNone(result)

    def test_checkRetmax_error(self) :
        stderr = StringIO.StringIO()
        retmax = 10001
        with self.assertRaises(SystemExit) :
            mod._checkRetmax(retmax, stderr)

### ** Test _fileLinesToList

class TestFileLinesToList(unittest.TestCase) :

### *** Test

    def test_fileLinesToList_empty(self) :
        inputFile = os.path.join(CUR_DIR, "test_genbank/test_fileLinesToList/input_empty")
        result = mod._fileLinesToList(inputFile)
        expected = []
        self.assertEqual(result, expected)

    def test_fileLinesToList_blank(self) :
        inputFile = os.path.join(CUR_DIR, "test_genbank/test_fileLinesToList/input_blank")
        result = mod._fileLinesToList(inputFile)
        expected = []
        self.assertEqual(result, expected)

    def test_fileLinesToList_000(self) :
        inputFile = os.path.join(CUR_DIR, "test_genbank/test_fileLinesToList/input1")
        result = mod._fileLinesToList(inputFile)
        expected = ["toto", "tata", "titi", "tutu"]
        self.assertEqual(result, expected)

    def test_fileLinesToList_001(self) :
        inputFile = os.path.join(CUR_DIR, "test_genbank/test_fileLinesToList/input2")
        result = mod._fileLinesToList(inputFile)
        expected = ["toto", "toto tata  tutu", "titi"]
        self.assertEqual(result, expected)

    def test_fileLinesToList_error(self) :
        inputFile = "nonExistent"
        with self.assertRaises(IOError) :
            result = mod._fileLinesToList(inputFile)

### ** Test _processArgsToLogic_search

class TestProcessArgsToLogic_search(unittest.TestCase) :

### *** setUp

    def setUp(self) :
        self.parser = mod._makeParser_search()
        self.stdout = StringIO.StringIO()
        self.stderr = StringIO.StringIO()
        
### *** Test

    def test_no_query_no_listId(self) :
        args = self.parser.parse_args([])
        with self.assertRaises(SystemExit) :
            result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
            
    def test_query_listId(self) :
        args = self.parser.parse_args(["--query", "toto", "--listId", "tata"])
        with self.assertRaises(SystemExit) :
            result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)

    def test_no_email_query(self) :
        args = self.parser.parse_args(["--query", "toto"])
        with self.assertRaises(SystemExit) :
            result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)

    def test_no_email_listId(self) :
        args = self.parser.parse_args(["--listId", "toto"])
        with self.assertRaises(SystemExit) :
            result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)

    def test_query_000(self) :
        args = self.parser.parse_args(["--query", "toto", "--email", "myEmail"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.actionFlags.get("DoGenbankSearch", False))

    def test_query_001(self) :
        args = self.parser.parse_args(["--query", "toto", "--email", "myEmail"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertFalse(result.actionFlags.get("DoGetList", False))
        
    def test_listId_000(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email", "myEmail"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.actionFlags.get("DoGetList", False))

    def test_listId_001(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email", "myEmail"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertFalse(result.actionFlags.get("DoGenbankSearch", False))

    def test_download_000(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email", "myEmail"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertFalse(result.download)

    def test_download_001(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email",
                                       "myEmail", "--download"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.download)

    def test_forceDownload_000(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email",
                                       "myEmail", "-f"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.download)

    def test_forceDownload_000b(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email",
                                       "myEmail", "-f"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.forceDownload)

    def test_forceDownload_001(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email",
                                       "myEmail", "-d", "-f"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.download)

    def test_forceDownload_001b(self) :
        args = self.parser.parse_args(["--listId", "toto", "--email",
                                       "myEmail", "-d", "-f"])
        result = mod._processArgsToLogic_search(args, self.stdout, self.stderr)
        self.assertTrue(result.forceDownload)

### ** Test _parseDocSumXML

class TestGenbankParseDocSumXML(unittest.TestCase) :

### *** setUp and tearDown

    def setUp(self) :
        with open(os.path.join(CUR_DIR,
                               "test_genbank",
                               "test_parseDocSumXML",
                               "docSum_noResults.xml"), "r") as fo :
            self.xml_noResults = fo.read()
        with open(os.path.join(CUR_DIR,
                               "test_genbank",
                               "test_parseDocSumXML",
                               "docSum_results.xml"), "r") as fo :
            self.xml_results = fo.read()

### *** Test

    def test_parseDocSumXML_noResults(self) :
        result = mod._parseDocSumXML(self.xml_noResults)
        expected = []
        self.assertEqual(result, expected)
        
    def test_parseDocSumXML_results_000(self) :
        result = mod._parseDocSumXML(self.xml_results)
        expected = 4
        self.assertEqual(len(result), expected)

    def test_parseDocSumXML_results_001(self) :
        result = mod._parseDocSumXML(self.xml_results)
        expected = set(["Status", "Comment", "Extra", "CreateDate", "Title",
                        "TaxId", "ReplacedBy", "Caption", "Length", "Flags",
                        "UpdateDate", "Gi"])
        self.assertEqual(set(result[0].keys()), expected)

### ** Test writeDocSums

class TestGenbankWriteDocSums(unittest.TestCase) :

### *** setUp and tearDown

    def setUp(self) :
        def loadXMLtoDict(filename) :
            file_path = os.path.join(CUR_DIR, "test_genbank/test_parseDocSumXML",
                                     filename)
            with open(file_path, "r") as fi :
                file_content = fi.read()
            return mod._parseDocSumXML(file_content)
        self.docsumsNoResults = loadXMLtoDict("docSum_noResults.xml")
        self.docsumsResults = loadXMLtoDict("docSum_results.xml")
        self.handle = StringIO.StringIO()
        
### *** Test

    def test_writeDocSums_noResults_000(self) :
        result = mod.writeDocSums(self.docsumsNoResults, self.handle)
        self.assertIsNone(result)

    def test_writeDocSums_noResults_000b(self) :
        result = mod.writeDocSums(self.docsumsNoResults, self.handle)
        output = self.handle.getvalue()
        expected = ""
        self.assertEqual(output, expected)

    def test_writeDocSums_results_000(self) :
        result = mod.writeDocSums(self.docsumsResults, self.handle)
        self.assertIsNone(result)

    def test_writeDocSums_results_000b(self) :
        result = mod.writeDocSums(self.docsumsResults, self.handle)
        output = self.handle.getvalue()
        with open(os.path.join(CUR_DIR, "test_genbank",
                               "test_parseDocSumXML",
                               "docSum_results.table"), "r") as fi :
                  expected = fi.read()
        self.assertEqual(output, expected)

### ** Test _recordIsWGS

class TestRecordIsWGS(unittest.TestCase) :

### *** setUp and tearDown

    def setUp(self) :
        self.recordNoWGS = (""
        "Dummy record\n"
        "\n"
        "This is a dummy record\n"
        "No WGS at beginning of line\n")
        self.recordWGS =  (""
        "Dummy record\n"
        "\n"
        "This is a dummy record\n"
        "WGS at beginning of line\n")
        self.recordWGStwice =  (""
        "Dummy record\n"
        "\n"
        "This is a dummy record\n"
        "WGS at beginning of line\n"
        "WGS again on another line\n")

### *** Test

    def test_recordIsWGS_noWGS(self) :
        result = mod._recordIsWGS(self.recordNoWGS)
        self.assertFalse(result)

    def test_recordIsWGS_WGS(self) :
        result = mod._recordIsWGS(self.recordWGS)
        expected = "WGS at beginning of line"
        self.assertEqual(result, expected)

    def test_recordIsWGS_error(self) :
        with self.assertRaises(Exception) :
            result = mod._recordIsWGS(self.recordWGStwice)

### ** Test _makeWGSurl

class TestMakeWGSurl(unittest.TestCase) :

### *** Test

    def test_makeWGSurl_bad_start_000(self) :
        wgs = ""
        with self.assertRaises(Exception) :
            result = mod._makeWGSurl(wgs)

    def test_makeWGSurl_bad_start_001(self) :
        wgs = " WGS "
        with self.assertRaises(Exception) :
            result = mod._makeWGSurl(wgs)
    
    def test_makeWGSurl_000(self) :
        wgs = "WGS         LBDC01000001-LBDC01000059"
        expected = "http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=LBDC01.1.gbff.gz"
        result = mod._makeWGSurl(wgs)
        self.assertEqual(result, expected)

    def test_makeWGSurl_001(self) :
        wgs = "WGS         AVVN01000001-AVVN01000071"
        expected = "http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=AVVN01.1.gbff.gz"
        result = mod._makeWGSurl(wgs)
        self.assertEqual(result, expected)

### ** Test _downloadWGS

class TestDownloadWGS(unittest.TestCase) :

# Great help from http://www.saltycrane.com/blog/2012/11/using-pythons-gzip-and-stringio-compress-data-memory/
    
### *** setUp and tearDown

    def setUp(self) :
        self.oldUrlopen = mod.urllib2.urlopen
        self.calledUrl = None
        def f(url) :
            self.calledUrl = url
            s = url + " going through compression/decompression"
            o = StringIO.StringIO()
            g = mod.gzip.GzipFile(fileobj = o, mode = "w")
            g.write(s)
            g.close()
            o.seek(0)
            return o
        mod.urllib2.urlopen = f
        self.oldStdout = sys.stdout
        sys.stdout = StringIO.StringIO()
        
    def tearDown(self) :
        mod.urllib2.urlopen = self.oldUrlopen
        sys.stdout = self.oldStdout
    
### *** Test

    def test_downloadWGS_url(self) :
        url = "testUrl"
        result = mod._downloadWGS(url)
        self.assertEqual(self.calledUrl, "testUrl")
        

    def test_downloadWGS_000(self) :
        url = "testUrl"
        result = mod._downloadWGS(url)
        expected = "testUrl going through compression/decompression"
        self.assertEqual(result, expected)

    def test_downloadWGS_error_000(self) :
        url = "testUrl"
        def f(x) :
            return StringIO.StringIO("toto")
        mod.urllib2.urlopen = f
        result = mod._downloadWGS(url)
        self.assertIsNone(result)

    def test_downloadWGS_error_001(self) :
        url = "testUrl"
        def f(x) :
            return StringIO.StringIO("toto")
        mod.urllib2.urlopen = f
        result = mod._downloadWGS(url)
        expected = "Not a gzipped file\n"
        sys.stdout.seek(0)
        self.assertEqual(sys.stdout.read(), expected)
