Command-line scripts
====================

pygenbank-search
----------------

``pygenbank-search`` is a tool to perform searches on GenBank and to retrieve
GenBank records, either as document summaries or as full records.

The user has to provide an email address for use of the Entrez resource.

Type ``pygenbank-search --help`` for detailed usage.

See http://www.ncbi.nlm.nih.gov/books/NBK49540/ for more details on GenBank
search queries.

Examples
********

Please use your own email address as the ``--email`` argument.

* Search GenBank and retrieve document summaries::

    pygenbank-search --query "hemoglobin AND mammal" --retmax 10000 --email "name@address" > mySearch
    less -S mysearch
  
* Search GenBank and retrieve full records::

    mkdir gbResults # Records will be saved here
    pygenbank-search -q "myoglobin AND sperm whale" -e "name@address" -d -o gbResults

* Specify a length range in the GenBank query::

    pygenbank-search -q "carcinus maenas" -e "name@address" -r 10000 > mySearch
    pygenbank-search -q "carcinus maenas AND 1000:100000[SLEN]" -e "name@address" -r 10000 > mySearchLength

* Specify a taxon in the GenBank query::

    MY_QUERY="complete genome AND staphylococcus aureus [PORGN]"
    pygenbank-search -q "$MY_QUERY" -e "name@address" -r 10000 > mySearch

* More complex query to get all complete *Staphylococcus aureus* genomes (up to
  10000)::

    MY_QUERY="complete genome AND staphylococcus aureus [PORGN] AND 1000000:10000000 [SLEN]"
    echo $MY_QUERY
    pygenbank-search -q "$MY_QUERY" -e "name@address" -r 10000 > mySearch

pygenbank-extract-CDS
---------------------

``pygenbank-extract-CDS`` is a tool to extract CDS summaries from GenBank
records and to produce fasta file with unique amino-acid sequences if needed
(e.g. to prepare a clustering analysis).

Type ``pygenbank-extract-CDS --help`` for detailed usage.

Examples
********

* Get CDS summaries for all GenBank files in the current directory::

    pygenbank-extract-CDS *.gb > mySummaries

How to profile command-line script execution
--------------------------------------------

`cprofilev` is a convenient tool to visualize the results of a profiling run
of a Python script::

  sudo pip install cprofilev

`cprofilev` can be used to profile the execution of a script this way::

  python -m cprofilev myScript.py [args]

The ouput is visible at the address ``http://localhost:4000``.

Using `cprofilev` with the command-line scripts
***********************************************

`pygenbank-search` and `pygenbank-extract-CDS` use entry points in the
`genbank.py` module, and cannot be called directly with the Python interpreter
to use the `cprofilev` module at the same time (at least I didn't find a way to
do it for now).

To solve this problem, there is a bit of code added at the end of the
`genbank.py` module to make it callable from the python interpreter. The module
can then be called with::

  python genbank.py search [args]
  python genbank.py extract-CDS [args]

where `[args]` are passed to the corresponding `_main_...` functions. For
example::

  python genbank.py search -q "hemocyanin" > summaries
  python genbank.py extract-CDS --help

Note that the full path to `genbank.py` must be provided (so here we assume we
are running the profiling run from within the module folder).
  
To perform a profiling run::

  python -m cprofilev genbank.py extract-CDS -u toto.fasta *.gb > summaries

and then visit ``http://localhost:4000`` with a web browser.
