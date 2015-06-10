|docs| |build status| |coverage|

Description
===========

Python tool to provide a simple interface to NCBI's GenBank database.

The code is written for Python 2.7.

PyGenBank is built on top of Biopython and can be used to perform a search on
GenBank and download records.

Installation
============

**To install the module and the command line tools**, type::

  sudo pip install git+https://github.com/matthieu-bruneaux/pygenbank

You can test if the installation worked with::

  pygenbank-search -h
  pygenbank-extract-CDS -h

and from Python::

  import genbank as gb
  dir(gb)
  
**To remove the module and the command line tools**, type::

  sudo pip uninstall pygenbank 
   
Documentation
=============

Documentation for the project is available on `Read the Docs <http://pygenbank.readthedocs.org/en/latest/>`_.

.. |docs| image:: https://readthedocs.org/projects/pygenbank/badge/?version=latest
   :target: http://pygenbank.readthedocs.org/en/latest/
   :alt: 'Docs'
.. |build status| image:: https://travis-ci.org/matthieu-bruneaux/pygenbank.svg?branch=master
   :target: https://travis-ci.org/matthieu-bruneaux/pygenbank?branch%3Dmaster
   :alt: 'Build status'
.. |coverage| image:: https://coveralls.io/repos/matthieu-bruneaux/pygenbank/badge.svg?branch=master
   :target: https://coveralls.io/r/matthieu-bruneaux/pygenbank?branch%3Dmaster
   :alt: 'Coverage'
