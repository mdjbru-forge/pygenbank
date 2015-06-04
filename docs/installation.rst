Installation
============

Easy install with pip
---------------------

**To install the module and the command line tools**, type::

  sudo pip install git+https://github.com/matthieu-bruneaux/pygenbank

You can test if the installation worked with::

  pygenbank-search -h
  pygenbank-extract-CDS -h

and from Python::

  import genbank as gb
  dir(gb)
  
**To remove the module and the command line tools**, type::

  sudo pip uninstall pyGenBank 
   
Installation from a cloned repository
-------------------------------------

First, clone the GitHub repository::

  git clone https://github.com/matthieu-bruneaux/pygenbank
  cd pygenbank
  
A Makefile is provided with the Python project folder. Just type `make` to have
a short summary of the different targets available.

To install the module::

  sudo make install

To remove the module::
     
  sudo make uninstall

You can also run some tests and regenerate the documentation::
    
  make tests
  make doc
