Installation
============

Easy install with pip
---------------------

Local install (no sudo rights, recommended)
+++++++++++++++++++++++++++++++++++++++++++

To install the module and the command line tools, without sudo rights, type
(thanks to this `post from Kaz'hack
<http://kazhack.org/?post/2014/12/12/pip-gem-install-without-sudo>`_ for
installing without sudo rights)::

  pip install --user --upgrade git+https://github.com/matthieu-bruneaux/pygenbank

This will install the executables in ``~/.local/bin/`` and the Python module in
``~/.local/lib/python2.7/site-packages/``.

To run the executables from the command line, you might need to add
``~/.local/bin/`` to the ``$PATH`` variable in your ``~/.bashrc`` file::

  # Lines to add in your ~/.bashrc file, if needed
  PATH=$PATH:~/.local/bin
  export PATH

System-wide install (with sudo rights, use with care!)
++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you want to install the module and the command line tools system-wide
(i.e. accessible for anyone), you can run ``pip`` with sudo rights (this is not
the recommended way to install!)::

  sudo pip install --upgrade git+https://github.com/matthieu-bruneaux/pygenbank

Test installation
+++++++++++++++++

You can test if the installation worked with::

  pygenbank-search -h
  pygenbank-extract-CDS -h

and from Python::

  import genbank as gb
  dir(gb)

Uninstall
+++++++++
  
To remove the module and the command line tools, type::

  pip uninstall pygenbank -y

if pygenbank was installed without sudo rights. If it was installed with sudo
rights, then you will need to use::

  sudo pip uninstall pygenbank-y
   
Installation from a cloned repository
-------------------------------------

First, clone the GitHub repository::

  git clone https://github.com/matthieu-bruneaux/pygenbank
  cd pygenbank
  
A Makefile is provided with the Python project folder. Just type `make` to have
a short summary of the different targets available.

To install the module::

  make install

This will install the module for the current user (using ``pip
install --user -e .``). See above how to add ``~/.local/bin`` to the ``PATH``
variable if the command line tools are not accessible after this install.
  
To remove the module::
     
  make uninstall

You can also run some tests and regenerate the documentation::
    
  make tests
  make doc
