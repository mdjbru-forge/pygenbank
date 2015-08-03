### * Variables

PYTHON_MODULE=genbank
PYTHON_MODULE_EGG=$(PYTHON_MODULE).egg-info

COVERED_PACKAGES=$(PYTHON_MODULE)
SPHINX_DOC_FOLDER=docs/

# ### ** Main script
# PYTHON=python
# MODULE_NAME=pydep
# MYSCRIPT_NOPY=pydep
# MYSCRIPT=$(MODULE_NAME)/$(MYSCRIPT_NOPY).py

# ### ** Tests
# TEST_DIR=tests
# TEST_SCRIPT=tests.py

# ### ** Examples
# EXAMPLE_MOD_NOPY=exampleModule

### * Help

Help:
	@echo "Makefile for the $(PYTHON_MODULE) Python module                   "
	@echo "                                                                  "
	@echo "Type: \"make <target>\" where <target> is one of the following:   "
	@echo "                                                                  "
	@echo "  test             Run the tests with coverage output             "
	@echo "  doc              Run Sphinx to make the docs                    "
	@echo "  clean            Remove generated doc, tests and pyc files      "
	@echo "                                                                  "
	@echo "  install          Install the module and command-line tools      "
	@echo "  uninstall        Uninstall the module                           "

### * Main targets

### ** test
tests: test
test:
	nosetests tests/ --with-coverage --cover-package=$(COVERED_PACKAGES) --cover-html \
	  --with-html --html-file=tests/nosetests.html
	@echo -e "\nThe coverage results are accessible from cover/index.html"
	@echo "The html version of the test results are accessible from tests/nosetests.html"

### ** doc
docs: doc
doc:
	sphinx-apidoc -o $(SPHINX_DOC_FOLDER) ./ setup.py
	cd $(SPHINX_DOC_FOLDER); make html
	@echo -e "\nThe documentation is accessible from $(SPHINX_DOC_FOLDER)_build/html/index.html"

### ** clean
clean:
	rm -f .coverage
	rm -fr cover
	rm -f tests/nosetests.html
	rm -fr docs/_build
	find . -name \*.pyc -type f -delete

### ** install
install:
	rm -fr $(PYTHON_MODULE_EGG)
	pip install --user .

### ** uninstall
uninstall:
	pip uninstall -y $(PYTHON_MODULE)
