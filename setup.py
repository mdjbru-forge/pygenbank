# Written using resources from:
# https://packaging.python.org/en/latest/distributing.html#working-in-development-mode
# https://github.com/pypa/sampleproject/blob/master/setup.py
# https://docs.python.org/2/distutils/examples.html

from setuptools import setup

setup(name = "pygenbank",
      version = "0.0.1",
      py_modules = ["genbank"],
      entry_points = {
          "console_scripts" : [
              "pygenbank-search=genbank:_main_search",
              "pygenbank-extract-CDS=genbank:_main_extract_CDS"
          ]
      }
)
