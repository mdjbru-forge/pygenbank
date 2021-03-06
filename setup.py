# Written using resources from:
# https://packaging.python.org/en/latest/distributing.html#working-in-development-mode
# https://github.com/pypa/sampleproject/blob/master/setup.py
# https://docs.python.org/2/distutils/examples.html

from setuptools import setup

setup(name = "pygenbank",
      version = "0.0.3a",
      py_modules = ["genbank"],
      entry_points = {
          "console_scripts" : [
              "pygenbank=genbank:main"
          ]
      }
)
