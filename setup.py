#!/usr/bin/env python3

import sys

if sys.version_info[0] < 3:
    sys.exit('Sorry, Python < 3.x is not supported')

# Try using setuptools first, if it's installed
from setuptools import setup, find_packages

# Need to add all dependencies to setup as we go!
setup(name='bcmap',
      packages=find_packages(),
      version='0.0.1',
      description="map barcodes to full sequences from a high-throughput mapping experiment",
      long_description=open("README.md").read(),
      author='Michael J. Harms',
      author_email='harmsm@gmail.com',
      url='https://github.com/harmslab/bcmap',
      download_url='',
      zip_safe=False,
      install_requires=["numpy"],
      classifiers=['Programming Language :: Python'])
