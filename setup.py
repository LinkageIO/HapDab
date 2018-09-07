#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
import os
from setuptools.command.develop import develop
from setuptools.command.install import install
from subprocess import check_call

import io
import re

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
    name = 'hapdab',
    version = find_version('hapdab','__init__.py'),

    description = 'An library for managing and analyzing genotypes and haplotypes',
    url = 'http://linkage.io',
    author = 'Rob Schaefer',
    license = "Copyright Linkage Analytics 2016-2018. Available under the MIT License",
    author_email = 'rob@linkage.io',


    classifiers=[
	# How mature is this project? Common values are
	#   3 - Alpha
	#   4 - Beta
	#   5 - Production/Stable
	'Development Status :: 4 - Beta',

	# Indicate who your project is intended for
	'Intended Audience :: Developers',
	'Topic :: Software Development :: Build Tools',

	# Pick your license as you wish (should match "license" above)
	 'License :: OSI Approved :: MIT License',

	# Specify the Python versions you support here. In particular, ensure
	# that you indicate whether you support Python 2, Python 3 or both.
	'Programming Language :: Python :: 3',
	'Programming Language :: Python :: 3.6',
    ],
    keywords='data storage biology freeze', 
    project_urls={
        'Documentation' : 'http://linkage.io',
        'Source' : 'https://github.com/LinkageIO/HapDab',
        'Tracker' : 'https://github.com/LinkageIO/Minus80/HapDab'
    },



    packages = find_packages(),
    scripts = [
    ],
    ext_modules = [],
    cmdclass = {
    },

    package_data = {
        '':['*.cyx'],    
        'beagle':'include/beagle/beagle.08Jun17.d8b.jar'
    },
    python_requires='>=3.6',
    setup_requires=[
        'setuptools>=18.0',
    ],
    install_requires = [		
        'minus80>=0.2.2',
        'locuspocus>=0.1.2',
        'scipy>=1.1.0',
        'pandas>=0.23.4',
        'numpy>=1.15.0',
        #'pysam>=0.15.0',
        #'cassandra-driver>=3.14.0',
        'aiofiles>=0.4.0'
    ],
    include_package_data=True,

    entry_points='''
        [console_scripts]
        hapdab=hapdab.cli.hapdab:cli
    ''',
    

)
