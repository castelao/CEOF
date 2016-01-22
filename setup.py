#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Read and comply:
#    https://python-packaging-user-guide.readthedocs.org/en/latest/tutorial.html#creating-your-own-project

here = path.abspath(path.dirname(__file__))

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
        'numpy',
        'scipy'
]

version_file = open(path.join(here, 'VERSION'))
version = version_file.read().strip()

setup(
    name='CEOF',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version=version,

    description='Apply the Complex Empirical Orthogonal Function Analysis',
    long_description=readme + '\n\n' + history,

    # The project's main homepage.
    url='https://github.com/castelao/pyCEOF',

    # Author details
    author='Guilherme Castelao',
    author_email='guilherme@castelao.net',

    install_requires=requirements,
    license='3-clause BSD',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: BSD License',

        'Programming Language :: Python :: 2.7',
    ],

    keywords='EOF, PCA, complex empirical orthogonal function, complex principal component analysis',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=[
        'ceof',
    ],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages.
    # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)
