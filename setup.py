#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
from setuptools import setup
from setuptools import find_packages

from varalign import __version__, __license__


def gather_dependencies():
    with open('requirements.txt', 'r') as f_in:
        return [l for l in f_in.read().rsplit(os.linesep)
                if l and not l.startswith("#")]
DEPENDENCIES = gather_dependencies()


setup(
    name='VarAlign',
    version=__version__,
    packages=find_packages(exclude=['tests', 'tests.*']),
    # data_files=[('config', ['varalign/config.txt'])],
    package_data={'varalign': ['config.txt']},
    include_package_data=True,
    scripts=['bin/filter_swiss'],
    install_requires=DEPENDENCIES,
    url='https://github.com/stuartmac/VarAlign/',
    license=__license__,
    author='Stuart MacGowan',
    author_email='s.macgowan@dundee.ac.uk',
    description='This package is used to map and aggregate variants in multiple sequence alignments.'
)
