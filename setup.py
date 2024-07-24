#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

from varalign import __version__, __license__


def gather_dependencies():
    with open("requirements.txt", "r") as f_in:
        return [line.strip() for line in f_in if line and not line.startswith("#")]


DEPENDENCIES = gather_dependencies()

setup(
    name="VarAlign",
    version=__version__,
    packages=find_packages(exclude=["tests", "tests.*"]),
    package_data={"varalign": ["config.txt"]},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "filter_swiss=varalign.cli.filter_swiss:main",
            "varalign=varalign.cli.varalign:main",
            "index_pfam=varalign.cli.index_pfam:main",
        ],
    },
    install_requires=DEPENDENCIES,
    url="https://github.com/stuartmac/VarAlign/",
    license=__license__,
    author="Stuart MacGowan",
    author_email="s.macgowan@dundee.ac.uk",
    description="This package is used to map and aggregate variants in multiple sequence alignments.",
)
