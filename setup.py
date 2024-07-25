#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup

from varalign import __license__, __version__


def gather_dependencies():
    with open("requirements.txt", "r") as f_in:
        return [
            line.strip() for line in f_in if line and not line.startswith("#")
        ]


DEPENDENCIES = gather_dependencies()


def read_readme():
    with open("README.md", "r", encoding="utf-8") as f:
        return f.read()


setup(
    name="VarAlign",
    version=__version__,
    packages=find_packages(
        include=["varalign", "varalign.*"], exclude=["tests", "tests.*"]
    ),
    package_data={
        "varalign": [
            "config.txt",
            "data/*",
            "lib/aacon/*",
        ],
    },
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
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.6",
)
