# VarAlign
VarAlign is a python module that aggregates human genetic variation over the columns of a multiple sequence alignment.

## Table of Contents

- [Installing](#installing)
  - [Enabling structural analysis](#enabling-structural-analysis)
- [Configuration](#configuration)

## Installing

Installing VarAlign (uses Conda)

```sh
# Download
$ git clone https://github.com/stuartmac/VarAlign.git

# Set up conda environment with all requirements
# You may need to add Bioconda channels:
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$
$ cd VarAlign
$ conda env create -f environment.yml
$ source activate varalign-env-py3

# Install VarAlign
$ pip install .

# Run tests (some will probably fail, in particular ProIntVar isn't installed yet!)
$ python -m unittest discover tests
```

### Enabling structural analysis
*\*requires ProIntVar and Arpeggio*

Arpeggio needs to go in a seperate environment because it requires Python 2.
```
# Setup environment with arpeggio dependencies
$ conda create -n arpeggio python=2 pip numpy biopython openbabel
$ source activate arpeggio

# Get patched Arpeggio
# Remember to leave the VarAlign folder if you're following this in order!
$ git clone https://bitbucket.org/biomadeira/arpeggio
$ cd arpeggio
$ python arpeggio.py -h
```

Install and configure ProIntVar. (NB. ProIntVar requires Python 3.)
```
# Install ProIntVar (https://github.com/bartongroup/ProIntVar)
$ source activate varalign-env-py3
$ git clone https://github.com/bartongroup/ProIntVar.git
$ cd ProIntVar

# Patch ProIntVar and install
$ git apply /path/to/VarAlign/ProIntVar.patch
$ pip install .

# Configure ProIntVar
$ ProIntVar-config-setup prointvar_config.ini

# *** Edit the following values in prointvar_config.ini ***
# arpeggio_bin = /path/to/arpeggio/arpeggio.py
# python_exe = /path/to/anaconda/envs/arpeggio/bin/python
# python_path = /path/to/anaconda/envs/arpeggio/python/lib/site-packages/

$ ProIntVar-config-load prointvar_config.ini

# Check it works (rerun VarAlign tests)
# If you ran the tests earlier you may see a FileExists error for .../VarAlign/tests/tmp, delete this and try again.
$ cd path/to/VarAlign/
$ python -m unittest discover tests
```


## Configuration

VarAlign uses a configuration file to set key paths and parameters that need to be set before you run. Priority is given to a config file
that is present in the execution directory, letting you keep parameters beside results, but a global config file can also be used.

Setting up a config file in the working directory
```sh
$ cd /path/to/desired/working/dir/

# Get a copy of the template config file shipped with VarAlign
$ cp /path/to/VarAlign/varalign/config.txt ./

# Edit the settings as you require

# Testing that the new values are correctly loaded by VarAlign
$ python
>>> from varalign.config import defaults
>>> defaults.gnomad
'./sample_swissprot_PF00001.18_full.vcf.gz'
```


## Run the pipeline
I recommend you download a Pfam alignment that has at least a few human sequences and then try:

`varalign --species HUMAN <YOUR_ALIGNMENT>`
