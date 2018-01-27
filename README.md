# VarAlign
VarAlign is a python module that aggregates human genetic variation over the columns of a multiple sequence alignment.

## Installing

Installing VarAlign (uses Conda)

```sh
# Download
$ git clone https://github.com/stuartmac/VarAlign.git

# Set up conda environment with all requirements
$ cd VarAlign
$ conda env create -f environment.yml
$ source activate varalign-env

# Install VarAlign
$ pip install .

# Check it works
$ cd tests/data
$ python ../../varalign/align_variants.py sample_swissprot_PF00001.18_full.sto
```

### Enabling structural analysis (requires ProIntVar and Arpeggio)

Arpeggio needs to go in a seperate environment because it requires Python 2.
```
# Setup environment with arpeggio dependencies
$ conda create -n arpeggio python=2 pip numpy biopython openbabel
$ source activate arpeggio

# Get patched Arpeggio
$ git clone https://bitbucket.org/biomadeira/arpeggio
$ cd arpeggio
$ python arpeggio.py -h
```

Install and configure ProIntVar. (NB. ProIntVar requires Python 3.)
```
# Install ProIntVar (https://github.com/bartongroup/ProIntVar)
$ source activate varalign-env-py3
$ git clone https://github.com/bartongroup/ProIntVar-Core.git
$ cd ProIntVar-Core
$ pip install .

# Configure ProIntVar
$ ProIntVar-Core-config-setup prointvar_config.ini

# *** Edit the following values in prointvar_config.ini ***
# arpeggio_bin = /path/to/arpeggio/arpeggio.py
# python_exe = /path/to/anaconda/envs/arpeggio/bin/python
# python_path = /path/to/anaconda/envs/arpeggio/python/lib/site-packages/

$ ProIntVar-Core-config-load prointvar_config.ini

# Check it works (after running `align_variants.py` as above)
$ cd path/to/VarAlign/tests/data
$ python ../../varalign/prointvar_analysis.py sample_swissprot_PF00001.18_full.sto
```


## Configuration

VarAlign uses a configuration file to set key paths and parameters that need to be set before you run. Priority is given to a config file
that is present in the execution directory, letting you keep parameters beside results, but a global config file can also be used.

Setting up a config file in the working directory
```sh
$ cd /path/to/desired/working/dir/

# Get a copy of the template config.ini file shipped with ProIntVar
$ cp /path/to/VarAlign/varalign/config.txt ./

# Edit the settings as you require

# Testing that the new values are correctly loaded by ProIntVar
$ python
>>> from varalign.config import defaults
>>> defaults.gnomad
'./sample_swissprot_PF00001.18_full.vcf.gz'
```
