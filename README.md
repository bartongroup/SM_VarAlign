# VarAlign
VarAlign is a python module that aggregates human genetic variation over the columns of a multiple sequence alignment.

## Installing

Installing VarAlign (uses Conda)

```sh
# Download
$ git clone https://github.com/stuartmac/VarAlign.git

# set up a conda environment with all requirements
$ cd VarAlign
$ conda env create -f environment.yml
$ source activate varalign-env

# install VarAlign
$ pip install .
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
