# VarAlign
VarAlign is a python module that aggregates human genetic variation over the columns of a multiple sequence alignment.

## Installing

Setting up a conda environment
```sh
$ conda create -n varalign python pip
$ source activate varalign
```

Installing VarAlign

```sh
$ wget https://github.com/stuartmac/VarAlign/new/master.zip -O VarAlign.zip
$ unzip VarAlign.zip

# alternatively
$ git clone https://github.com/stuartmac/VarAlign/new/master.git

# installing requirements
$ cd VarAlign
$ pip install -r requirements.txt
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
