# STRONG - Strain Resolution ON Graphs

## Overview

STRONG resolves strains on assembly graphs by resolving variants on core COGs using co-occurrence across multiple samples.

## Prerequisites

The following pieces of software should be installed on your machine before attempting to install STRONG
    - conda (miniconda)
    - cmake, zlib, GNU readline, G++
    
For a standard Ubuntu 16.04 distribution. The above packages would be installed as:

```
    sudo apt-get update
    sudo apt-get -y install libbz2-dev libreadline-dev cmake g++ zlib1g zlib1g-dev
```

Python is also need for the conda install we recommend Python 3.7.
To install miniconda follow the instructions (here)[https://docs.conda.io/en/latest/miniconda.html]:

## Installation

STRONG can be installed anywhere but for the below we assume it will be placed in ~/repos that you 
may have to create in your home dir:
```
cd ~/repos
```

We begin by cloning STRONG recursively:

```
git clone --recurse-submodules https://github.com/chrisquince/STRONG.git
```

STRONG contains [DESMAN](https://github.com/chrisquince/DESMAN) and [BayesPaths](https://github.com/chrisquince/BayesPaths) as submodules.

If you need to update in future:

```
cd STRONG
git submodule foreach git pull origin master
```

We recommend that you first compile the SPAdes and COG tools executables outside of conda:

```
cd ./SPAdes/assembler

./spades_compile.sh

./build_cog_tools.sh 

cd ../..
```

The full list of requirements is listed in the file conda_env.yaml we recommend mamba for install. This can be 
itself installed through conda by:
```
conda install -c conda-forge mamba
```

Then we use mamba to resolve the STRONG environment from with the STRONG home directory:

```
mamba env create -f conda_env.yaml
```

This should take 5 - 10 minutes with mamba.


Once the STRONG environment has been installed activate it with the following command :

```
conda activate STRONG
```


It is also necessary to install the BayesPaths executable with the STRONG conda:

```
cd BayesPaths
python ./setup.py install
```

BayesPaths use precompiled executables in the runfg_source directory. These are only compatible 
with Linux x86-64 and on other platforms they will require compilation from source see
the [BayesPaths repo](https://github.com/chrisquince/BayesPaths]) for details. 


Finally we will also need a version of the COG database installed. We make this available for download 
and again we recommend placing it in a directory ~/Database but it could be placed anywhere:

```
mkdir Database
wget https://strongtest.s3.climb.ac.uk/rpsblast_cog_db.tar.gz
tar -xvzf rpsblast_cog_db.tar.gz
```

## Native installation (Not supported yet)

STRONG has a lot of required software, this is an attempt to demonstrate how to install all of them to avoid the 
conda recipe above. 


## Quick start

First we will download a fairly simple synthetic test data set from known microbial strains into another directory 
~/STRONG_Runs that we will use for STRONG output:

```
mkdir ~/STRONG_Runs
cd  ~/STRONG_Runs
wget https://strongtest.s3.climb.ac.uk/Test.tar.gz
tar -xvzf Test.tar.gz
```

We are now ready to run STRONG from within the COG_pipe directory. Using the following command:

```
cd ~/repos/STRONG/COG_pipe
python3 ./start.py --config config.yaml ~/STRONG_Runs/TestResults --threads 32 --dryrun
```

Optionally pass snakemake parameters with the -s option e.g. '-s --dryrun'

## Config file

```
# ------ Samples ------
samples: '*' # specify a list samples to use or '*' to use all samples

# ------ Resources ------
threads : 8 # single task nb threads

# ------ Assembly parameters ------ 
data: /mnt/gpfs/Hackathon/Test  # path to data folder

# ----- Annotation database -----
cog_database: /home/sebr/seb/Database/rpsblast_cog_db/Cog # COG database

# ----- SPAdes tools dependency -----
soft: /home/sergei/cog_tools

# ----- Binning parameters ------
concoct_contig_size: 1000
read_length: 150
assembly: 
    assembler: spades
    k: [77]
    mem: 2000
    threads: 24
    dir: /home/sergei/cog_tools/spades/bin

# ----- BayesPaths parameters ------
bayespaths:
    nb_strains: 16

# ----- DESMAN parameters ------
desman:
    execution: 1
    nb_haplotypes: 10
    nb_repeat: 5
    min_cov: 1

# -----  MAGAnalysis ------
maganalysis: 
    execution: 0

# -----  Evaluation ------
evaluation:
    execution: 1
    genomes: "/mnt/gpfs/Hackathon/Test/Eval" # path to refferences genomes 
```

## Pipeline
# Assembly and binning 
![alt tag](./Figures/Dag_rules1.png)
# BayesAGraohsSVA
![alt tag](./Figures/Dag_rules2.png)
# Desman 
to be uploaded
# MAGanalysis
![alt tag](./Figures/Dag_rules5.png)
# Evaluation
![alt tag](./Figures/Dag_rules6.png)