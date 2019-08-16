# STRONG - Strain Resolution ON Graphs

## Overview

STRONG resolves strain on assembly graphs by resolving variants on core COGs using co-occurrence across multiple samples.

## Installation

Requires recursive cloning:

```
git clone --recurse-submodules https://github.com/chrisquince/STRONG.git
```

To update

```
git submodule foreach git pull origin master
```

## Quick start

Run from within the COG_pipe directory. Using the following command:

```
python3 ./start.py --config config.yaml output_dir --threads 32
```

Optionally pass snakemake parameters e.g. '--dryrun'

## Config file

```
# ------ Resssources ------ 
threads : 8 # single task nb threads
# ------ Assembly parameters ------ 
data: /mnt/gpfs/Hackathon/Test  # path to data folder
# ---- Annotation database -----
cog_database: /home/sebr/seb/Database/rpsblast_cog_db/Cog # COG database 
# ---- Spade tools dependency -----
soft: /home/sergei/cog_tools2
# ----Binning parameters ------
concoct_contig_size: 1000
read_length: 150
assembly: 
    assembler: spades
    k: [77]
    mem: 2000
    threads: 24
    dir: /home/sergei/cog_tools2/spades/bin
    groups: ['*'] # specify a group of sample to coassemble
# ---- Bayespaths parameters ------
bayespaths:
    nb_strains: 16
# ---- Desman parameters ------
desman:
    execution: 1
    nb_haplotypes: 10
    nb_repeat: 5
    min_cov: 1
# ---- Maganalysis ------
maganalysis: 
    execution: 0
```

## Pipeline

![alt tag](./Figures/Dag1.png)

![alt tag](./Figures/Dag2.pdf)
