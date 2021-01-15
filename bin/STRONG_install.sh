#!/bin/bash

# STRONG folder
STRONG_dir= $dirname($(dirname $pwd/$0))

# SPAdes install
$STRONG_dir/SPAdes/assembler/build_cog_tools.sh && echo "SPAdes install succesfully"

# mamba
mamba -h || conda install -c conda-forge mamba

# conda env 
mamba env create -f $STRONG_dir/conda_env.yaml

conda activate STRONG

# install bayespath
python $STRONG_dir/BayesPaths/setup.py install

# install desman
python $STRONG_dir/DESMAN/setup.py install

# correct R lapack library
ln -s $CONDA_PREFIX/lib/R/modules/lapack.so $CONDA_PREFIX/lib/R/modules/libRlapack.so

# concoct_refine
PATH_concoctR=$(which concoct_refine)
sed -i 's/values/to_numpy/g' $PATH_concoctR
sed -i 's/as_matrix/to_numpy/g' $PATH_concoctR
sed -i 's/int(NK), args.seed, args.threads)/ int(NK), args.seed, args.threads, 500)/g' $PATH_concoctR



