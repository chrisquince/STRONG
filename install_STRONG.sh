#!/bin/bash

# STRONG folder
{
set -e
STRONG_dir=$(dirname $(dirname "$(pwd)/$0"))
LOG=$STRONG_dir/install.log

# SPAdes install
echo "SPAdes install ongoing"
cd $STRONG_dir/SPAdes/assembler
./build_cog_tools.sh &>$LOG && echo "SPAdes install succesfully"

echo "creating STRONG env"
# mamba
echo "creating STRONG env" >>$LOG
mamba -h &>>$LOG || conda install -y -c conda-forge mamba &>>$LOG

# conda env 
mamba env create -f $STRONG_dir/conda_env.yaml &>>$LOG

CONDA_PATH=$(dirname $(dirname $(which conda)))
echo $CONDA_PATH
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate STRONG


# install bayespath
echo "Installing BayesPaths"
cd  $STRONG_dir/BayesPaths
python ./setup.py install &>>$LOG

# install desman
echo "Installing DESMAN"
cd $STRONG_dir/DESMAN
python ./setup.py install &>>$LOG

# correct R lapack library
ln -fs $CONDA_PREFIX/lib/R/modules/lapack.so $CONDA_PREFIX/lib/R/modules/libRlapack.so

# concoct_refine
PATH_concoctR=$(which concoct_refine)
sed -i 's/values/to_numpy/g' $PATH_concoctR
sed -i 's/as_matrix/to_numpy/g' $PATH_concoctR
sed -i 's/int(NK), args.seed, args.threads)/ int(NK), args.seed, args.threads, 500)/g' $PATH_concoctR
} || { echo -e "\033[0;31msomething sinister just happened, check the log at :\n$LOG\033[0m"; exit 1;}
# check install
$STRONG_dir/SnakeNest/scripts/check_on_dependencies.py

