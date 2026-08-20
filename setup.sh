#!/bin/bash

# Download AF2 parameters
mkdir -p ./src/params
cd ./src/params
echo "Downloading AlphaFold parameters..."
wget -c https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
tar -xf alphafold_params_2021-07-14.tar
rm alphafold_params_2021-07-14.tar
cd ../../

#Python packages
CONDA_CMD=$(command -v mamba || command -v conda) # pick mamba if installed, otherwise default to conda
echo "Creating environment using $CONDA_CMD..."
$CONDA_CMD env create -f environment.yml 2>/dev/null || $CONDA_CMD env update -f environment.yml --prune
CONDA_BASE=$($CONDA_CMD info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate EBM

## HHblits
echo "Setting up HH-suite..."
git clone https://github.com/soedinglab/hh-suite.git
mkdir -p hh-suite/build && cd hh-suite/build
cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_INSTALL_PREFIX=. ..
make -j 4
make install
cd ../..
conda deactivate

### Download uniclust30_2018_08
echo "Downloading Uniclust30 database..."
cd data
wget -c https://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz
tar -zxvf uniclust30_2018_08_hhsuite.tar.gz
rm -f uniclust30_2018_08_hhsuite.tar.gz
cd ..

echo "Setup completed!"
