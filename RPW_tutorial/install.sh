#!/bin/bash

# Install Python dependencies for notebook
# Make sure the NASA CDF software is installed before running this script (see install_cdf.sh)


# Get virtualenv name (rpw_tuto by default)
venv=${1:-rpw_tuto}
# Create it
echo "Creating virtualenv ${venv} ..."
mkdir -p ~/.virtualenvs
python3 -m venv ~/.virtualenvs/${venv}
# Load it
source ~/.virtualenvs/${venv}/bin/activate

# Update pip
python3 -m pip install pip -U

# Install jupyter lab and jupyter notebook
python3 -m pip install jupyterlab
python3 -m pip install notebook

# Add Virtualenv as Python Kernel
echo "Adding virtualenv ${venv} as Python Kernel ..."
ipython kernel install --name "${venv}" --user

# Install poetry
python3 -m pip install poetry

# Install Python dependencies
poetry install --no-root
