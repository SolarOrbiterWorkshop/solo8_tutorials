#!/bin/bash

# Can be used to start jupyter notebook for the RPW tutorial
# To execute this script:
#       /bin/bash start_tuto.sh
#
# Make sure that dependencies are installed before running it (see README for details).

# Get virtualenv name (rpw_tuto by default)
venv=${1:-rpw_tuto}
# Load it
if [[ -f ~/.virtualenvs/${venv}/bin/activate ]];then
  source ~/.virtualenvs/${venv}/bin/activate
else
  echo "ERROR - Virtualenv ${venv} is not defined!"
  exit 1
fi

# Load NASA CDF environment
if [[ -f ./cdf38_1-dist/bin/definitions.B ]];then
  source ./cdf38_1-dist/bin/definitions.B
fi

# Start jupyter server
jupyter notebook &

deactivate
