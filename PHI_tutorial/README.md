# SO/PHI Data Tutorial @ Solo8

This is a repo that contains the SO/PHI tutorial for the Data Analysis Day at the Solar Orbiter 8th Workshop in Belfast UK.

The tutorial will be run through the jupyter notebook: `belfast_data_tutorial.ipnb`

## Execution

A zip of the data file folder can be downloaded here (570MB): https://owncloud.gwdg.de/index.php/s/WzTfaBYacfVSNpZ

Password: solo8

The necessary python packages can be found in the `requirements.txt` file.

Quick-start with conda:

First `cd` into the PHI_tutorial folder
```bash=
conda config --add channels conda-forge
conda config --set channel_priority strict
conda config --set channel_priority strict
conda create --name phi_tutorial_env --file requirements.txt
conda activate phi_tutorial_env
python -m ipykernel install --user --name phi_tutorial_env
jupyter notebook
```
Please then make sure that the kernel is set to `phi_tutorial_env`

## Author

### Jonas Sinjan (PhD Student at Max Planck Institute for Solar System Research, Goettingen, Germany)

#### For questions you can reach me by email: sinjan@mps.mpg.de

#### *with many thanks to Daniele Calchetti, Gherardo Valori and the SO/PHI team*

<img src="./philogo-1.png" width="220" align="left"/>


