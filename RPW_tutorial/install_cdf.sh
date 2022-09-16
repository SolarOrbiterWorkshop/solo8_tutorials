#!/bin/bash

# Download and install NASA CDF software

if ! [[ -f cdf38_1-dist-all.tar.gz ]];then
  CDF_URL="https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/cdf38_1-dist-all.tar.gz"
  wget $CDF_URL
fi
tar -xf cdf38_1-dist-all.tar.gz
cd ./cdf38_1-dist
make OS=linux ENV=gnu CURSES=yes FORTRAN=no UCOPTIONS=-O2 SHARED=yes all
make install
cd ..
