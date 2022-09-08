FROM python:3.8
ENV PYTHONUNBUFFERED=1

# Use bash as default shell
ENV SHELL /bin/bash

# install pip and poetry
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache poetry

# Make sure to have notebook jupyterlab installed
RUN python3 -m pip install --no-cache-dir notebook jupyterlab

# install libs needed to install nasa cdf
RUN apt-get update && apt-get install -y gfortran ncurses-dev

# create user with a home directory
ARG NB_USER="rpw_user"
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --shell /bin/bash \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}

# install cdf
ENV CDF_URL="https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/cdf38_1-dist-all.tar.gz"
RUN wget $CDF_URL
RUN tar -xf cdf38_1-dist-all.tar.gz
WORKDIR ${HOME}/cdf38_1-dist
RUN make OS=linux ENV=gnu CURSES=yes FORTRAN=no UCOPTIONS=-O2 SHARED=yes all
RUN make install

WORKDIR ${HOME}
USER root

# avoid virtualenv creation
RUN poetry config virtualenvs.create false --local
RUN poetry install --no-root

RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

# Add CDF source file to .bashrc
RUN echo "source $HOME/cdf38_1-dist/bin/definitions.B" >> ~/.bashrc
