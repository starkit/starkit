#!/bin/bash

# create local env
conda create -n python36 python=3.6 --yes

source activate python36

conda install -c conda-forge doctr --yes

doctr deploy . --built-docs docs/_build/html 