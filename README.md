# Look-Ahead Study

This repository contains all code outside pyABC for analysis and figure
generation for the manuscript
`A Wall-time Minimizing Parallelization Strategy for Sequential Approximate Bayesian Computation`.

## Installation

The study was performed on various Linux systems with Anaconda version TODO,
using python 3.8.5 TODO. The anaconda environment can be downloaded and
installed via:

    wget https://repo.anaconda.com/archive/Anaconda3-2020-07-Linux-x86_64.sh
    bash Anaconda3-2020-07-Linux-x86_64.sh

A virtual environment can be installed via:

    conda create -n la_study python=3.8

and activated via:

    conda activate la_study

Core dependencies of the study are:

    conda install redis==5.0.3
    pip install pyabc==0.10.12

Further dependencies can be found in the respective model subfolders.
In particular, the tumor model requires an installation of
[tumor2d](https://github.com/icb-dcm/tumor2d), and various virus models
require an installation of [morpheus](https://morpheus.gitlab.io).

## Structure

* The `models` folder contains the various test models. Each model consists of
  a model definition, scripts to execute the analysis, and scripts to generate
  the manuscript and supplementary figures.
* The `Batch_pyABC` folder includes the version of Emad's shell scripts to run
  a redis-server based analysis on a cluster environment.
