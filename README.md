This project contains the Python package `puf_quasideal_model`, which provides utilities for exploiting the quasi-ideal PUF model developed by the Group of Electronic Design (GDE) research group at the University of Zaragoza, Spain.

## Requirements
The software is Python-based. A comprehensive list of dependencies can be found in the `requirements.txt` file

## Installation
To get started, follow these steps:

1. Clone this repository and navigate to its directory:
    ```bash
    git clone [] 
    cd []
    ```
2. Create a Python virtual environment:
    ```bash
    python -m venv [venv_name]
    ```
3. Install the dependencies listed in `requirements.txt`:
    ```bash
    python -m pip install -r requirements.txt
    ```
4. Install the `puf_quasideal_model` package:
    ```bash
    python setup.py install
    ```

## Usage

### Simulation

After installation, this package provides two different console scripts: `pqm-simulation` and `pqm-dks-distribution`. These scripts accept various command-line options that control the simulation parameters. The first script generates a *.txt* file containing the simulated intra- and inter-Hamming distance distributions, while the second script generates a *.txt* file containing the Kolmogorov-Smirnov distance distribution obtained from both intra- and inter-Hamming distance binomial fits.

### Analysis

This project includes a Jupyter notebook `Tutorial/Tutorial.ipynb` which guides the user through the entire analysis workflow.