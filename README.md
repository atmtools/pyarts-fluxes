# pyarts-fluxes

Version: 0.8

Python module for calculating radiative fluxes with ARTS.
The module is an easy-to-use wrapper to calculate radiative fluxes with ARTS.
The idea behind is to prepare a basic setup so that the user can easily calculate radiative fluxes with ARTS without having to deal with the actual ARTS simulation setup.
Expierenceed users can still access the ARTS workspace and modify it as they like.

Get ARTS (pyarts): https://radiativetransfer.org/getarts/

## Requirements

pyarts-fluxes requires the following Python packages:

- pyarts >=2.6.16
- numpy >=2.0.0

## Installation

To install pyarts-fluxes, clone the repository and run the setup script:

```bash
git clone https://github.com/atmtools/pyarts-fluxes
cd pyarts-fluxes
python -m pip install --user -e .
```

## Usage

See the examples in the examples folder for usage.
