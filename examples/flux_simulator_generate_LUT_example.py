#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:21:58 2024

@author: Manfred Brath

Example script to demonstrate the usage of the flux_simulator_module for generating lookuptables.

"""

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import FluxSimulator as fsm
from pyarts.plots import arts_lookup
from pyarts import arts


# %% generate example atmosphere needed for standard LUT
# This atmosphere is not intended to be fully realistic, but to be simply
# an example for the calculation of LUT.

# set pressure grid
pressure_profile = np.linspace(1000e2, 1e2, 80)

# create water vapor profile
# Water vapor is simply define by a 1st order
# polynomial in log-log space
# log h2o = a + b * log pressure
b = 4
a = -6 - b * 4
logH2O = a + b * np.log10(pressure_profile)
H2O_profile = 10**logH2O

# create temperature profile
# Temperature is simply define by a 1st order
# polynomial of log pressure
# T = a + b * log pressure
# For pressure < 100 hPa, the temperature is set to 200 K
b = 100
a = 200 - b * 4
temperature_profile = a + b * np.log10(pressure_profile)
temperature_profile[pressure_profile < 100e2] = (
    200  # set temperature to 200 K below 100 hPa
)

# CO2 vmr value
CO2 = 400e-6  # [vmr]

# set vmr_field
vmr = np.zeros((2, len(pressure_profile), 1, 1))
vmr[0, :, 0, 0] = H2O_profile
vmr[1, :, 0, 0] = CO2

# plot atmosphere profiles
fig, ax = plt.subplots(1, 2)
ax[0].semilogy(temperature_profile, pressure_profile / 100, label="Temperature")
ax[0].set_xlabel("Temperature [K]")
ax[0].set_ylabel("Pressure [hPa]")
ax[0].invert_yaxis()

ax[1].loglog(H2O_profile, pressure_profile / 100, label="Water vapor")
ax[1].loglog(vmr[1, :, 0, 0], pressure_profile / 100, label="CO$_2$")
ax[1].set_xlabel("vmr")
ax[1].set_ylabel("Pressure [hPa]")
ax[1].invert_yaxis()
ax[1].legend()


# %% generate lut

# set frequency grid
min_wvn = 10
max_wvn = 3210
n_freq_lw = 200
wvn = np.linspace(min_wvn, max_wvn, n_freq_lw)
f_grid_lw = arts.convert.kaycm2freq(wvn)


# setup ARTS
flux_simulator_LUT = fsm.FluxSimulator("TESTLUT_SW")
flux_simulator_LUT.set_frequency_grid(f_grid_lw)
flux_simulator_LUT.set_species(
    ["H2O, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350", "CO2, CO2-CKDMT252"]
)

# Standard LUT
flux_simulator_LUT.get_lookuptable_profile(
    pressure_profile, temperature_profile, vmr, recalc=True
)
LUT = deepcopy(flux_simulator_LUT.ws.abs_lookup.value)

# Wide LUT
flux_simulator_LUT.get_lookuptableWide(recalc=True)
LUT_wide = deepcopy(flux_simulator_LUT.ws.abs_lookup.value)


# plot
arts_lookup.plot_arts_lookup(LUT)
arts_lookup.plot_arts_lookup(LUT_wide)
