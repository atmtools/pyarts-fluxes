#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:21:58 2024

@author: Manfred Brath

Example script to demonstrate the usage of the flux_simulator_module for batch
spectral radiance simulations.

"""

import numpy as np
import FluxSimulator as fsm
from pyarts import arts, xml


# Load test data
atms = xml.load("../atmdata/batch_atmospheres/atms_short.xml")
auxs = xml.load("../atmdata/batch_atmospheres/auxs_short.xml")
# Load with more profiles than needed for demoing start/end_index
# atms = xml.load("../atmdata/batch_atmospheres/atms.xml")
# auxs = xml.load("../atmdata/batch_atmospheres/auxs.xml")


# set frequency
min_wavelength_sw = 3e-7  # [m]
max_wavelength_sw = 5e-6  # [m]
n_freq_sw = 200

wvl = np.linspace(min_wavelength_sw, max_wavelength_sw, n_freq_sw)  # [m]
f_grid_sw = arts.convert.wavelen2freq(wvl[::-1])


# some data preparations
surface_altitudes = [aux_i[1] for aux_i in auxs]
surface_tempratures = [aux_i[0] for aux_i in auxs]
geographical_positions = [[aux_i[4], aux_i[5]] for aux_i in auxs]
sun_positions = [[1.495978707e11, 0.0, -120.0] for aux_i in auxs]
refls = [[0.3] for i in range(len(auxs))]


# =============================================================================
# the simulation
# =============================================================================

# setup ARTS
# the setup_name determines the LUT cache path (cache/<setup_name>); if no
# LUT.xml is found there, a new one is calculated
FluxSimulator_batch = fsm.FluxSimulator("BATCH_Radiance_Test_v2")
FluxSimulator_batch.set_frequency_grid(f_grid_sw)
FluxSimulator_batch.emission = 0
FluxSimulator_batch.gas_scattering = True
FluxSimulator_batch.set_species(
    [
        "H2O, H2O-SelfContCKDMT430, H2O-ForeignContCKDMT430",
        "O2-*-1e12-1e99,O2-CIAfunCKDMT100",
        "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
        "CO2, CO2-CKDMT252",
        "O3",
        "O3-XFIT",
    ]
)


# the angular grid is shared by all profiles of the batch
# N_za defines the number of zenith angles, N_aa the number of azimuth angles
N_za = 20
N_aa = 41

results = FluxSimulator_batch.radiance_simulator_batch(
    atms,
    surface_tempratures,
    surface_altitudes,
    refls,
    geographical_positions,
    sun_positions,
    N_za=N_za,
    N_aa=N_aa,
    end_index=5,
)


# =============================================================================
# show the results
# =============================================================================

print("result keys:", list(results.keys()))
print(
    "radiance (clearsky) for",
    len(results["array_of_spectral_radiance_clearsky"]),
    "profiles"
)
print(
    "radiance array shape of profile 0 (freq x altitude x zen x azi):",
    results["array_of_spectral_radiance_clearsky"][0].shape,
)
print("zenith angle grid [deg]:", results["zenith_angle"])
print("azimuth angle grid [deg]:", results["azimuth_angle"])
print("altitude grid of profile 0 [m]:", results["array_of_altitude"][0])
