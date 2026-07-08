#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:21:58 2024

@author: Manfred Brath

Example script to demonstrate the usage of the flux_simulator_module for batch simulations.

"""

import numpy as np
import FluxSimulator as fsm
from pyarts import arts, xml


# Load test data
atms = xml.load("../atmdata/batch_atmospheres/atms_short.xml")
auxs = xml.load("../atmdata/batch_atmospheres/auxs_short.xml")


# set frequency
min_wavelength_sw = 3e-7  # [m]
max_wavelength_sw = 5e-6  # [m]
n_freq_sw = 200

wvl = np.linspace(min_wavelength_sw, max_wavelength_sw, n_freq_sw)  # [m]
f_grid_sw = arts.convert.wavelen2freq(wvl[::-1])


# some data preparations
surface_altitudes = [aux_i[1] for aux_i in auxs]
surface_temperatures = [aux_i[0] for aux_i in auxs]
geographical_positions = [[aux_i[4], aux_i[5]] for aux_i in auxs]
sun_positions = [[1.495978707e11, 0.0, -120.0] for aux_i in auxs]
refls = [[0.3] for i in range(len(auxs))]


# =============================================================================
# the simulation
# =============================================================================

# setup ARTS
FluxSimulator_batch = fsm.FluxSimulator("BATCH_Test")
FluxSimulator_batch.set_frequency_grid(f_grid_sw)
FluxSimulator_batch.emission = 0
FluxSimulator_batch.gas_scattering = True
FluxSimulator_batch.set_species(
    [
        "H2O, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350",
        "O2-*-1e12-1e99,O2-CIAfunCKDMT100",
        "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
        "CO2, CO2-CKDMT252",
        "O3",
        "O3-XFIT",
    ]
)

FluxSimulator_batch.define_particulate_scatterer(
    "LWC", "pnd_agenda_CGLWC", "MieSpheres_H2O_liquid", ["mass_density"]
)
FluxSimulator_batch.define_particulate_scatterer(
    "RWC", "pnd_agenda_CGRWC", "MieSpheres_H2O_liquid", ["mass_density"]
)
FluxSimulator_batch.define_particulate_scatterer(
    "IWC",
    "pnd_agenda_CGIWC",
    "HexagonalColumn-ModeratelyRough.modified",
    ["mass_density"],
)
FluxSimulator_batch.define_particulate_scatterer(
    "SWC",
    "pnd_agenda_CGSWC_tropic",
    "10-PlateAggregate-ModeratelyRough.modified",
    ["mass_density"],
)
FluxSimulator_batch.define_particulate_scatterer(
    "GWC", "pnd_agenda_CGGWC", "Droxtal-SeverelyRough.modified", ["mass_density"]
)


results = FluxSimulator_batch.flux_simulator_batch(
    atms,
    surface_temperatures,
    surface_altitudes,
    refls,
    geographical_positions,
    sun_positions,
    end_index=5,
)
