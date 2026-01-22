#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:33:01 2022

@author: Manfred Brath
"""

from pyarts.workspace import arts_agenda

# =============================================================================
# gas scattering agendas
# =============================================================================


def create_gas_scattering_agendas_in_WS(ws):

    gas_scattering_agenda_list = [
        "gas_scattering_agenda__Rayleigh",
        "gas_scattering_agenda__IsotropicConstant",
    ]

    for gas_scattering_agenda in gas_scattering_agenda_list:
        ws.AgendaCreate(gas_scattering_agenda)

    return ws, gas_scattering_agenda_list


# gas scattering agenda
@arts_agenda
def gas_scattering_agenda__Rayleigh(ws):
    ws.Ignore(ws.rtp_vmr)
    ws.gas_scattering_coefAirSimple()
    ws.gas_scattering_matRayleigh()


@arts_agenda
def gas_scattering_agenda__IsotropicConstant(ws):
    ws.Ignore(ws.rtp_vmr)
    ws.gas_scattering_coefXsecConst(ConstXsec=4.65e-31)  # 31
    ws.gas_scattering_matIsotropic()


# =============================================================================
#  surface agenda definitions
# =============================================================================


def create_surface_agendas_in_WS(ws):

    iy_surface_agenda_list = [
        "iy_surface_agenda__FlatUserReflectivity",
        "iy_surface_agenda__FlatFresnel",
        "iy_surface_agenda__Lambertian",
    ]

    for iy_surface_agenda in iy_surface_agenda_list:
        ws.AgendaCreate(iy_surface_agenda)

    return ws, iy_surface_agenda_list


# surface scattering agenda for a flat surface with userdefined reflectivity
@arts_agenda
def iy_surface_agenda__FlatUserReflectivity(ws):
    ws.iySurfaceInit()
    ws.Ignore(ws.dsurface_rmatrix_dx)
    ws.Ignore(ws.dsurface_emission_dx)
    ws.iySurfaceFlatReflectivity()
    ws.iySurfaceFlatReflectivityDirect()


# surface scattering agenda for a flat surface with Fresnel reflectivity
@arts_agenda
def iy_surface_agenda__FlatFresnel(ws):
    ws.iySurfaceInit()
    ws.Ignore(ws.dsurface_rmatrix_dx)
    ws.Ignore(ws.dsurface_emission_dx)
    ws.iySurfaceFlatRefractiveIndex()
    ws.iySurfaceFlatRefractiveIndexDirect()


# surface scattering agenda for a lambertian surface with user defined reflectivity
@arts_agenda
def iy_surface_agenda__Lambertian(ws):
    ws.iySurfaceInit()
    ws.Ignore(ws.dsurface_rmatrix_dx)
    ws.Ignore(ws.dsurface_emission_dx)
    ws.iySurfaceLambertian()
    ws.iySurfaceLambertianDirect()


def set_iy_surface_agendas(ws):
    ws.iy_surface_agenda__FlatUserReflectivity = iy_surface_agenda__FlatUserReflectivity
    ws.iy_surface_agenda__FlatFresnel = iy_surface_agenda__FlatFresnel
    ws.iy_surface_agenda__Lambertian = iy_surface_agenda__Lambertian

    return ws


# =============================================================================
# pnd_agendas
# =============================================================================


#### Seifert and Beheng, 2006 two moment scheme  ####
# ==========================================================================


def create_pnd_agendas_SB06_in_WS(ws):

    pnd_agenda_list = [
        "pnd_agenda_SB06LWC",
        "pnd_agenda_SB06IWC",
        "pnd_agenda_SB06RWC",
        "pnd_agenda_SB06SWC",
        "pnd_agenda_SB06GWC",
        "pnd_agenda_SB06HWC",
    ]

    for pnd_agenda in pnd_agenda_list:
        ws.AgendaCreate(pnd_agenda)

    return ws, pnd_agenda_list


# LWC
@arts_agenda
def pnd_agenda_SB06LWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "mass", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdSeifertBeheng06(hydrometeor_type="cloud_water", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# IWC
@arts_agenda
def pnd_agenda_SB06IWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "mass", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdSeifertBeheng06(hydrometeor_type="cloud_ice", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Rain
@arts_agenda
def pnd_agenda_SB06RWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "mass", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdSeifertBeheng06(hydrometeor_type="rain", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Snow
@arts_agenda
def pnd_agenda_SB06SWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "mass", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdSeifertBeheng06(hydrometeor_type="snow", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Graupel
@arts_agenda
def pnd_agenda_SB06GWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "mass", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdSeifertBeheng06(hydrometeor_type="graupel", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Hail
@arts_agenda
def pnd_agenda_SB06HWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "mass", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdSeifertBeheng06(hydrometeor_type="hail", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


def set_pnd_agendas_SB06(ws):
    ws.pnd_agenda_SB06LWC = pnd_agenda_SB06LWC
    ws.pnd_agenda_SB06IWC = pnd_agenda_SB06IWC
    ws.pnd_agenda_SB06RWC = pnd_agenda_SB06RWC
    ws.pnd_agenda_SB06SWC = pnd_agenda_SB06SWC
    ws.pnd_agenda_SB06GWC = pnd_agenda_SB06GWC
    ws.pnd_agenda_SB06HWC = pnd_agenda_SB06HWC

    return ws


#### Milbrandt and Yau, 2005 two moment scheme  ####
# ==========================================================================


def create_pnd_agendas_MY05_in_WS(ws):

    pnd_agenda_list = [
        "pnd_agenda_MY05LWC",
        "pnd_agenda_MY05IWC",
        "pnd_agenda_MY05RWC",
        "pnd_agenda_MY05SWC",
        "pnd_agenda_MY05GWC",
        "pnd_agenda_MY05HWC",
    ]

    for pnd_agenda in pnd_agenda_list:
        ws.AgendaCreate(pnd_agenda)

    return ws, pnd_agenda_list


# LWC
@arts_agenda
def pnd_agenda_MY05LWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdMilbrandtYau05(hydrometeor_type="cloud_water", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# IWC
@arts_agenda
def pnd_agenda_MY05IWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdMilbrandtYau05(hydrometeor_type="cloud_ice", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Rain
@arts_agenda
def pnd_agenda_MY05RWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdMilbrandtYau05(hydrometeor_type="rain", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Snow
@arts_agenda
def pnd_agenda_MY05SWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdMilbrandtYau05(hydrometeor_type="snow", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Graupel
@arts_agenda
def pnd_agenda_MY05GWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdMilbrandtYau05(hydrometeor_type="graupel", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


# Hail
@arts_agenda
def pnd_agenda_MY05HWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.psdMilbrandtYau05(hydrometeor_type="hail", t_min=0, t_max=999, picky=0)
    ws.pndFromPsdBasic()


def set_pnd_agendas_MY05(ws):
    ws.pnd_agenda_MY05LWC = pnd_agenda_MY05LWC
    ws.pnd_agenda_MY05IWC = pnd_agenda_MY05IWC
    ws.pnd_agenda_MY05RWC = pnd_agenda_MY05RWC
    ws.pnd_agenda_MY05SWC = pnd_agenda_MY05SWC
    ws.pnd_agenda_MY05GWC = pnd_agenda_MY05GWC
    ws.pnd_agenda_MY05HWC = pnd_agenda_MY05HWC

    return ws


#### Cosmo Style Scheme  ####
# ==========================================================================
# For all ice species we assume hard ice spheres
"""
The psds are not exactly the ones from the cosmo.
Graupel and snow should be similar to COSMO graupel.
Cloud ice and cloud liquid are taken from Geer and Baordo (2014).
RWC is simply one that was alreay existing in ARTS."""


def create_pnd_agendas_CG_in_WS(ws):

    pnd_agenda_list = [
        "pnd_agenda_CGLWC",
        "pnd_agenda_CGIWC",
        "pnd_agenda_CGRWC",
        "pnd_agenda_CGSWC_tropic",
        "pnd_agenda_CGSWC_midlatitude",
        "pnd_agenda_CGGWC",
    ]

    for pnd_agenda in pnd_agenda_list:
        ws.AgendaCreate(pnd_agenda)

    return ws, pnd_agenda_list


# LWC
@arts_agenda
def pnd_agenda_CGLWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)    
    ws.psdModifiedGammaMass(n0=-999, mu=2, la=2.13e5, ga=1, t_min=0, t_max=400)
    ws.pndFromPsdBasic()


# IWC
@arts_agenda
def pnd_agenda_CGIWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)
    ws.psdModifiedGammaMass(n0=-999, mu=2, la=2.05e5, ga=1, t_min=0, t_max=400)
    ws.pndFromPsdBasic()


# RWC
@arts_agenda
def pnd_agenda_CGRWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)
    ws.psdAbelBoutle12(t_min=0, t_max=400)
    ws.pndFromPsdBasic()


# SWC tropic
@arts_agenda
def pnd_agenda_CGSWC_tropic(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)
    ws.psdFieldEtAl07(regime="TR", t_min=0, t_max=400)
    ws.pndFromPsdBasic()


# SWC tropic
@arts_agenda
def pnd_agenda_CGSWC_midlatitude(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)
    ws.psdFieldEtAl07(regime="ML", t_min=0, t_max=400)
    ws.pndFromPsdBasic()


# GWC
@arts_agenda
def pnd_agenda_CGGWC(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)
    ws.psdModifiedGammaMass(n0=4e6, mu=0, la=-999, ga=1, t_min=0, t_max=400)
    ws.pndFromPsdBasic()


def set_pnd_agendas_CG(ws):
    ws.pnd_agenda_CGLWC = pnd_agenda_CGLWC
    ws.pnd_agenda_CGIWC = pnd_agenda_CGIWC
    ws.pnd_agenda_CGRWC = pnd_agenda_CGRWC
    ws.pnd_agenda_CGSWC_tropic = pnd_agenda_CGSWC_tropic
    ws.pnd_agenda_CGSWC_midlatitude = pnd_agenda_CGSWC_midlatitude
    ws.pnd_agenda_CGGWC = pnd_agenda_CGGWC

    return ws

# =============================================================================
# additional pnd agendas

def create_additional_pnd_agendas_in_WS(ws):

    additional_pnd_agenda_list = [
        "pnd_agenda_Delanoe14",
        "pnd_agenda_AbelBoutle12",
        "pnd_agenda_FieldEtAl07ML",
        "pnd_agenda_FieldEtAl07TR",
        "pnd_agenda_McFarquaharHeymsfield97",
        "pnd_agenda_AerosolApproximatelyLogNormal",
    ]

    for pnd_agenda in additional_pnd_agenda_list:
        ws.AgendaCreate(pnd_agenda)

    return ws, additional_pnd_agenda_list


@arts_agenda
def pnd_agenda_Delanoe14(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_volume_equ", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dveq", species_index=ws.agenda_array_index)    
    ws.psdDelanoeEtAl14( t_min=0, t_max=400, n0star=-999.0, Dm=-999.0 )
    ws.pndFromPsdBasic()

@arts_agenda
def pnd_agenda_AbelBoutle12(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)    
    ws.psdAbelBoutle12(t_min=0, t_max=400)
    ws.pndFromPsdBasic()

# Field et al 2007
@arts_agenda
def pnd_agenda_FieldEtAl07ML(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)    
    ws.psdFieldEtAl07(regime="ML", t_min=0, t_max=400)
    ws.pndFromPsdBasic()    

@arts_agenda
def pnd_agenda_FieldEtAl07TR(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)    
    ws.psdFieldEtAl07(regime="TR", t_min=0, t_max=400)
    ws.pndFromPsdBasic()

@arts_agenda
def pnd_agenda_McFarquaharHeymsfield97(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_volume_equ", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dveq", species_index=ws.agenda_array_index)    
    ws.psdMcFarquaharHeymsfield97(t_min=0, t_max=400)
    ws.pndFromPsdBasic()

@arts_agenda
def pnd_agenda_AerosolApproximatelyLogNormal(ws):
    ws.ExtractFromMetaSingleScatSpecies(
        ws.pnd_size_grid, ws.scat_meta, "diameter_max", ws.agenda_array_index
    )
    ws.Copy(ws.psd_size_grid, ws.pnd_size_grid)
    ws.ScatSpeciesSizeMassInfo(x_unit="dmax", species_index=ws.agenda_array_index)
    ws.psdModifiedGammaMass(n0=-999, mu=4.527, la=2.82e10, ga=1, t_min=0, t_max=400)
    ws.pndFromPsdBasic()



def set_additional_pnd_agendas(ws):
    ws.pnd_agenda_Delanoe14 = pnd_agenda_Delanoe14
    ws.pnd_agenda_AbelBoutle12 = pnd_agenda_AbelBoutle12
    ws.pnd_agenda_FieldEtAl07ML = pnd_agenda_FieldEtAl07ML
    ws.pnd_agenda_FieldEtAl07TR = pnd_agenda_FieldEtAl07TR
    ws.pnd_agenda_McFarquaharHeymsfield97 = pnd_agenda_McFarquaharHeymsfield97
    ws.pnd_agenda_AerosolApproximatelyLogNormal = pnd_agenda_AerosolApproximatelyLogNormal

    return ws


# =============================================================================
# dobatch calculation agendas
# =============================================================================

@arts_agenda
def dobatch_calc_agenda_allsky(ws):

    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()

    # set cloudbox to full atmosphere
    ws.cloudboxSetFullAtm()

    # No jacobian calculations
    ws.jacobianOff()

    # no sensor
    ws.sensorOff()

    # set surface altitude
    ws.Extract(ws.z_surface, ws.array_of_z_surface, ws.ybatch_index)

    # set surface skin temperature
    ws.Extract(ws.DummyVariable, ws.vector_of_T_surface, ws.ybatch_index)
    ws.Tensor3SetConstant(ws.surface_props_data, 1, 1, 1, ws.DummyVariable)
    ws.Copy(ws.surface_skin_t, ws.DummyVariable)

    # set surface reflectivity
    ws.Extract(
        ws.surface_scalar_reflectivity,
        ws.array_of_surface_scalar_reflectivity,
        ws.ybatch_index,
    )

    # set geographical position
    ws.VectorExtractFromMatrix(ws.lat_true, ws.matrix_of_Lat, ws.ybatch_index, "row")
    ws.VectorExtractFromMatrix(ws.lon_true, ws.matrix_of_Lon, ws.ybatch_index, "row")

    # set sun position
    ws.Extract(ws.sun_pos, ws.array_of_sun_positions, ws.ybatch_index)
    ws.Extract(ws.sun_dist, ws.sun_pos, 0)
    ws.Extract(ws.sun_lat, ws.sun_pos, 1)
    ws.Extract(ws.sun_lon, ws.sun_pos, 2)

    ws.sunsChangeGeometry(
        index=ws.sun_index,
        distance=ws.sun_dist,
        latitude=ws.sun_lat,
        longitude=ws.sun_lon,
    )

    # setoutput for console
    ws.StringSet(ws.Text, "Allsky: DObatch Index")
    ws.Print(ws.Text, 0)
    ws.Print(ws.ybatch_index, 0)
    ws.Extract(ws.suns_do, ws.ArrayOfSuns_Do, ws.ybatch_index)
    ws.StringSet(ws.Text, "suns_do")
    ws.Print(ws.Text, 0)
    ws.Print(ws.suns_do, 0)

    # set cloudbox to full atm and calculate pndfield
    ws.cloudboxSetFullAtm()
    ws.pnd_fieldCalcFromParticleBulkProps()

    # do checks
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.scat_data_checkedCalc()

    # the actual rt simulation
    ws.spectral_irradiance_fieldDisort(
        nstreams=ws.NstreamIndex, Npfct=-1, emission=ws.EmissionIndex
    )

    # calculate fluxes
    ws.RadiationFieldSpectralIntegrate(
        ws.irradiance_field,
        ws.f_grid,
        ws.spectral_irradiance_field,
        quadrature_weights=ws.quadrature_weights,
    )

    # reset not needed quantities to save memory
    ws.Tensor7SetConstant(ws.spectral_radiance_field, 0, 0, 0, 0, 0, 0, 0, 0.0)
    ws.Tensor5SetConstant(ws.spectral_irradiance_field, 0, 0, 0, 0, 0, 0.0)
    ws.Tensor5SetConstant(ws.radiance_field, 0, 0, 0, 0, 0, 0.0)


@arts_agenda
def dobatch_calc_agenda_clearsky(ws):

    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()

    # set cloudbox to full atmosphere
    ws.cloudboxSetFullAtm()

    # No jacobian calculations
    ws.jacobianOff()

    # no sensor
    ws.sensorOff()

    # set surface altitude
    ws.Extract(ws.z_surface, ws.array_of_z_surface, ws.ybatch_index)

    # set surface skin temperature
    ws.Extract(ws.DummyVariable, ws.vector_of_T_surface, ws.ybatch_index)
    ws.Tensor3SetConstant(ws.surface_props_data, 1, 1, 1, ws.DummyVariable)
    ws.Copy(ws.surface_skin_t, ws.DummyVariable)

    # set surface reflectivity
    ws.Extract(
        ws.surface_scalar_reflectivity,
        ws.array_of_surface_scalar_reflectivity,
        ws.ybatch_index,
    )

    # set geographical position
    ws.VectorExtractFromMatrix(ws.lat_true, ws.matrix_of_Lat, ws.ybatch_index, "row")
    ws.VectorExtractFromMatrix(ws.lon_true, ws.matrix_of_Lon, ws.ybatch_index, "row")

    # set sun position
    ws.Extract(ws.sun_pos, ws.array_of_sun_positions, ws.ybatch_index)
    ws.Extract(ws.sun_dist, ws.sun_pos, 0)
    ws.Extract(ws.sun_lat, ws.sun_pos, 1)
    ws.Extract(ws.sun_lon, ws.sun_pos, 2)

    ws.sunsChangeGeometry(
        index=ws.sun_index,
        distance=ws.sun_dist,
        latitude=ws.sun_lat,
        longitude=ws.sun_lon,
    )

    ws.StringSet(ws.Text, "Clearsky: DObatch Index")
    ws.Print(ws.Text, 0)
    ws.Print(ws.ybatch_index, 0)
    ws.Extract(ws.suns_do, ws.ArrayOfSuns_Do, ws.ybatch_index)
    ws.StringSet(ws.Text, "suns_do")
    ws.Print(ws.Text, 0)
    ws.Print(ws.suns_do, 0)

    ws.cloudboxSetFullAtm()
    ws.pnd_fieldZero()

    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.scat_data_checkedCalc()

    ws.spectral_irradiance_fieldDisort(
        nstreams=ws.NstreamIndex, Npfct=-1, emission=ws.EmissionIndex
    )

    ws.RadiationFieldSpectralIntegrate(
        ws.irradiance_field,
        ws.f_grid,
        ws.spectral_irradiance_field,
        quadrature_weights=ws.quadrature_weights,
    )

    # reset not needed quantities to save memory
    ws.Tensor7SetConstant(ws.spectral_radiance_field, 0, 0, 0, 0, 0, 0, 0, 0.0)
    ws.Tensor5SetConstant(ws.spectral_irradiance_field, 0, 0, 0, 0, 0, 0.0)
    ws.Tensor5SetConstant(ws.radiance_field, 0, 0, 0, 0, 0, 0.0)


# =============================================================================
# aux functions
# =============================================================================


def create_agendas_in_ws(ws, pnd_agendas=False):
    ws, gas_scattering_agenda_list = create_gas_scattering_agendas_in_WS(ws)
    ws, surface_agenda_list = create_surface_agendas_in_WS(ws)

    pnd_agenda_list = []
    if pnd_agendas:
        ws, pnd_agenda_SB06_list = create_pnd_agendas_SB06_in_WS(ws)
        ws, pnd_agenda_MY05_list = create_pnd_agendas_MY05_in_WS(ws)
        ws, pnd_agenda_CG_list = create_pnd_agendas_CG_in_WS(ws)

        pnd_agenda_list = (
            pnd_agenda_SB06_list + pnd_agenda_MY05_list + pnd_agenda_CG_list
        )

    return ws, gas_scattering_agenda_list, surface_agenda_list, pnd_agenda_list
