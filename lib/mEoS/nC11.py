#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class nC11(MEoS):
    """Multiparameter equation of state for n-undecane"""
    name = "undecane"
    CASNumber = "1120-21-4"
    formula = "CH3-9(CH2)-CH3"
    synonym = ""
    rhoc = unidades.Density(236.791383074)
    Tc = unidades.Temperature(638.8)
    Pc = unidades.Pressure(1990.4, "kPa")
    M = 156.30826  # g/mol
    Tt = unidades.Temperature(247.541)
    Tb = unidades.Temperature(468.934)
    f_acent = 0.539
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 15

    Fi1 = {"ao_log": [1, -120.4274],
           "pow": [-3, -2, -1, 0, 1, 2],
           "ao_pow": [-3.515339, 28.27708, -136.8378, -46.40384, 107.1876, 
                      1.419929],
           "tau*logtau": -31.81246, 
           "tau*logdeñta": 0, 
           "ao_exp": [], "titao": [], 
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for undecane of Aleksandrov et al. (2011)",
        "__doi__": {"autor": "Aleksandrov, I.S., Gerasimov, A.A., and Grigor’ev, B.A.",
                    "title": "Using fundamental equations of state for calculating the thermodynamic properties of normal undecane", 
                    "ref": "Thermal Engineering, 58(8):691-698, 2011",
                    "doi": "10.1134/S0040601511080027"}, 

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 700., "Pmax": 500000.0, "rhomax": 4.97, 
        "Pmin": 0.0004461, "rhomin": 4.962, 

        "nr1": [-0.66172706, 1.3375396, -2.5608399, 0.1067891, 0.28873614e-3,
                0.49587209e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [0.55407101e-7, 0.99754712, 1.5774025, 0.13108354e-2,
                -0.59326961, -0.93001876e-1, -0.17960228, -0.22560853e-1],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-9.3961, 4.4531, -5.2658, -4.7352],
        "exp": [1, 1.5, 2.2, 4.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [4.5273, -7.5714, 13.920, -13.464, 5.8411],
        "exp": [0.46, 0.84, 1.25, 1.7, 2.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-4.3093, -3.4358, -17.473, -58.573, -133.83],
        "exp": [0.466, 1.02, 2.4, 5.3, 11.4]}
