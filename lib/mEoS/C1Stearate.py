#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class C1Stearate(MEoS):
    """Multiparameter equation of state for methyl stearate"""
    name = "methyl stearate"
    CASNumber = "112-61-8"
    formula = "C19H38O2"
    synonym = ""
    rhoc = unidades.Density(237.101584226)
    Tc = unidades.Temperature(775.0)
    Pc = unidades.Pressure(1239.0, "kPa")
    M = 298.50382  # g/mol
    Tt = unidades.Temperature(311.84)
    Tb = unidades.Temperature(629.56)
    f_acent = 1.02
    momentoDipolar = unidades.DipoleMoment(1.54, "Debye")
    id = 39

    CP1 = {"ao": 0.0,
           "an": [247.115], "pow": [-0.0916606],
           "ao_exp": [276.94, 408.997, 472.702],
           "exp": [556.17, 1311.85, 2825.71],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl estearate of Huber et al. (2009).",
        "__doi__": {"autor": "Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, L.S., and Bruno, T.J.",
                    "title": "Model for the Thermodynamic Properties of a Biodiesel Fuel", 
                    "ref": "Energy Fuels, 2009, 23 (7), pp 3790–3797",
                    "doi": "10.1021/ef900159g"}, 
            
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 2.86, 
        "Pmin": 0.00000576, "rhomin": 2.85, 

        "nr1": [0.3959635e-1, 0.2466654e1, -0.3895950e1, -0.1167375, 0.4127229e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 1.25, 1.65, 0.8],

        "nr2": [-0.1403734e1, -0.6465264, 0.1934675e1, -0.1608124e1, -0.1113813e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [3.1, 3.4, 2.3, 3.8, 1.2],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2125325e1, -0.7772671, -0.4183684],
        "d3": [1, 1, 3],
        "t3": [3.2, 3.8, 3.8],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.14597e2, 0.13836e2, -0.14484e2, -0.51856e1, -0.27076e1],
        "exp": [1.0, 1.5, 2.12, 4.7, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.11202e2, 0.78636e2, -0.12554e3, 0.72942e2, -0.11524e2],
        "exp": [0.439, 0.59, 0.73, 0.9, 1.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.18187e2, 0.81619e2, -0.90210e2, -0.52888e3, 0.11270e4, -0.84453e3],
        "exp": [0.71, 1.3, 1.5, 6.0, 7.0, 8.0]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2010)",
               "__doc__": """Perkins, R.A. and Huber, M.L., unpublished work, 2010.""",

               "Tref": 775.0, "kref": 1,
               "no": [-0.27125e-3, 0.259365e-2, 0.350241e-1, -0.902273e-2],
               "co": [0, 1, 2, 3],

               "Trefb": 775.0, "rhorefb": 0.7943*M, "krefb": 1.,
               "nb": [-0.410106e-1, 0.328443e-1, -0.418506e-2, 0.0, 0.0,
                      0.606657e-1, -0.498407e-1, 0.121752e-1, 0.0, 0.0],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 8.75e-10, "Tcref": 1162.5}

    _thermal = thermo0,
