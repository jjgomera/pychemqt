#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class C1Linolenate(MEoS):
    """Multiparameter equation of state for methyl linolenate"""
    name = "methyl linolenate"
    CASNumber = "301-00-8"
    formula = "C19H32O2"
    synonym = ""
    rhoc = unidades.Density(247.798121314)
    Tc = unidades.Temperature(772.0)
    Pc = unidades.Pressure(1369.0, "kPa")
    M = 292.45618  # g/mol
    Tt = unidades.Temperature(218.65)
    Tb = unidades.Temperature(629.13)
    f_acent = 1.14
    momentoDipolar = unidades.DipoleMoment(1.54, "Debye")
    id = 39

    CP1 = {"ao": 0.0,
           "an": [79.5913], "pow": [0.214648],
           "ao_exp": [290.379, 81.4323, 474.881],
           "exp": [1213.24, 578.752, 2799.79],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl linoleate of Huber et al. (2009).",
        "__doi__": {"autor": "Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, L.S., and Bruno, T.J.",
                    "title": "Model for the Thermodynamic Properties of a Biodiesel Fuel", 
                    "ref": "Energy Fuels, 2009, 23 (7), pp 3790–3797",
                    "doi": "10.1021/ef900159g"}, 
            
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 3.29, 
        "Pmin": 1.e-17, "rhomin": 3.28, 

        "nr1": [0.4070829e-1, 0.2412375e1, -0.3756194e1, -0.1526466, 0.4682918e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.15, 1.24, 1.6, 1.28],

        "nr2": [-0.1470958e1, -0.76455, 0.1908964e1, -0.1629366e1, -0.1242073e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.9, 3.15, 2.16, 2.8, 1.4],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2180707e1, -0.7537264, -0.4347781],
        "d3": [1, 1, 3],
        "t3": [2.5, 3.0, 3.1],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.14239e2, 0.85361e1, -0.29678e2, 0.29106e2, -0.16922e2],
        "exp": [1.0, 1.5, 2.86, 3.5, 4.6]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.16258e2, -0.19465e3, 0.14231e4, -0.22420e4, 0.10012e4],
        "exp": [0.681, 1.26, 1.58, 1.7, 1.8]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.11463e2, 0.45192e2, -0.65779e2, -0.18386e4, 0.40689e4, -0.25124e4],
        "exp": [0.65, 1.55, 1.8, 6.6, 7.2, 7.8]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2010)",
               "__doc__": """Perkins, R.A. and Huber, M.L., unpublished work, 2010.""",

               "Tref": 772.0, "kref": 1,
               "no": [-0.27125000e-3, 0.25936500e-2, 0.35024100e-1, -0.90227300e-2],
               "co": [0, 1, 2, 3],

               "Trefb": 772.0, "rhorefb": 0.8473*M, "krefb": 1.,
               "nb": [-0.41010600e-1, 0.32844300e-1, -0.41850600e-2, 0.0, 0.0,
                      0.60665700e-1, -0.49840700e-1, 0.12175200e-1, 0.0, 0.0],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 8.75e-10, "Tcref": 1158.0}

    _thermal = thermo0,
