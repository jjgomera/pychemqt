#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class DME(MEoS):
    """Multiparameter equation of state for dimethylether"""
    name = "dimethylether"
    CASNumber = "115-10-6"
    formula = "CH3-O-CH3"
    synonym = "R-170"
    rhoc = unidades.Density(273.6465336)
    Tc = unidades.Temperature(400.378)
    Pc = unidades.Pressure(5336.8, "kPa")
    M = 46.06844  # g/mol
    Tt = unidades.Temperature(131.66)
    Tb = unidades.Temperature(248.368)
    f_acent = 0.196
    momentoDipolar = unidades.DipoleMoment(1.301, "Debye")
    id = 133

    CP1 = {"ao": 4.039,
           "an": [], "pow": [],
           "ao_exp": [2.641, 2.123, 8.992, 6.191],
           "exp": [361, 974, 1916, 4150],
           "ao_hyp": [], "hyp": []}

    Fi1 = {"ao_log": [1, 3.039],
           "pow": [0, 1],
           "ao_pow": [-1.928925, 3.150284],
           "ao_exp": [2.641, 2.123, 8.992, 6.191],
           "titao": [361/Tc, 974/Tc, 1916/Tc, 4150/Tc], 
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for DME of Wu et al. (2011).",
        "__doi__": {"autor": "Wu, J., Zhou, Y., and Lemmon, E.W.",
                    "title": "An Equation of State for the Thermodynamic Properties of Dimethyl Ether", 
                    "ref": "J. Phys. Chem. Ref. Data 40, 023104 (2011)",
                    "doi":  "10.1063/1.3582533"}, 
        "R": 8.314472,
        "cp": CP1,
        "ref": {"Tref": 273.15, "Pref": 1.0, "ho": 23242, "so": 131.3883}, 

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 40000.0, "rhomax": 19.15, 
        "Pmin": 0.0022, "rhomin": 19.15, 

        "nr1":  [0.29814139e-1, 0.14351700e1, -0.26496400e1, -0.29515532,
                 0.17035607],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.4366, 1.011, 1.137, 0.45],

        "nr2": [-0.94642918, -0.99250514e-1, 0.11264071e1, -0.76936548,
                -0.20717696e-1, 0.24527037],
        "d2": [1, 3, 2, 2, 7, 1],
        "t2": [2.83, 1.5, 1.235, 2.675, 0.7272, 1.816],
        "c2": [2, 2, 1, 2, 1, 1],
        "gamma2": [1]*6,

        "nr3": [0.11863438e1, -0.49398368, -0.16388716, -0.27583584e-1],
        "d3": [1, 1, 3, 3],
        "t3": [1.783, 3.779, 3.282, 1.059],
        "alfa3": [0.965336, 1.508580, 0.963855, 9.726430],
        "beta3": [1.287190, 0.806235, 0.777942, 197.681000],
        "gamma3": [1.277720, 0.430750, 0.429607, 1.138490],
        "epsilon3": [0.672698, 0.924246, 0.750815, 0.800022],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for DME of Ihmels and Lemmon (2007)",
        "__doi__": {"autor": "Ihmels, E.C. and Lemmon, E.W.",
                    "title": "Experimental densities, vapor pressures, and critical point, and a fundamental equation of state for dimethyl ether", 
                    "ref": " Fluid Phase Equilibria (Impact Factor: 2.24). 10/2007; 260(1):36-48",
                    "doi":  "10.1016/j.fluid.2006.09.016"}, 
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 1350.0, "Pmax": 1000000.0, "rhomax": 73.96, 
        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1":  [1.22690, -2.47245, 0.119889, 0.0000354],
        "d1": [1, 1, 3, 8],
        "t1": [0.21, 1.0, 0.5, 1.0],

        "nr2": [0.567139, 0.166649, -0.078412, -0.289066, -0.031272, -0.065607],
        "d2": [2, 1, 5, 1, 4, 3],
        "t2": [1.4, 3.1, 1.5, 5.0, 5.9, 3.7],
        "c2": [1, 1, 1, 2, 2, 2],
        "gamma2": [1]*6}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.063157], "exp": [1.2595]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70647e1, 0.17709e1, -0.21544e1, -0.22098e1],
        "exp": [1, 1.5, 2.6, 5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.13362e-1, 0.75723e1, -0.10108e2, 0.52885e1, 0.53918e-1],
        "exp": [0.188, 0.53, 0.73, 0.93, 4.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.408630565713e1, -0.441287217274e1, -0.122912228629e2,
               -0.393407197305e2, -0.159976585292e1, -0.864097952455e2],
        "exp": [0.486, 1.39, 2.7, 5.7, 10., 12.]}
