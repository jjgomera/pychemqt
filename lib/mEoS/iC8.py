#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class iC8(MEoS):
    """Multiparameter equation of state for isooctane"""
    name = "isooctane"
    CASNumber = "540-84-1"
    formula = "(CH3)2CHCH2C(CH3)3"
    synonym = ""
    rhoc = unidades.Density(242.1644624)
    Tc = unidades.Temperature(544)
    Pc = unidades.Pressure(2572.0, "kPa")
    M = 114.22852  # g/mol
    Tt = unidades.Temperature(165.77)
    Tb = unidades.Temperature(372.358)
    f_acent = 0.303
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 82

    CP1 = {"ao": 10.76,
           "an": [], "pow": [],
           "ao_exp": [15.48, 34.42, 21.42],
           "exp": [775, 1900, 5100],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isooctane of Blackham and Lemmon (2011).",
        "__doi__": {"autor": "Blackham, T.M. and Lemmon, E.W.",
                    "title": "", 
                    "ref": "to be published in Int. J. Thermophys., 2011.",
                    "doi": ""}, 
            
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 1000000.0, "rhomax": 6.97, 
        "Pmin": 0.00001796, "rhomin": 6.96, 

        "nr1": [0.568901e-1, 0.196155e1, -0.281164e1, -0.815112, 0.326583],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 0.75, 1.11, 0.55],

        "nr2": [-0.160893e1, -0.454734, 0.108306e1, -0.722876, -0.434052e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.2, 3.7, 1.53, 2.1, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*6, 
  
        "nr3": [0.196648e1, -0.465082, -0.409398, 0.232131e-1],
        "d3": [1, 1, 3, 3, ],
        "t3": [0.88, 1.1, 2.75, 1.0],
        "alfa3": [0.75, 1.13, 0.87, 4.73],
        "beta3": [0.59, 1.45, 0.5, 10.52],
        "gamma3": [1.44, 0.68, 0.51, 0.8],
        "epsilon3": [0.66, 0.9, 0.54, 0.18],
        "nr4": []}

    eq = helmholtz1,
    _PR = -0.0058658

    _surface = {"sigma": [0.0476182, 5.992036e-17], "exp": [1.1914, 2.1914]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.7985, 8.1280, -7.3106, -3.9392, -1.6732],
        "exp": [1, 1.5, 1.6, 4.0, 16.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [1.1535, 1.3709, 0.38804],
        "exp": [0.286, 0.54, 3.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-2.5793, -6.4934, -18.631, -54.123, -123.58],
        "exp": [0.366, 1.11, 3.0, 6.4, 14.0]}

    thermo0 = {"eq": 5, "omega": 3,
               "__name__": "Chung (1988)",
               "__doi__": {"autor": "T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E.",
                           "title": "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties", 
                           "ref": "Ind. Eng. Chem. Res., 1988, 27 (4), pp 671–679",
                           "doi": "10.1021/ie00076a024"}, 
               "w": 0.3, "mur": 0.0, "k": 0.0}

    _viscosity = thermo0,
    _thermal = thermo0, 
