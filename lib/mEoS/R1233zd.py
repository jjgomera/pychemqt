#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R1233zd(MEoS):
    """Multiparameter equation of state for R1233zd"""
    name = "1-chloro-3,3,3-trifluoroprop-1-ene"
    CASNumber = "102687-65-0"
    formula = "CHCl=CH-CF3"
    synonym = "R1233zd"
    rhoc = unidades.Density(476.31109204)
    Tc = unidades.Temperature(438.75)
    Pc = unidades.Pressure(3570.9, "kPa")
    M = 130.4961896  # g/mol
    Tt = unidades.Temperature(195.15)
    Tb = unidades.Temperature(291.47)
    f_acent = 0.305
    momentoDipolar = unidades.DipoleMoment(1.44, "Debye")
#    id = 671

    CP1 = {"ao": 4.0,
           "an": [], "pow": [],
           "ao_exp": [8.962, 11.94],
           "exp": [400, 1900],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1233zd(E) of Mondejar et al. (2013).",
        "__doi__": {"autor": "Mondejar, M.E., McLinden, M.O., Lemmon, E.W.",
                    "title": "Thermodynamic Properties of Trans-1-chloro-3,3,3-trifluoropropene (R1233zd(E)): Vapor Pressure, p-rho-T Data, Speed of Sound Measurements and Equation of State", 
                    "ref": "to be submitted to J. Chem. Eng. Data, 2013.",
                    "doi": ""}, 
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 550.0, "Pmax": 100000.0, "rhomax": 11.41, 
        "Pmin": 0.25, "rhomin": 11.41, 

        "nr1": [0.03920261, 1.639052, -1.997147, -0.6603372, 0.1498682],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1., 0.24, 0.83, 1.17, 0.6],

        "nr2": [-1.408791, -0.7920426, 0.8549678, -0.5301920, -0.01408562],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.2, 2.88, 1.1, 2., 1.07],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,
 
        "nr3": [1.335117, -0.5441797, -0.05862723, -0.04123614, -0.6619106],
        "d3": [1, 1, 3, 2, 3],
        "t3": [1.27, 1.94, 2., 1.5, 1.],
        "alfa3": [1.215, 1.5, 1.1, 2.52, 4.55],
        "beta3": [1.27, 0.82, 0.94, 20., 32.],
        "gamma3": [1.32, 0.82, 0.66, 0.66, 1.39],
        "epsilon3": [0.77, 0.976, 1.08, 0.62, 0.61],
        "nr4": []}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.6021, 2.3265, -1.9771, -4.8451, -4.8762],
        "exp": [1.0, 1.5, 2.0, 4.3, 14.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [2.13083, 0.583568, 0.247871, 0.472173],
        "exp": [0.355, 0.9, 3.5, 8.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.0152, -6.5621, -19.427, -62.650, -181.64],
        "exp": [0.397, 1.2, 3.1, 6.6, 15.0]}
