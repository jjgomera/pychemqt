#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R21(MEoS):
    """Multiparameter equation of state for R21"""
    name = "dichlorofluoromethane"
    CASNumber = "75-43-4"
    formula = "CHCl2F"
    synonym = "R21"
    rhoc = unidades.Density(526.0138)
    Tc = unidades.Temperature(451.48)
    Pc = unidades.Pressure(5181.2, "kPa")
    M = 102.9227  # g/mol
    Tt = unidades.Temperature(142.8)
    Tb = unidades.Temperature(282.01)
    f_acent = 0.2061
    momentoDipolar = unidades.DipoleMoment(1.37, "Debye")
    id = 642

    CP1 = {"ao": 0.2376576/8.31451*102.92,
           "an": [0.12714330e-2/8.31451*102.92, 0.32413520e-6/8.31451*102.92,
                  -0.24924280e-8/8.31451*102.92, 0.17172080e-11/8.31451*102.92],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-21 of Platzer et al. (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., and Maurer, G.",
                    "title": "Thermophysical properties of refrigerants", 
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""}, 
        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP", 
        
        "Tmin": 200.0, "Tmax": 473.19, "Pmax": 137900.0, "rhomax": 15.36, 
        "Pmin": 0.6828e-4, "rhomin": 16.519, 

        "nr1": [-.44386484873e2, .926505600935e1, -.551709104376, .504676623431,
                -.732431415692, -.868403860387, .146234705555, -.280576335053,
                0.864743656093, -0.270767233732e1, 0.330476390706e1,
                -0.210878239171, 0.449531449589, 0.120779813143,
                -0.277297953777, 0.305441291172e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.443864848730e2, -0.926505600935e1, 0.551709104376,
                0.121128809552e1, 0.167119476587, -0.504876793028e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.07470252]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.0673], "exp": [1.23]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70336e1, 0.15672e1, -0.33932e1, 0.17582e1, -0.86765e1],
        "exp": [1.0, 1.5, 3.0, 7.0, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.33546, 0.18208e2, -0.26400e2, 0.10586e2],
        "exp": [0.09, 0.78, 0.92, 1.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.38213, -0.55559e1, -0.15886e2, -0.44766e2, -0.27606e3],
        "exp": [0.09, 0.667, 2.5, 6.0, 15.0]}
