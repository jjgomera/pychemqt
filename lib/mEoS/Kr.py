#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Kr(MEoS):
    """Multiparameter equation of state for krypton

    >>> kripton=Kr(T=300, P=0.1)
    >>> print "%0.1f %0.3f %0.3f %0.3f %0.5f %0.4f %0.4f %0.2f" % (kripton.T, kripton.rho, kripton.u.kJkg, kripton.h.kJkg, kripton.s.kJkgK, kripton.cv.kJkgK, kripton.cp.kJkgK, kripton.w)
    300.0 3.366 123.073 152.779 1.12967 0.1490 0.2492 222.65
    """
    name = "krypton"
    CASNumber = "7439-90-9"
    formula = "Kr"
    synonym = ""
    rhoc = unidades.Density(909.2083)
    Tc = unidades.Temperature(209.48)
    Pc = unidades.Pressure(5525.0, "kPa")
    M = 83.798  # g/mol
    Tt = unidades.Temperature(115.775)
    Tb = unidades.Temperature(119.73)
    f_acent = -0.00089
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    # id = 971
    id = 1

    CP1 = {"ao": 2.5,
           "an": [], "pow": [], "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for krypton of Lemmon and Span (2006).",
        "__doc__": u"""Lemmon, E.W. and Span, R., "Short Fundamental Equations of State for 20 Industrial Fluids," J. Chem. Eng. Data, 51:785-850, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 200000.0, "rhomax": 33.42, 
        "Pmin": 73.5, "rhomin": 29.2, 

        "nr1": [0.83561, -2.3725, 0.54567, 0.014361, 0.066502, 0.00019310],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.16818, -0.033133, -0.15008, -0.022897, -0.021454, 0.0069397],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for krypton of Polt et al. (1992).",
        "__doc__": u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 780.0, "Pmax": 375000.0, "rhomax": 33.55, 
        "Pmin": 73.476, "rhomin": 29.249, 

        "nr1": [-0.402218741560, 0.679250544381, -0.1878869802860,
                0.603399982935, -0.177297564389e1, 0.581208430222,
                -0.733585469788, 0.164651929067, -0.319923148922e-1,
                0.333278228743, 0.219652478083e-1, 0.751994891628e-1,
                -0.212109737251, -0.645185506524e-2, 0.409175610200e-1,
                0.169416098754e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.402218741560, -0.679250544381, 0.187886980286,
                0.108265263587, -0.137102675805, -0.110549803007],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2, 2, 2, 2, 2, 2],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.0431], "exp": [1.2]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [6.273], "expt1": [0], "expd1": [1],
                   "a2": [6.485, 13.48, -82.51, -170.4],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.7, 2.7]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 101.325,
                "Tmin": Tt, "Tmax": 800.0,
                "a1": [-2345.757, 1.080476685], "exp1": [0, 1.6169841],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 73.197,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-11.5616], "exp2": [1],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.59697e1, 0.12673e1, -0.95609, -0.35630e2, 0.56884e2],
        "exp": [1.0, 1.5, 2.95, 9.3, 10.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.20593e2, -0.65490e2, 0.94407e2, -0.69678e2, 0.22810e2],
        "exp": [0.62, 0.84, 1.07, 1.34, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.64163e1, 0.89956e1, -0.10216e2, -0.13477e2, -0.21152e3, 0.21375e3],
        "exp": [0.525, 0.77, 1.04, 3.2, 8.3, 9.0]}
