#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Butene_1(MEoS):
    """Multiparameter equation of state for 1-butene

    >>> buteno=Butene_1(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.5f %0.5f %0.2f" % (buteno.T, buteno.rho, buteno.h.kJkg, buteno.s.kJkgK, buteno.cv.kJkgK, buteno.cp.kJkgK, buteno.w)
    300.0 2.31151 -4.169 0.00742 1.39641 1.56581 217.19
    """
    name = "butene"
    CASNumber = "106-98-9"
    formula = "CH3-CH2-CH=CH2"
    synonym = ""
    rhoc = unidades.Density(237.89)
    Tc = unidades.Temperature(419.29)
    Pc = unidades.Pressure(4005.1, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(87.8)
    Tb = unidades.Temperature(266.84)
    f_acent = 0.192
    momentoDipolar = unidades.DipoleMoment(0.339, "Debye")
    id = 24

    CP1 = {"ao": 3.9197,
           "an": [], "pow": [],
           "ao_exp": [2.9406, 6.5395, 14.5395, 5.8971],
           "exp": [274, 951, 2127, 5752],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for 1-butene of Lemmon and Ihmels (2005)",
        "__doc__":  u"""Lemmon, E.W., Ihmels, E.C. Thermodynamic properties of the butenes Part II. Short fundamental equations of state. Fluid Phase Equilibria 228 – 229 (2004), 173 – 187.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 525., "Pmax": 70000.0, "rhomax": 14.59, 
        "Pmin": 0.0000000008, "rhomin": 14.58, 

        "nr1": [0.78084, -2.8258, 0.99403, 0.017951, 0.088889, 0.00024673],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.22846, -0.074009, -0.22913, -0.062334, -0.025385, 0.011040],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.0566], "exp": [1.25]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.71727e1, 0.26360e1, -0.20781e1, -0.28860e1, -0.13041e1],
        "exp": [1, 1.5, 2, 4.35, 16.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.16857e2, -0.46280e2, 0.53727e2, -0.23314e2, 0.18889e1],
        "exp": [0.547, 0.73, 0.92, 1.14, 2.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31106e1, -0.63103e1, -0.19272e2, -0.48739e2, -0.99898e2, -0.19001e3],
        "exp": [0.415, 1.27, 3.34, 7.0, 14.5, 28.0]}
