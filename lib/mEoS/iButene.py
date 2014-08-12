#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class iButene(MEoS):
    """Multiparameter equation of state for isobutene

    >>> buteno=iButene(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (buteno.T, buteno.rho, buteno.u.kJkg, buteno.h.kJkg, buteno.s.kJkgK, buteno.cv.kJkgK, buteno.cp.kJkgK, buteno.w)
    300.0 2.31141 406.68 449.941 1.64447 1.4647 1.6344 216.67
    """
    name = "isobutene"
    CASNumber = "115-11-7"
    formula = "CH2=C(CH3)2"
    synonym = ""
    rhoc = unidades.Density(233.963)
    Tc = unidades.Temperature(418.09)
    Pc = unidades.Pressure(4009.8, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(132.4)
    Tb = unidades.Temperature(266.15)
    f_acent = 0.193
    momentoDipolar = unidades.DipoleMoment(0.5, "Debye")
    id = 27

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [4.8924, 7.832, 7.2867, 8.7293],
           "exp": [399, 1270, 2005, 4017],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isobutene of Lemmon and Ihmels (2005)",
        "__doc__":  u"""Lemmon, E.W., Ihmels, E.C. Thermodynamic properties of the butenes Part II. Short fundamental equations of state. Fluid Phase Equilibria 228 – 229 (2004), 173 – 187.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1":  [0.77111, -2.7971, 1.0118, 0.02073, 0.085086, 0.0021968],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.20633, -0.078843, -0.23726, -0.080211, -0.027001, 0.013072],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.0551], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.68973e1, 0.12475e1, -0.25441e1, -0.29282e1, 0.15778e1],
        "exp": [1., 1.5, 3.16, 6.2, 7.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.62591e2, -0.20805e3, 0.33243e3, -0.29555e3, 0.11148e3],
        "exp": [0.65, 0.8, 0.98, 1.16, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31841e1, -0.64014e1, -0.93817e1, -0.11160e2, -0.52298e2, -0.12195e3],
        "exp": [0.431, 1.29, 3.3, 3.54, 7.3, 15.8]}
