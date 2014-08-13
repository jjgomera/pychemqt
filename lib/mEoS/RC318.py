#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class RC318(MEoS):
    """Multiparameter equation of state for RC318

    >>> rc318=RC318(T=500, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (rc318.T, rc318.rho, rc318.u.kJkg, rc318.h.kJkg, rc318.s.kJkgK, rc318.cv.kJkgK, rc318.cp.kJkgK, r365mfc.w)
    500.0 3.5926 288.51 316.34 1.0399 1.26519 1.32359 169.91
    """
    name = "octafluorocyclobutane "
    CASNumber = "406-58-6"
    formula = "cyclo-C4F8"
    synonym = "RC318"
    rhoc = unidades.Density(619.973)
    Tc = unidades.Temperature(388.38)
    Pc = unidades.Pressure(2777.5, "kPa")
    M = 200.0312  # g/mol
    Tt = unidades.Temperature(233.35)
    Tb = unidades.Temperature(267.175)
    f_acent = 0.3553
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 692

    CP1 = {"ao": 0.121,
           "an": [0.2903e-2, -0.25327e-5, 0.77191e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-C318 of Platzer et al. (1990)",
        "__doc__":  u"""Platzer, B., Polt, A., and Maurer, G., "Thermophysical properties of refrigerants," Berlin:  Springer-Verlag, 1990.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1": [-0.104729119796e1, 0.138034128822e1, -0.333625769594,
                0.109415569278e1, -0.268265237178e1, 0.173403063905e1,
                -0.163611246876e1, 0.304834499143, 0.102771551787,
                -0.232367895587e-1, 0.166151957803, -0.250103914479e-1,
                0.935680977639e-1, 0.431929127445e-1, -0.133439845861,
                0.255416632165e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.104729119796e1, -0.138034128822e1, 0.333625769594,
                -0.510485781618, 0.181840728111e1, -0.138530893970e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.99943992]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.05145], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.78467e1, 0.24555e1, -0.30824e1, -0.58263e1, 0.35483e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.30181, 0.29345e1, -0.13741e1, 0.14650e1, 0.16963],
        "exp": [0.11, 0.32, 0.57, 0.84, 2.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.24491e2, 0.53255e2, -0.38863e2, -0.24938e2, -0.90092e2],
        "exp": [0.61, 0.77, 0.92, 3.3, 7.5]}
