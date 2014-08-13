#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R114(MEoS):
    """Multiparameter equation of state for R114

    >>> r114=R114(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r114.T, r114.rho, r114.u.kJkg, r114.h.kJkg, r114.s.kJkgK, r114.cv.kJkgK, r114.cp.kJkgK, r114.w)
    300.0 1558.7841 -33.86 -33.79 0.1152 0.67117 0.91951 693.14
    """
    name = "1,2-dichloro-1,1,2,2-tetrafluoroethane"
    CASNumber = "76-14-2"
    formula = "CClF2CClF2"
    synonym = "R114"
    rhoc = unidades.Density(579.969)
    Tc = unidades.Temperature(418.83)
    Pc = unidades.Pressure(3257.0, "kPa")
    M = 170.921  # g/mol
    Tt = unidades.Temperature(180.63)
    Tb = unidades.Temperature(276.741)
    f_acent = 0.25253
    momentoDipolar = unidades.DipoleMoment(0.658, "Debye")
    id = 231

    CP1 = {"ao": 0.97651380e-1,
           "an": [0.32408610e-2, -0.58953640e-5, 0.67379290e-8, -0.35463640e-11],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-114 of Platzer et al. (1990).",
        "__doc__":  u"""Platzer, B., Polt, A., and Maurer, G., "Thermophysical properties of refrigerants," Berlin:  Springer-Verlag, 1990.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1": [-0.340776521414, 0.323001398420, -0.424950537596e-1,
                0.107938879710e1, -0.199243619673e1, -0.155135133506,
                -0.121465790553, -0.165038582393e-1, -0.186915808643,
                0.308074612567, 0.115861416115, 0.276358316589e-1, 0.108043243088,
                0.460683793064e-1, -0.174821616881, 0.317530854287e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.340776521414, -0.323001398420, 0.424950537596e-1,
                -0.166940100976e1, 0.408693082002e1, -0.241738963889e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.21103865]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.05084], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.72195e1, 0.16357e1, -0.14576e1, -0.69580e1, 0.57181e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.43023, 0.22722e2, -0.27118e2, 0.13247e2, -0.90529e1],
        "exp": [0.095, 0.93, 1.1, 2.0, 3.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.46609, -0.68355e1, -0.16715e3, 0.15805e5, -0.31859e5, 0.21548e5],
        "exp": [0.09, 0.76, 4.0, 6.5, 7.0, 8.0]}
