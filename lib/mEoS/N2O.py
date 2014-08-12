#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class N2O(MEoS):
    """Multiparameter equation of state for nitrous oxide

    >>> n2o=N2O(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (n2o.T, n2o.rho, n2o.u.kJkg, n2o.h.kJkg, n2o.s.kJkgK, n2o.cv.kJkgK, n2o.cp.kJkgK, n2o.w)
    300.0 1.77389 413.83 470.200 2.43074 0.6915 0.8848 267.85
    """
    name = "nitrous oxide"
    CASNumber = "10024-97-2"
    formula = "N2O"
    synonym = "R-744A"
    rhoc = unidades.Density(452.011)
    Tc = unidades.Temperature(309.52)
    Pc = unidades.Pressure(7245.0, "kPa")
    M = 44.0128  # g/mol
    Tt = unidades.Temperature(182.33)
    Tb = unidades.Temperature(184.68)
    f_acent = 0.1613
    momentoDipolar = unidades.DipoleMoment(0.1608, "Debye")
    id = 110

    CP1 = {"ao": 3.5,
           "an": [], "pow": [],
           "ao_exp": [2.1769, 1.6145, 0.48393],
           "exp": [879, 2372, 5447],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for nitrous oxide of Lemmon and Span (2006)",
        "__doc__":  u"""Lemmon, E.W., Span, R. Short fundamental equations of state for 20 industrial fluids. J. Chem. Eng. Data 51 (2006), 785 – 850.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1":  [0.88045, -2.4235, 0.38237, 0.068917, 0.00020367],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.13122, 0.46032, -0.0036985, -0.23263, -0.00042859, -0.042810,
                -0.023038],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [0.0745], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.69078e1, 0.26620e1, -0.22386e1, -0.38002e1, 0.76922],
        "exp": [1.0, 1.5, 1.9, 4.8, 5.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.67919e1, -0.16069e2, 0.25632e2, -0.20755e2, 0.71963e1],
        "exp": [0.47, 0.72, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31287e1, -0.77651e2, 0.21442e3, -0.47809e3, 0.75185e3, -0.46279e3],
        "exp": [0.409, 1.91, 2.33, 3.0, 3.6, 4.0]}
