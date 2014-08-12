#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Acetone(MEoS):
    """Multiparameter equation  of state for Acetone

    >>> acetona=Acetone(T=300, P=0.1)
    >>> print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.4f %0.1f" % (acetona.T, acetona.rho, acetona.u.kJkg, acetona.h.kJkg, acetona.s.kJkgK, acetona.cv.kJkgK, acetona.cp.kJkgK, acetona.w)
    300.0 782.6222 256.774 256.902 1.19826 1.5527 2.1475 1155.7
    """
    name = "acetone"
    CASNumber = "67-64-1"
    formula = "CH3COCH3"
    synonym = ""
    rhoc = unidades.Density(272.97)
    Tc = unidades.Temperature(508.1)
    Pc = unidades.Pressure(4700.0, "kPa")
    M = 58.07914  # g/mol
    Tt = unidades.Temperature(178.5)
    Tb = unidades.Temperature(329.22)
    f_acent = 0.3071
    momentoDipolar = unidades.DipoleMoment(2.88, "Debye")
    id = 140

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [3.7072, 7.0675, 11.012],
           "exp": [310, 3480, 1576],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for acetone of Lemmon and Span (2006)",
        "__doc__":  u"""Lemmon, E.W., Span, R. Short Fundamental Equations of State for 20 Industrial Fluids. J. Chem. Eng. Data 51 (2006), 785 – 850.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1":  [0.90041, -2.1267, -0.083409, 0.065683, 0.00016527],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.039663, 0.72085, 0.0092318, -0.17217, -0.14961, -0.076124,
                -0.018166],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eqi = helmholtz1,

    _surface = {"sigma": [0.07], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.76214e1, 0.17441e1, -0.20514e1, -0.26644e1, -0.69437],
        "exp": [1, 1.5, 2.57, 4.43, 15.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.11118e2, -0.29507e2, 0.35255e2, -0.14712e2, 0.95560],
        "exp": [0.456, 0.626, 0.8, 1.0, 2.47]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.25200e1, -0.66065e1, -0.25751e2, 0.78120e1, -0.53778e2, -0.11684e3],
        "exp": [0.36, 1.05, 3.2, 4.0, 6.5, 14.0]}
