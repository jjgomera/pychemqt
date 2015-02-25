#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R404A(MEoS):
    """Multiparameter equation of state for R404A
    44% R125, 4% R134a, 52% R143a

    >>> r404A=R404A(T=500, P=0.1)
    >>> print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (r404A.T, r404A.rho, r404A.h.kJkg, r404A.s.kJkgK, r404A.cv.kJkgK, r404A.cp.kJkgK, r404A.w)
    500.0 2.3535 205.479 0.481 1.08723 1.1738 213.92
    """
    name = "R404A"
    CASNumber = ""
    formula = "R125+R134a+R143a"
    synonym = "R404A"
    rhoc = unidades.Density(482.162772)
    Tc = unidades.Temperature(345.27)
    Pc = unidades.Pressure(3734.8, "kPa")
    M = 97.6038  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(226.93)
    f_acent = 0.293
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 62

    CP1 = {"ao": 0.,
           "an": [1.2744], "pow": [0.3],
           "ao_exp": [0.63075, 3.5979, 5.0335], "exp": [413, 804, 1727],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-404A of Lemmon (2003)",
        "__doc__":  u"""Lemmon, E.W., "Pseudo Pure-Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C," Int. J. Thermophys., 24(4):991-1006, 2003.""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 14.21, 
        "Pmin": 22.65, "rhomin": 14.21, 

        "nr1": [0.610984e1, -0.779453e1, 0.183377e-1, 0.262270, -0.351688e-2,
                0.116181e-1, 0.105992e-2],
        "d1": [1, 1, 1, 2, 2, 4, 6],
        "t1": [0.67, 0.91, 5.96, 0.7, 6, 0.3, 0.7],

        "nr2": [0.850922, -0.520084, -0.464225e-1, 0.62119, -.195505, .336159,
                -.376062e-1, -.636579e-2, -.758262e-1, -.221041e-1, .310441e-1,
                0.132798e-1, 0.689437e-1, -0.507525e-1, 0.161382e-1],
        "d2": [1, 1, 1, 2, 2, 3, 4, 7, 2, 3, 4, 4, 2, 3, 5],
        "t2": [1.7, 3.3, 7, 2.05, 4.3, 2.7, 1.8, 1.25, 12, 6, 8.7, 11.6, 13,
               17, 16],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*15}

    eq = helmholtz1,

    _surface = {"sigma": [0.06868, -0.04576], "exp": [1.3142, 2.3084]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.00026863, -6.5757, -4.1802, -7.9102],
        "exp": [0.1, 0.972, 3.8, 9]}
    _liquid_Pressure = {
        "eq": 5,
        "ao": [0.061067, -6.5646, -3.6162, -3.9771],
        "exp": [0.54, 0.965, 3.7, 9.0]}
