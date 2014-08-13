#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R407C(MEoS):
    """Multiparamenter equation of state for R407C
    23% R32, 25% R125, 52% R134a

    >>> r407C=R407C(T=500, P=0.1)
    >>> print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (r407C.T, r407C.rho, r407C.h.kJkg, r407C.s.kJkgK, r407C.cv.kJkgK, r407C.cp.kJkgK, r407C.w)
    500.0 2.0783 195.613 0.479 1.03122 1.1291 229.27
    """
    name = "R407C"
    CASNumber = ""
    formula = "R32+R125+R134a"
    synonym = "R407C"
    rhoc = unidades.Density(453.430936)
    Tc = unidades.Temperature(359.345)
    Pc = unidades.Pressure(4631.7, "kPa")
    M = 86.2036  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(229.52)
    f_acent = 0.363
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 62

    CP1 = {"ao": 0.,
           "an": [0.76575], "pow": [0.4],
           "ao_exp": [1.4245, 3.9419, 3.1209], "exp": [864, 1887, 4802],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-407C of Lemmon (2003)",
        "__doc__":  u"""Lemmon, E.W., "Pseudo Pure-Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C," Int. J. Thermophys., 24(4):991-1006, 2003.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [0.105880e1, -0.112018e1, 0.629064, -0.351953, 0.455978e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.241, 0.69, 2.58, 1.15, 0.248],

        "nr2": [-0.175725e1, -0.112009e1, 0.277353e-1, 0.898881, -0.117591e1,
                0.818591e-1, -0.794097e-1, -0.104047e-4, 0.233779, -0.291790,
                0.154776e-1, -0.314579e-1, -0.442552e-2, -0.101254e-1,
                0.915953e-2, -0.361575e-2],
        "d2": [1, 2, 2, 3, 3, 5, 5, 5, 1, 1, 4, 4, 2, 4, 5, 6],
        "t2": [2.15, 2.43, 5.3, 0.76, 1.48, 0.24, 2.86, 8., 3.3, 4.7, 0.45,
               8.4, 16.2, 26, 16, 8.7],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
        "gamma2": [1]*16}

    eq = helmholtz1,

    _surface = {"sigma": [0.064017], "exp": [1.2557]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.086077, -6.6364, -2.4648, -3.4776],
        "exp": [0.4, 0.965, 3.1, 5.]}
    _liquid_Pressure = {
        "eq": 5,
        "ao": [0.48722, -6.6959, -1.4165, -2.5109],
        "exp": [0.54, 0.925, 2.7, 4.7]}
