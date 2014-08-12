#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class C1Palmitate(MEoS):
    """Multiparameter equation of state for methyl palmitate

    >>> metilciclohexano=C1Palmitate(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "methyl palmitate"
    CASNumber = "112-39-0"
    formula = "C17H34O2"
    synonym = ""
    rhoc = unidades.Density(242.59424202)
    Tc = unidades.Temperature(755.0)
    Pc = unidades.Pressure(1350.0, "kPa")
    M = 270.45066  # g/mol
    Tt = unidades.Temperature(302.71)
    Tb = unidades.Temperature(602.3)
    f_acent = 0.91
    momentoDipolar = unidades.DipoleMoment(1.54, "Debye")
    id = 39

    CP1 = {"ao": 0.0,
           "an": [120.529],
           "pow": [0.0801627],
           "ao_exp": [345.62, 289.038, 301.639],
           "exp": [2952.37, 734.653, 1593.55],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl linoleate of Huber et al. (2009).",
        "__doc__":  u"""Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, L.S., and Bruno, T.J. "Model for the Thermodynamic Properties of a Biodiesel Fuel," Energy & Fuels, 23:3790-3797, 2009.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [0.4282821e-1, 0.2443162e1, -0.3757540e1, -0.1588526, 0.4055990e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.36, 1.22, 1.45, 0.7],

        "nr2": [-0.1524090e1, -0.7686167, 0.1799950e1, -0.1590967e1, -0.1267681e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [3.0, 3.9, 2.2, 2.9, 1.25],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2198347e1, -0.7737211, -0.4314520],
        "d3": [1, 1, 3],
        "t3": [2.6, 3.0, 3.2],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.13378e2, 0.12258e2, -0.12205e2, -0.58118e1, -0.25451e1],
        "exp": [1.0, 1.5, 2.04, 4.3, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.54203, 0.13191e2, -0.20107e2, 0.11328e2, -0.60761],
        "exp": [0.18, 0.5, 0.7, 0.9, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.11612e2, 0.16300e3, -0.47913e3, 0.72986e3, -0.48202e3, -0.18198e3],
        "exp": [0.65, 1.78, 2.15, 2.7, 3.1, 9.8]}

    thermo0 = {"eq": 1,
             "__name__": "Perkins (2010)",
             "__doc__": """Perkins, R.A. and Huber, M.L., unpublished work, 2010.""",

             "Tref": 755.0, "kref": 1,
             "no": [-0.27125e-3, 0.259365e-2, 0.350241e-1, -0.902273e-2],
             "co": [0, 1, 2, 3],

             "Trefb": 755.0, "rhorefb": 0.897*M, "krefb": 1.,
             "nb": [-0.410106e-1, 0.328443e-1, -0.418506e-2, 0.0, 0.0,
                    0.606657e-1, -0.498407e-1, 0.121752e-1, 0.0, 0.0],
             "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
             "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
             "cb": [0]*10,

             "critical": 3,
             "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
             "Xio": 0.194e-9, "gam0": 0.0496, "qd": 8.75e-10, "Tcref": 1132.5}

    _thermal = thermo0,


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=C1Palmitate(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
