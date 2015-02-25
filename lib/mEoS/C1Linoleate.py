#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class C1Linoleate(MEoS):
    """Multiparameter equation of state for methyl linoleate

    >>> metilciclohexano=C1Linoleate(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "methyl linoleate"
    CASNumber = "112-63-0"
    formula = "C19H34O2"
    synonym = ""
    rhoc = unidades.Density(238.051213304)
    Tc = unidades.Temperature(799.0)
    Pc = unidades.Pressure(1341.0, "kPa")
    M = 294.47206  # g/mol
    Tt = unidades.Temperature(238.1)
    Tb = unidades.Temperature(628.84)
    f_acent = 0.805
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 39

    CP1 = {"ao": 0.0,
           "an": [190.986],
           "pow": [0.020213],
           "ao_exp": [437.371, 287.222, 321.956],
           "exp": [3052.11, 746.631, 1624.33],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl linoleate of Huber et al. (2009).",
        "__doc__":  u"""Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, L.S., and Bruno, T.J. "Model for the Thermodynamic Properties of a Biodiesel Fuel," Energy & Fuels, 23:3790-3797, 2009.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 3.16, 
        "Pmin": 0.7e-14, "rhomin": 3.16, 

        "nr1": [0.3183187e-1, 0.1927286e1, -0.3685053e1, 0.8449312e-1],
        "d1": [4, 1, 1, 3],
        "t1": [1, 0.2, 1.2, 1.0],

        "nr2": [-0.9766643, -0.4323178, 0.2000470e1, -0.1752030e1, -0.1726895e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.2, 2.5, 1.8, 1.92, 1.47],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2116515e1, -0.7884271, -0.3811699],
        "d3": [1, 1, 3],
        "t3": [1.7, 2.3, 2.1],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.10946e2, 0.48849e1, -0.46773e1, -0.80201e1, -0.89572e1],
        "exp": [1.0, 1.5, 2.22, 3.6, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.22705e3, -0.66763e3, 0.72323e3, -0.49244e3, 0.21391e3],
        "exp": [0.83, 0.98, 1.17, 1.5, 1.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.85880e1, 0.14766e2, -0.24195e2, -0.37474e3, 0.32689e3, -0.19125e3],
        "exp": [0.568, 1.08, 1.4, 4.8, 5.0, 9.0]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2010)",
               "__doc__": """Perkins, R.A. and Huber, M.L., unpublished work, 2010.""",

               "Tref": 799.0, "kref": 1,
               "no": [-0.10904200e-3,  0.24054300e-2,  0.40736400e-1, -0.10592800e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 799.0, "rhorefb": 0.8084*M, "krefb": 1.,
               "nb": [-0.713126e-1, 0.466421e-1, -0.557406e-2, 0.0, 0.0,
                      0.989415e-1, -0.65785e-1, 0.128922e-1, 0.0, 0.0],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 8.75e-10, "Tcref": 1198.5}

    thermo1 = {"eq": 5, "omega": 3,
               "__doc__": """T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E. "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties" Ind. Eng. Chem. Res. 1998, 27, 671-679""",
               "__name__": "Chung (1988)",
               "w": 0.805, "mur": 0.0, "k": 0.0}

    _thermal = thermo0, thermo1


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=C1Linoleate(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
