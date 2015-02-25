#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R507A(MEoS):
    """Multiparameter equation of state for R507A
    50% R125, 50% R143a

    >>> r507A=R507A(T=300, P=0.1)
    >>> print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w)
    500.0 2.3837 204.140 0.475 1.07947 1.1649 212.51
    """
    name = "R507A"
    CASNumber = ""
    formula = "R125+R143a"
    synonym = "R507A"
    rhoc = unidades.Density(490.7370688)
    Tc = unidades.Temperature(343.765)
    Pc = unidades.Pressure(3704.9, "kPa")
    M = 98.8592  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(226.41)
    f_acent = 0.286
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 62
    # id = None

    CP1 = {"ao": 0.,
           "an": [1.5680], "pow": [0.25],
           "ao_exp": [0.95006, 4.1887, 5.5184], "exp": [364, 815, 1768],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-404A of Lemmon (2003)",
        "__doc__":  u"""Pseudo Pure-Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C," Int. J. Thermophys., 24(4):991-1006, 2003.""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 14.13, 
        "Pmin": 23.23, "rhomin": 14.13, 

        "nr1": [0.624982e1, -0.807855e1, 0.264843e-1, 0.286215, -0.507076e-2,
                0.109552e-1, 0.116124e-2],
        "d1": [1, 1, 1, 2, 2, 4, 6],
        "t1": [0.692, 0.943, 5.8, 0.77, 5.84, 0.24, 0.69],

        "nr2": [0.138469e1, -0.922473, -0.503562e-1, 0.822098, -0.277727,
                0.358172, -0.126426e-1, -0.607010e-2, -.815653e-1, -.233323e-1,
                .352952e-1, .159566e-1, .755927e-1, -.542007e-1, .170451e-1],
        "d2": [1, 1, 1, 2, 2, 3, 4, 7, 2, 3, 4, 4, 2, 3, 5],
        "t2": [2, 3, 7, 2.2, 4.3, 2.7, 1.2, 1.23, 12, 6, 8.5, 11.5, 13, 17, 16.2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*15}

    eq = helmholtz1,

    _surface = {"sigma": [0.06701, -0.04297], "exp": [1.3066, 2.3145]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.5459, 2.338, -2.237, -4.1535],
        "exp": [1, 1.5, 2.1, 4.7]}
    _liquid_Pressure = {
        "eq": 5,
        "ao": [-7.4853, 2.0115, -2.0141, -3.7763],
        "exp": [1, 1.5, 2., 4.6]}
