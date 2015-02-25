#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class oH2(MEoS):
    """Multiparamente equation of state for hydrogen (orto)

    >>> hidrogeno=oH2(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.4f %0.4f %0.2f" % (hidrogeno.T, hidrogeno.rho, hidrogeno.h.kJkg, hidrogeno.s.kJkgK, hidrogeno.cv.kJkgK, hidrogeno.cp.kJkgK, hidrogeno.w)
    300.0 0.08077 26.133 0.14170 10.7414 14.8676 1309.43
    """
    name = "ortohydrogen"
    CASNumber = "1333-74-0o"
    formula = "H2"
    synonym = "R-702o"
    rhoc = unidades.Density(31.1362)
    Tc = unidades.Temperature(33.32)
    Pc = unidades.Pressure(1310.65, "kPa")
    M = 2.01594  # g/mol
    Tt = unidades.Temperature(14.008)
    Tb = unidades.Temperature(20.4)
    f_acent = -0.219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [2.54151, -2.3661, 1.00365, 1.22447],
           "exp": [856, 1444, 2194, 6968],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ortohydrogen of Leachman et al. (2007)",
        "__doc__":  u"""Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W. Fundamental equations of state for parahydrogen, normal hydrogen, and orthohydrogen. J. Phys. Chem. Ref. Data, 38 (2009), 721 – 748.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 2000000.0, "rhomax": 38.2, 
        "Pmin": 7.461, "rhomin": 38.2, 

        "nr1": [-6.83148, 0.01, 2.11505, 4.38353, 0.211292, -1.00939, 0.142086],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.7333, 1, 1.1372, 0.5136, 0.5638, 1.6248, 1.829],

        "nr2": [-0.87696, 0.804927],
        "d2": [1, 3],
        "t2": [2.404, 2.105],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.710775, 0.0639688, 0.0710858, -0.087654, 0.647088],
        "d3": [2, 1, 3, 1, 1],
        "t3": [4.1, 7.658, 1.259, 7.589, 3.946],
        "alfa3": [1.169, 0.894, 0.04, 2.072, 1.306],
        "beta3": [0.4555, 0.4046, 0.0869, 0.4415, 0.5743],
        "gamma3": [1.5444, 0.6627, 0.763, 0.6587, 1.4327],
        "epsilon3": [0.6366, 0.3876, 0.9437, 0.3976, 0.9626],
        "nr4": []}

    eq = helmholtz1,

    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [2.0297, 0.0069], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [0.181, 0.021, -7.4],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 3]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.488684e1, 0.105310e1, 0.856947, -0.185355],
        "exp": [1.0, 1.5, 2.7, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.43911e1, -0.75872e1, 0.10402e2, -0.72651e1, 0.18302e1],
        "exp": [0.53, 0.93, 1.35, 1.8, 2.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31463e1, -0.16183e2, 0.31803e2, -0.21961e3, 0.43123e3, -0.25591e3],
        "exp": [0.491, 2.1, 2.9, 4.4, 5.0, 5.5]}
