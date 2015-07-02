#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class CF3I(MEoS):
    """Multiparameter equation of state for trifluoroiodomethane

    >>> trifluoroiodometano=CF3I(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (trifluoroiodometano.T, trifluoroiodometano.rho, trifluoroiodometano.u.kJkg, trifluoroiodometano.h.kJkg, trifluoroiodometano.s.kJkgK, trifluoroiodometano.cv.kJkgK, trifluoroiodometano.cp.kJkgK, trifluoroiodometano.w)
    300.0 8.02849 7.17 19.624 0.09457 0.3043 0.3514 118.60
    """
    name = "trifluoroiodomethane"
    CASNumber = "2314-97-8"
    formula = "CF3I"
    synonym = ""
    rhoc = unidades.Density(868.)
    Tc = unidades.Temperature(396.44)
    Pc = unidades.Pressure(3953., "kPa")
    M = 195.9103796  # g/mol
    Tt = unidades.Temperature(120.)
    Tb = unidades.Temperature(251.3)
    f_acent = 0.18
    momentoDipolar = unidades.DipoleMoment(0.92, "Debye")
    id = 645

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [6.2641], "exp": [694.1467],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for CF3I of McLinden and Lemmon (2013)",
        "__doi__": {"autor": "McLinden, M.O. and Lemmon, E.W.",
                    "title": "Thermodynamic Properties of R-227ea, R-365mfc, R-115, and R-13I1", 
                    "ref": "to be submitted to J. Chem. Eng. Data, 2013.",
                    "doi": ""}, 

        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR", 

        "Tmin": Tt, "Tmax": 420., "Pmax": 20000.0, "rhomax": 14.1, 
        "Pmin": 0.0004623, "rhomin": 14.05, 

        "nr1": [0.112191e1, -0.308087e1, 0.111307e1, -0.184885, 0.110971,
                0.325005e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.23, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.333357, -0.288288e-1, -0.371554, -0.997985e-1, -0.333205e-1,
                0.207882e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.05767], "exp": [1.298]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.8642, 1.7877, -1.0619, -2.1677],
        "exp": [1.0, 1.5, 1.9, 3.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [2.0711, 1.562, -2.599, 1.7177],
        "exp": [0.38, 1.3, 1.9, 2.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.0987, -6.8771, -19.701, -46.86, -100.02],
        "exp": [0.41, 1.33, 3.5, 7.4, 16.0]}
