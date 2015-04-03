#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R236fa(MEoS):
    """Multiparameter equation of state for R236fa

    >>> r236fa=R236fa(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r236fa.T, r236fa.rho, r236fa.u.kJkg, r236fa.h.kJkg, r236fa.s.kJkgK, r236fa.cv.kJkgK, r236fa.cp.kJkgK, r236fa.w)
    300.0 7.6962 302.73 315.72 1.4786 0.75013 0.79926 116.42
    """
    name = "1,1,1,3,3,3-hexafluoropropane"
    CASNumber = "690-39-1"
    formula = "CF3CH2CF3"
    synonym = "R236fa"
    rhoc = unidades.Density(551.2945)
    Tc = unidades.Temperature(398.07)
    Pc = unidades.Pressure(3200.0, "kPa")
    M = 152.0393  # g/mol
    Tt = unidades.Temperature(179.52)
    Tb = unidades.Temperature(271.71)
    f_acent = 0.37721
    momentoDipolar = unidades.DipoleMoment(1.982, "Debye")
    id = 671
    # id = 1873

    CP1 = {"ao": 53.4662555,
           "an": [1, 2], "pow": [0.228092134, 0.352999168e-4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-236fa of Outcalt and McLinden (1995).",
        "__doc__":  u"""Outcalt, S.L. and McLinden, M.O., "An equation of state for the thermodynamic properties of R236fa," NIST report to sponsor (U.S. Navy, David Taylor Model Basin) under contract N61533-94-F-0152, 1995.""",
        "R": 8.314471,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 40000.0, "rhomax": 11.30, 
        "Pmin": 0.162, "rhomin": 11.29, 

        "b": [None, -0.661121874831e-1, 0.861763902745e1, -0.233732255968e3,
              0.437486232843e5, -0.539677761508e7, -0.757588552002e-2,
              0.107379563512e2, -0.106626588551e5, -0.103047455432e6,
              -0.194868091617e-2, 0.438365228107e1, -0.111207843880e4,
              -0.263710051508, 0.477521163113e2, 0.197804035098e4,
              -0.485710898935e1, 0.144821196401, -0.221059322936e2,
              0.926270169913, 0.577920666161e7, -0.985511065626e9,
              0.197199808018e6, 0.319420123094e10, 0.792946107314e4,
              -0.693606295610e6, 0.849836259084e2, 0.209702051124e7,
              0.110600369167e1, 0.953714711849e2, -0.881815206562e-2,
              0.973194908842e1, -0.935516922205e3]}

    eq = MBWR,

    _surface = {"sigma": [0.05389], "exp": [1.249]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.78785e1, 0.15884e1, -0.48864e1, -0.50273e1, 0.89900e1],
        "exp": [1.0, 1.5, 3.1, 8.0, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.12320e2, -0.27579e2, 0.40114e2, -0.35461e2, 0.13769e2],
        "exp": [0.579, 0.77, 0.97, 1.17, 1.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.44507e1, -.37583e1, -.20279e2, -.26801e3, .50171e3, -.34917e3],
        "exp": [0.506, 1.16, 2.8, 7.0, 8.0, 9.0]}
