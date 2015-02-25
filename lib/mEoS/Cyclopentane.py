#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Cyclopentane(MEoS):
    """Multiparameter equation of state for cyclopropane

    >>> cyc5=Cyclopentane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (cyc5.T, cyc5.rho, cyc5.u.kJkg, cyc5.h.kJkg, cyc5.s.kJkgK, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "cyclopropane"
    CASNumber = "287-92-3"
    formula = "C5H10"
    synonym = ""
    rhoc = unidades.Density(267.907678)
    Tc = unidades.Temperature(511.69)
    Pc = unidades.Pressure(4515.0, "kPa")
    M = 70.1329  # g/mol
    Tt = unidades.Temperature(179.722)
    Tb = unidades.Temperature(322.40)
    f_acent = 0.195
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 36

    CP1 = {"ao": 3.263,
           "an": [], "pow": [],
           "ao_exp": [2.151, 19.55, 14.45, 3.594],
           "exp": [179.0, 1336.0, 2911.0, 6420.0],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclopentane of Gedanitz et al. (2008).",
        "__doc__":  u"""Gedanitz, H., Davila, M.J., Lemmon, E.W.  unpublished equation, 2008.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 12.2, 
        "Pmin": 0.0089, "rhomin": 12.1, 

        "nr1": [0.4909331e-1, 0.1244679e1, -0.1990222e1, -0.5245596, 0.1764215],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.23, 0.94, 1.08, 0.53],

        "nr2": [-0.1066798e1, -0.5028152, 0.8484762, -0.4547443, -0.2767817e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.67, 1.8, 1.3, 2.5, 1.0],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.9455318, -0.3014822, -0.1675668, -0.637707],
        "d3": [1, 1, 3, 3],
        "t3": [0.87, 1.4, 2.4, 1.3],
        "alfa3": [1.023, 1.383, 0.996, 7.038],
        "beta3": [1.7, 1.55, 1.07, 87.17],
        "gamma3": [1.1, 0.64, 0.5, 1.26],
        "epsilon3": [0.713, 0.917, 0.688, 0.748]}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.64587e1, -0.12785e1, 0.35385e1, -0.24397e1, -0.29460e1],
        "exp": [1.0, 1.5, 1.83, 2.3, 4.3]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.27702e1, -0.17446, -0.11490e1, 0.13386e1, -0.36109e-1],
        "exp": [0.44, 0.7, 1.7, 2.0, 4.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.32453e1, -0.58511e1, -0.19268e2, -0.11351e3, 0.33499e3, -0.34655e3],
        "exp": [0.458, 1.22, 3.28, 7.7, 10.0, 11.0]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=Cyclopentane(T=300., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)

