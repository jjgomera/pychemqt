#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R227ea(MEoS):
    """Multiparameter equation of state for R227ea

    >>> r227ea=R227ea(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r227ea.T, r227ea.rho, r227ea.u.kJkg, r227ea.h.kJkg, r227ea.s.kJkgK, r227ea.cv.kJkgK, r227ea.cp.kJkgK, r227ea.w)
    300.0 7.0020 71.96 86.24 0.4638 0.75948 0.81568 122.13
    """
    name = "1,1,1,2,3,3,3-heptafluoropropane"
    CASNumber = "431-89-0"
    formula = "CF3CHFCF3"
    synonym = "R227ea"
    rhoc = unidades.Density(594.25)
    Tc = unidades.Temperature(374.9)
    Pc = unidades.Pressure(2925.0, "kPa")
    M = 170.02886  # g/mol
    Tt = unidades.Temperature(146.35)
    Tb = unidades.Temperature(256.81)
    f_acent = 0.357
    momentoDipolar = unidades.DipoleMoment(1.456, "Debye")
    id = 671
    # id = 1872

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [11.43, 12.83], "exp": [403, 1428],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-227ea of Lemmon et al. (2007)",
        "__doc__":  u"""Lemmon, E.W., McLinden, M.O., and Meier, K. to be published in J. Chem. Eng. Data, 2007.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [2.024341, -2.605930, 0.4957216, -0.8240820, 0.06543703],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.34, 0.77, 0.36, 0.9, 1],

        "nr2": [-1.02461, .6247065, .2997521, -.353917, -1.232043, -.8824483],
        "d2": [1, 3, 6, 6, 2, 3],
        "t2": [2.82, 2.1, 0.9, 1.13, 3.8, 2.75],
        "c2": [1, 1, 1, 1, 2, 2],
        "gamma2": [1]*6,

        "nr3": [0.1349661, -0.2662928, 0.1764733, 0.01536163, -0.004667185,
                -11.70854, 0.9114512],
        "d3": [1, 2, 1, 1, 4, 2, 1],
        "t3": [1.5, 1.5, 2.5, 5.4, 4, 1, 3.5],
        "alfa3": [0.83, 2.19, 2.44, 3.65, 8.88, 8.23, 2.01],
        "beta3": [1.72, 5.2, 2.31, 1.02, 5.63, 50.9, 1.56],
        "gamma3": [0.414, 1.051, 1.226, 1.7, 0.904, 1.42, 0.926],
        "epsilon3": [1.13, 0.71, 1.2, 1.7, 0.546, 0.896, 0.747]}

    eq = helmholtz1,

    _surface = {"sigma": [0.048731, 0.016959, -0.023732], "exp": [1.26, 1.76, 2.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.77961e1, 0.21366e1, -0.26023e1, -0.57444e1, 0.23982e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.20032e1, 0.49235, 0.13738, 0.21057, -0.12834],
        "exp": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.2135e1, -.68425e1, -.21447e2, -.20457e3, .51795e3, -.45908e3],
        "exp": [0.324, 1.03, 3.0, 7.4, 9.0, 10.0]}
