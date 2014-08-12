#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class neoC5(MEoS):
    """Multiparameter equation of state for neopentane

    >>> neopentano=neoC5(T=500, P=50)
    >>> print "%0.1f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (neopentano.T, neopentano.rho, neopentano.h.kJkg, neopentano.s.kJkgK, neopentano.cv.kJkgK, neopentano.cp.kJkgK, neopentano.w)
    500.0 508.46 263.82 0.12415 2.5649 2.9723 795.9
    """
    name = "neopentane"
    CASNumber = "463-82-1"
    formula = "C(CH3)4"
    synonym = ""
    rhoc = unidades.Density(235.93)
    Tc = unidades.Temperature(433.74)
    Pc = unidades.Pressure(3196.0, "kPa")
    M = 72.14878  # g/mol
    Tt = unidades.Temperature(256.6)
    Tb = unidades.Temperature(282.65)
    f_acent = 0.1961
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 9

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [14.422, 12.868, 17.247, 12.663],
           "exp": [710, 1725, 3280, 7787],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [11.7618, -20.1101, 33.1688, 0],
           "hyp": [0.635392636*Tc, 1.977271641*Tc, 4.169371131*Tc, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for neopentane of Lemmon and Span (2006).",
        "__doc__":  u"""Lemmon, E.W. and Span, R., "Short Fundamental Equations of State for 20 Industrial Fluids," J. Chem. Eng. Data, 51:785-850, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.1136, -3.1792, 1.1411, -0.10467, 0.11754, 0.00034058],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.29553, -0.074765, -0.31474, -0.099401, -0.039569, 0.023177],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for neopentane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP2,

        "nr1": [-0.146552261671e1, 0.199230626557e1, -0.500821886276,
                0.119809758161e1, -0.363135896710e1, 0.312770556886e1,
                -0.237405105853e1, 0.473735725047, 0.101500881659, 0.184937708516,
                -0.290527628579e-1, -0.258919377284e-1, 0.748831217999e-1,
                0.216569936506e-1, -0.100375687935, 0.234924630013e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.146552261671e1, -0.199230626557e1, 0.500821886276,
                -0.834410647812, 0.262918341468e1, -0.188136966583e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.968832]*6}

    eq = helmholtz1, helmholtz2

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70262e1, 0.20090e1, -0.19932e1, -0.28503e1, -0.53760],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.56080e1, -0.13549e2, 0.29912e2, -0.28143e2, 0.89021e1],
        "exp": [0.45, 0.7, 1.0, 1.25, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.25177e1, -0.63565e1, -0.11985e3, 0.43740e3, -0.10749e4, 0.74007e3],
        "exp": [0.366, 1.14, 4.0, 5.0, 6.0, 6.5]}
