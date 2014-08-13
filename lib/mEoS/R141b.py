#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R141b(MEoS):
    """Multiparameter equation of state for R141b

    >>> r141B=R141b(T=500, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r141B.T, r141B.rho, r141B.u.kJkg, r141B.h.kJkg, r141B.s.kJkgK, r141B.cv.kJkgK, r141B.cp.kJkgK,r141B.w)
    500.0 2.8320 599.85 635.16 2.2987 0.92633 0.99970 194.56
    """
    name = "1,1-dichloro-1-fluoroethane"
    CASNumber = "1717-00-6"
    formula = "CCl2FCH3"
    synonym = "R141b"
    rhoc = unidades.Density(458.56)
    Tc = unidades.Temperature(477.5)
    Pc = unidades.Pressure(4212.0, "kPa")
    M = 116.94962  # g/mol
    Tt = unidades.Temperature(169.68)
    Tb = unidades.Temperature(305.20)
    f_acent = 0.2195
    momentoDipolar = unidades.DipoleMoment(2.014, "Debye")
    id = 236
    # id = 1633

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [6.8978, 7.8157, 3.2039], "exp": [502, 1571, 4603],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-141b of Lemmon and Span (2006)",
        "__doc__":  u"""Lemmon, E.W. and Span, R., "Short Fundamental Equations of State for 20 Industrial Fluids," J. Chem. Eng. Data, 51:785-850, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.1469, -3.6799, 1.3469, 0.083329, 0.00025137],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.32720, 0.46946, -0.029829, -0.31621, -0.026219, -0.078043,
                -0.020498],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [0.06087], "exp": [1.235]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73784e1, 0.52955e1, -0.46639e1, -0.31122e1, -0.18972e1],
        "exp": [1.0, 1.5, 1.7, 4.2, 9.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.10443e2, -0.24726e2, 0.27718e2, -0.11220e2, 0.75848],
        "exp": [0.49, 0.68, 0.88, 1.1, 2.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31177e1, -0.68872e1, -0.18566e2, -0.40311e2, -0.95472e1,
               -0.12482e3],
        "exp": [0.398, 1.33, 3.3, 6.7, 7.0, 14.0]}
