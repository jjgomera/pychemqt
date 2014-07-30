#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class COS(MEoS):
    """Ecuación de estado de multiparametros para el sulfuro de carbonilo

    >>> cos=COS(T=500, P=0.1)
    >>> print "%0.1f %0.4f %0.2f %0.2f %0.4f %0.5f %0.5f %0.2f" % (cos.T, cos.rho, cos.u.kJkg, cos.h.kJkg, cos.s.kJkgK, cos.cv.kJkgK, cos.cp.kJkgK, cos.w)
    500.0 1.4482 444.73 513.78 1.9750 0.67596 0.81579 288.36
    """
    name="carbonyl sulfide"
    CASNumber="463-58-1"
    formula="COS"
    synonym=""
    rhoc=unidades.Density(445.1565)
    Tc=unidades.Temperature(378.77)
    Pc=unidades.Pressure(6370.0, "kPa")
    M=60.0751      #g/mol
    Tt=unidades.Temperature(134.3)
    Tb=unidades.Temperature(222.99)
    f_acent=0.0978
    momentoDipolar=unidades.DipoleMoment(0.7152, "Debye")
    id=219

    CP1={  "ao": 3.5,
                "an": [], "pow": [],
                "ao_exp": [2.1651, 0.93456, 1.0623, 0.34269],
                "exp": [768, 1363, 3175, 12829],
                "ao_hyp": [], "hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbonyl sulfide of Lemmon and Span (2006)",
        "__doc__":  u"""Lemmon, E.W., Span, R. Short fundamental equations of state for 20 industrial fluids. J. Chem. Eng. Data 51 (2006), 785 – 850.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1":  [0.94374, -2.5348, 0.59058, -0.021488, 0.082083, 0.00024689],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.21226, -0.041251, -0.22333, -0.050828, -0.028333, 0.016983],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq=helmholtz1,

    _surface={"sigma": [ 0.06], "exp": [1.26]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.67055e1, 0.34248e1, -0.26677e1, -0.24717e1], "exp": [1., 1.5, 1.78, 4.8]}
    _liquid_Density={ "eq": 1, "ao": [0.76592e1, -0.19226e2, 0.27883e2, -0.23637e2, 0.99803e1], "exp": [0.515, 0.767, 1.034, 1.4, 1.7]}
    _vapor_Density={ "eq": 3, "ao": [-0.32494e1, -0.71460e1, 0.35026e2, -0.34039e2, -0.64206e2, -0.15225e3], "exp": [0.423, 1.464, 5.3, 4.1, 7.0, 17.0]}

