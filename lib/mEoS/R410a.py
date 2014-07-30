#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class R410A(MEoS):
    """Ecuación de estado de multiparametros para el R410A
    50% R32, 50% R125

    >>> aire=R410A(T=500, P=0.1)
    >>> print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w)
    500.0 1.7492 191.122 0.449 0.98159 1.0975 252.59
    """
    name="R410A"
    CASNumber=""
    formula="R32+R125"
    synonym="R410A"
    rhoc=unidades.Density(459.0300696)
    Tc=unidades.Temperature(344.494)
    Pc=unidades.Pressure(4901.2, "kPa")
    M=72.5854       #g/mol
    Tt=unidades.Temperature(200.0)
    Tb=unidades.Temperature(221.71)
    f_acent=0.296
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=62

    CP1={  "ao": 0.,
                "an": [2.8749],
                "pow": [0.1],
                "ao_exp": [2.0623, 5.9751, 1.5612],
                "exp": [697, 1723, 3875],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-410A of Lemmon (2003)",
        "__doc__":  u"""Pseudo Pure-Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C," Int. J. Thermophys., 24(4):991-1006, 2003.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [0.987252, -0.103017e1, 0.117666e1, -0.138991, 0.302373e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.44, 1.2, 2.97, 2.95, 0.2],

        "nr2": [-0.253639e1, -0.196680e1, -0.830480, 0.172477, -0.261116, -0.745473e-1, 0.679757, -0.652431, 0.553849e-1, -0.710970e-1, -0.875332e-3, 0.200760e-1, -0.139761e-1, -0.185110e-1, 0.171939e-1, -0.482049e-2],
        "d2": [1, 2, 3, 5, 5, 5, 1, 1, 4, 4, 9, 2, 2, 4, 5, 6],
        "t2": [1.93, 1.78, 3., 0.2, 0.74, 3, 2.1, 4.3, 0.25, 7, 4.7, 13, 16, 25, 17, 7.4],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3],
        "gamma2": [1]*16}

    eq=helmholtz1,

    _surface={"sigma": [0.06443], "exp": [1.245]}
    _vapor_Pressure={ "eq": 5, "ao": [-7.4411, 1.9883, -2.4924, -3.2633], "exp": [1, 1.6, 2.4, 5]}
    _liquid_Pressure={ "eq": 5, "ao": [-7.2818, 2.5093, -3.2695, -2.8022], "exp": [1, 1.8, 2.4, 4.9]}


