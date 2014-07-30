#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class Trans_2_butene(MEoS):
    """Ecuación de estado de multiparametros para el trans-2-buteno

    >>> buteno=Trans_2_butene(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.5f %0.5f %0.2f" % (buteno.T, buteno.rho, buteno.h.kJkg, buteno.s.kJkgK, buteno.cv.kJkgK, buteno.cp.kJkgK, buteno.w)
    300.0 2.31798 -4.715 0.00633 1.44222 1.61487 216.41
    """
    name="trans-butene "
    CASNumber="624-64-6"
    formula="CH3-CH=CH-CH3"
    synonym=""
    rhoc=unidades.Density(236.376)
    Tc=unidades.Temperature(428.61)
    Pc=unidades.Pressure(4027.3, "kPa")
    M=56.10632      #g/mol
    Tt=unidades.Temperature(167.6)
    Tb=unidades.Temperature(274.03)
    f_acent=0.21
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=26

    CP1={  "ao": 3.9988,
                "an": [], "pow": [],
                "ao_exp": [5.3276, 13.29, 9.6745, 0.40087],
                "exp": [362, 1603, 3729, 4527],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for trans-butene of Lemmon and Ihmels (2005)",
        "__doc__":  u"""Lemmon, E.W., Ihmels, E.C. Thermodynamic properties of the butenes Part II. Short fundamental equations of state. Fluid Phase Equilibria 228 – 229 (2004), 173 – 187.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1":  [0.81107, -2.8846, 1.0265, 0.016591, 0.086511, 0.00023256],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.22654, -0.072182, -0.24849, -0.071374, -0.024737, 0.011843],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq=helmholtz1,

    _vapor_Pressure={ "eq": 5, "ao": [-0.76226e1, 0.79421e1, -0.69631e1, -0.65517e1, 0.39584e1], "exp": [1.0, 1.5, 1.65, 4.8, 5.3]}
    _liquid_Density={ "eq": 1, "ao": [0.12452e2, -0.34419e2, 0.52257e2, -0.42889e2, 0.15463e2], "exp": [0.52, 0.73, 0.97, 1.24, 1.5]}
    _vapor_Density={ "eq": 3, "ao": [-0.31276e1, -0.60548e1, -0.18243e2, -0.60842e2, 0.13595e3, -0.18270e3], "exp": [0.412, 1.24, 3.2, 7.0, 10.0, 11.0]}

