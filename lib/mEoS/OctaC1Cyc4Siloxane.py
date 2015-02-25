#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class OctaC1Cyc4Siloxane(MEoS):
    """Multiparameter equation of state for octamethylcyclotetrasiloxane

    >>> metilciclohexano=OctaC1Cyc4Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "octamethylcyclotetrasiloxane"
    CASNumber = "556-67-2"
    formula = "C8H24O4Si4"
    synonym = "D4"
    rhoc = unidades.Density(307.0335906736056)
    Tc = unidades.Temperature(586.49)
    Pc = unidades.Pressure(1332.0, "kPa")
    M = 296.61576  # g/mol
    Tt = unidades.Temperature(290.25)
    Tb = unidades.Temperature(448.504)
    f_acent = 0.592
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 39
    # id=1430

    CP1 = {"ao": -18.256,
           "an": [1427.2e-3, -990.20e-6, 300.0e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for octamethylcyclotetrasiloxane of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": 300.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.21, 
        "Pmin": 0.0696, "rhomin": 3.2, 

        "nr1": [1.05392408, -2.22981918, 0.77573923, -0.6937405, 0.18721557,
                0.42193330e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.70301835, 0.47851888e-1, -0.8025348, -0.18968872,
                -0.22211781e-1, 0.60103354e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.87935e1, 0.27204e1, -0.48174e1, -0.69086e1],
        "exp": [1.0, 1.5, 2.2, 4.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.14563e1, -0.94215, 0.45065e1, -0.27688e1, 0.8745],
        "exp": [0.24, 0.5, 0.75, 1.0, 2.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.16204e1, -0.57888e1, -0.24291e2, 0.53567e2, -0.12135e3, -0.10976e4],
        "exp": [0.31, 0.78, 2.5, 4.4, 5.0, 15.0]}

    visco0 = {"eq": 5, "omega": 3,
              "__doc__": """T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E. "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties" Ind. Eng. Chem. Res. 1998, 27, 671-679""",
              "__name__": "Chung (1988)",
              "w": 0.592, "mur": 0.0, "k": 0.0}

    _viscosity = visco0,
#    _thermal=visco0,


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=OctaC1Cyc4Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
