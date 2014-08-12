#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class DodecaC1Cyc6Siloxane(MEoS):
    """Multiparameter equation of state for dodecamethylcyclohexasilosane

    >>> metilciclohexano=DodecaC1Cyc6Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "dodecamethylcyclohexasiloxane"
    CASNumber = "540-97-6"
    formula = "C12H36Si6O6"
    synonym = "D6"
    rhoc = unidades.Density(279.0957298413672)
    Tc = unidades.Temperature(645.78)
    Pc = unidades.Pressure(961.0, "kPa")
    M = 444.924  # g/mol
    Tt = unidades.Temperature(270.2)
    Tb = unidades.Temperature(518.11)
    f_acent = 0.736
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 39
    # id=1674

    CP1 = {"ao": 468.7,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [425104546.6, 0, 3151243909.0, 0],
           "hyp": [786.8, 0, 1792.1, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.69156186, -3.37962568, 0.38609039, 0.64598995e-1,
                0.10589012, 0.45456825e-4],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.74169279, -0.88102648e-1, -0.17373336, -0.10951368,
                -0.62695695e-1, 0.37459986e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.96557e1, 0.62155, 0.17863e1, -0.10496e2, -0.84102e1],
        "exp": [1.0, 1.5, 1.72, 3.18, 11.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.42563e2, -0.15707e3, 0.29502e3, -0.24191e3, 0.65145e2],
        "exp": [0.537, 0.68, 0.85, 1.0, 1.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.20930e1, -0.94442e1, -0.44731e2, -0.57898e2, -0.35144e2,
               -0.29661e3],
        "exp": [0.338, 1.02, 3.46, 7.1, 7.4, 15.0]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=DodecaC1Cyc6Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
