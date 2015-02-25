#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class HexaC1_2Siloxane(MEoS):
    """Multiparameter equation of state for hexamethyldisiloxane

    >>> metilciclohexano=HexaC1_2Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "hexamethyldisiloxane"
    CASNumber = "107-46-0"
    formula = "C6H18OSi2"
    synonym = "MM"
    rhoc = unidades.Density(304.4043888253152)
    Tc = unidades.Temperature(518.7)
    Pc = unidades.Pressure(1939.0, "kPa")
    M = 162.37752  # g/mol
    Tt = unidades.Temperature(204.93)
    Tb = unidades.Temperature(373.401)
    f_acent = 0.418
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 39
    # id=1376

    CP1 = {"ao": 51.894,
           "an": [741.34e-3, -416.10e-6, 70.e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": 273.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 5.21, 
        "Pmin": 0.00269, "rhomin": 5.2, 

        "nr1": [1.01686012, -2.19713029, 0.75443188, -0.68003426, 0.19082162,
                0.10530133e-2],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.6284595, 0.30903042e-1, -0.83948727, -0.20262381,
                -0.35131597e-1, 0.25902341e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.86671e1, 0.11649e2, -0.11484e2, -0.53256e1],
        "exp": [1.0, 1.5, 1.65, 4.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.14533e2, -0.49804e2, 0.83748e2, -0.70321e2, 0.24283e2],
        "exp": [0.584, 0.8, 1.02, 1.26, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.35719e1, -0.14740e3, 0.40699e3, -0.69676e3, 0.12541e4, -0.91199e3],
        "exp": [0.373, 2.15, 2.6, 3.3, 4.2, 4.6]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=HexaC1_2Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
