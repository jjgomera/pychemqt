#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class TetradecaC1_6Siloxane(MEoS):
    """Ecuación de estado de multiparametros para el hexamethyldisiloxane

    >>> metilciclohexano=TetradecaC1_6Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name="tetradecamethylhexasiloxane"
    CASNumber="107-52-8"
    formula="C14H42O5Si6"
    synonym="MD4M"
    rhoc=unidades.Density(285.6576532213632)
    Tc=unidades.Temperature(653.2)
    Pc=unidades.Pressure(877.0, "kPa")
    M=458.99328      #g/mol
    Tt=unidades.Temperature(214.15)
    Tb=unidades.Temperature(533.9)
    f_acent=0.836
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=39

    CP1={ "ao": -20.071,
                "an": [2228.5e-3, -1311.4e-6, 286.2e-9],
                "pow": [1, 2, 3],
                "ao_exp": [],
                "exp": [],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD4M of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.18492421, -1.87465636, -0.65713510e-1, -0.61812689, 0.19535804, 0.50678740e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],
 
        "nr2": [1.23544082, 0.49462708e-1, -0.73685283, -0.19991438, -0.55118673e-1, 0.28325885e-1 ],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq=helmholtz1, 

    _vapor_Pressure={ "eq": 5, "ao": [-0.10532e2, 0.33939e1, -0.89744e1, -0.56150e1], "exp": [1.0, 1.5, 2.75, 5.1]}
    _liquid_Density={ "eq": 1, "ao": [0.10453e1, 0.55476, 0.44536e1, -0.76179e1, 0.46237e1], "exp": [0.235, 0.6, 0.95, 1.35, 1.7]}
    _vapor_Density={ "eq": 3, "ao": [-0.10890e1, -0.84374e1, -0.35615e2, -0.73478e3, 0.19915e4, -0.16317e4], "exp": [0.231, 0.8, 2.9, 7.7, 9.0, 10.0]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=TetradecaC1_6Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
