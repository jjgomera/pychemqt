#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class DecaC1Cyc5Siloxane(MEoS):
    """Ecuación de estado de multiparametros para el hexamethyldisiloxane

    >>> metilciclohexano=DecaC1Cyc5Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name="decamethylcyclopentasiloxane"
    CASNumber="541-02-6"
    formula="C10H30O5Si5"
    synonym="D5"
    rhoc=unidades.Density(292.570762680819)
    Tc=unidades.Temperature(619.23462341)
    Pc=unidades.Pressure(1160.0, "kPa")
    M=370.7697      #g/mol
    Tt=unidades.Temperature(226.0)
    Tb=unidades.Temperature(484.05)
    f_acent=0.658
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=39
#    id=1671

    CP1={ "ao": -34.898,
                "an": [1861.5e-3, -1403.4e-6, 500.0e-9],
                "pow": [1, 2, 3],
                "ao_exp": [],
                "exp": [],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.40844725, -2.29248044, 0.42851607, -0.73506382, 0.16103808, 0.29643278e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.82412481, 0.15214274, -0.68495890, -0.55703624e-1, 0.13055391e-1, -0.31853761e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq=helmholtz1, 

    _vapor_Pressure={ "eq": 5, "ao": [-0.99967e1, 0.70091e1, -0.72265e1, -0.62938e1], "exp": [1.0, 1.5, 1.87, 3.8]}
    _liquid_Density={ "eq": 1, "ao": [0.303988e3, -0.110342e4, 0.134359e4, -0.705243e3, 0.164540e3], "exp": [0.57, 0.65, 0.73, 0.84, 0.96]}
    _vapor_Density={ "eq": 3, "ao": [-0.37577e1, -0.47669e1, -0.24233e2, -0.29872e3, 0.34441e3, -0.32498e3], "exp": [0.459, 1.02, 2.6, 6.7, 7.7, 11.0]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=DecaC1Cyc5Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
