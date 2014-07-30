#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class DecaC1_4Siloxane(MEoS):
    """Ecuación de estado de multiparametros para el hexamethyldisiloxane

    >>> metilciclohexano=DecaC1_4Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name="decamethyltetrasiloxane"
    CASNumber="141-62-8"
    formula="C10H30Si4O3"
    synonym="MD2M"
    rhoc=unidades.Density(284.1716396620275)
    Tc=unidades.Temperature(599.40)
    Pc=unidades.Pressure(1227.0, "kPa")
    M=310.685      #g/mol
    Tt=unidades.Temperature(205.2)
    Tb=unidades.Temperature(467.51)
    f_acent=0.668
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=39
#    id=1837

    CP1={ "ao": 331.9,
                "an": [],
                "pow": [],
                "ao_exp": [],
                "exp": [],
                "ao_hyp": [329620742.8, 0, 2556558319.0, 0],"hyp": [795.1, 0, 1813.8, 0]}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD2M of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.33840331, -2.62939393, 0.4398383, -0.53496715, 0.1818844, 0.40774609e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],
  
        "nr2": [1.13444506, 0.5774631e-1, -0.5917498, -0.11020225, -0.34942635e-1, 0.7646298e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq=helmholtz1, 

    _vapor_Pressure={ "eq": 5, "ao": [-0.10029e2, 0.44434e1, -0.36753e1, -0.68925e1, -0.32211e1], "exp": [1.0, 1.5, 2.06, 3.5, 10.0]}
    _liquid_Density={ "eq": 1, "ao": [0.12608e2, -0.32120e2, 0.33559e2, -0.11695e2, 0.76192], "exp": [0.48, 0.64, 0.8, 1.0, 2.6]}
    _vapor_Density={ "eq": 3, "ao": [-0.24684e1, -0.71262e1, -0.27601e2, -0.49458e2, -0.24106e2, -0.19370e3], "exp": [0.376, 0.94, 2.9, 5.9, 6.2, 13.0]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=DecaC1_4Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
