#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class OctaC1_3Siloxane(MEoS):
    """Multiparamenter equation of state for octamethyltrisiloxane

    >>> metilciclohexano=OctaC1_3Siloxane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "octamethyltrisiloxane"
    CASNumber = "107-51-7"
    formula = "C8H24O2Si3"
    synonym = "MDM"
    rhoc = unidades.Density(256.7394094963634)
    Tc = unidades.Temperature(564.09)
    Pc = unidades.Pressure(1415.0, "kPa")
    M = 236.531  # g/mol
    Tt = unidades.Temperature(187.2)
    Tb = unidades.Temperature(425.66)
    f_acent = 0.529
    momentoDipolar = unidades.DipoleMoment(1.079, "Debye")
    id = 39
    # id = 1893

    CP1 = {"ao": 275.1,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [266040871.9, 0, 2051643622.0, 0],
           "hyp": [802.6, 0, 1829.6, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MDM of Colonna et al. (2006).",
        "__doc__":  u"""Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W., Multiparameter Equations of State for Selected Siloxanes, Fluid Phase Equilibria, 244:193-211, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.19735372, -2.40380622, 0.3256564, -0.19971259, 0.11206277,
                0.15893999e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.51234323, -0.20660361e-1, -0.38978114, -0.1186931,
                -0.37203537e-1, 0.18359984e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.85589e1, 0.20278e1, -0.28501e1, -0.64397e1, -0.85460e1],
        "exp": [1.0, 1.5, 2.3, 4.0, 13.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.54145, -0.27650e-1, 0.41558e1, -0.19104e1, 0.67606],
        "exp": [0.12, 0.36, 0.6, 0.8, 2.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.16483e1, -0.71410e1, -0.23088e2, -0.70554e2, 0.19938e1, -0.20193e3],
        "exp": [0.296, 0.905, 2.8, 5.9, 12.0, 13.0]}

    visco0 = {"eq": 5, "omega": 3,
              "__doc__": """T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E. "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties" Ind. Eng. Chem. Res. 1998, 27, 671-679""",
              "__name__": "Chung (1988)",
              "w": 0.531, "mur": 0.0, "k": 0.0}

    _viscosity = visco0,
#    _thermal=visco0,


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=OctaC1_3Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
