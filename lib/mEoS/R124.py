#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R124(MEoS):
    """Multiparameter equation of state for R124

    >>> r124=R124(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r124.T, r124.rho, r124.u.kJkg, r124.h.kJkg, r124.s.kJkgK, r124.cv.kJkgK, r124.cp.kJkgK, r124.w)
    300.0 5.5897 286.68 304.57 1.5380 0.71520 0.77847 138.83
    """
    name = "1-chloro-1,2,2,2-tetrafluoroethane"
    CASNumber = "2837-89-0"
    formula = "CHClFCF3"
    synonym = "R124"
    rhoc = unidades.Density(560.)
    Tc = unidades.Temperature(395.425)
    Pc = unidades.Pressure(3624.295, "kPa")
    M = 136.4762  # g/mol
    Tt = unidades.Temperature(74.)
    Tb = unidades.Temperature(261.187)
    f_acent = 0.28810
    momentoDipolar = unidades.DipoleMoment(1.469, "Debye")
    id = 236
    # id = 1629

    CP1 = {"ao": 3.175638,
           "an": [14.77947/395.425, -5.2420986/395.425**2, 1.3381596/395.425**3],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 3.20532538,
           "an": [13.4403357/395.62, -2.32192933/395.62**2, -0.422826803/395.62**3],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-124 of de Vries et al. (1995).",
        "__doc__":  u"""de Vries, B., Tillner-Roth, R., and Baehr, H.D., "Thermodynamic Properties of HCFC 124," 19th International Congress of Refrigeration, The Hague, The Netherlands, International Institute of Refrigeration, IVa:582-589, 1995.""",
        "R": 8.314471,
        "cp": CP1,
        
        "Tmin": 120.0, "Tmax": 470.0, "Pmax": 40000.0, "rhomax": 13.5758, 
        "Pmin": 0.00000000032, "rhomin": 13.5758, 

        "nr1": [-0.1262962e-1, 0.2168373e1, -0.3330033e1, 0.1610361,
                -0.9666145e-4, 0.1191310e-1, -0.2880217e-2, 0.1681346e-2,
                0.1594968e-4],
        "d1": [1, 1, 1, 2, 2, 3, 5, 6, 8],
        "t1": [2, 0.5, 1, 0.5, 2.5, -1, 1, 0, -0.5],

        "nr2": [0.1289674, 0.1182213e-4, -0.4713997, -0.2412873, 0.6868066,
                -0.8621095e-1, 0.4728645e-5, 0.1487933e-1, -0.3001338e-1,
                0.1849606e-2, 0.4126073e-3],
        "d2": [2, 12, 1, 1, 1, 1, 15, 3, 3, 4, 9],
        "t2": [1.5, 1, 2.5, -0.25, 1, 5, 2, 15, 20, 15, 45],
        "c2": [1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4],
        "gamma2": [1]*11}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-124 of Younglove and McLinden (1994)",
        "__doc__":  u"""McLinden, M.O., Younglove, B.A., and Sandarusi, J., "Measurement of the PVT properties and formulation of an equation of state for refrigerant 124 (1-chloro-1,2,2,2-tetrafluoroethane)," 1994. (unpublished manuscript).""",
        "R": 8.314471,
        "cp": CP2,
        
        "Tmin": 120.0, "Tmax": 475.0, "Pmax": 36000.0, "rhomax": 13.98, 
        "Pmin": 0.00000000032, "rhomin": 14.54, 

        "b": [None, -0.195111839846e-1, 0.299978502039e1, -0.845849168162e2,
              0.146720754658e5, -0.232549336572e7, 0.938866046628e-3,
              -0.425069993257e1, 0.304859131600e4, 0.221314829910e7,
              -0.601971995213e-4, 0.100335188373e1, -0.468461812962e3,
              -0.927654315163e-2, -0.125426962519e2, -0.228534445089e4,
              0.168197835599e1, -0.537322295315e-1, 0.157915168095e2,
              -0.550297175283, -0.244349954189e7, -0.625153016263e8,
              -0.156149231820e6, 0.344268154495e10, -0.289212955106e4,
              0.108351996828e6, -0.404809912845e2, -0.220587292481e7,
              -0.564677367857, 0.175581172016e3, -0.762146322899e-3,
              -0.210617958917e1, 0.319236066221e2]}

    eq = helmholtz1, MBWR

    _surface = {"sigma": [0.05175], "exp": [1.197]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.75146e1, 0.37481e1, -0.30124e1, -0.37808e1, -0.53114],
        "exp": [1.0, 1.5, 1.68, 3.8, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.19127e1, 0.67778, -0.35129e-1, 0.30407, 0.69503e-1],
        "exp": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.28551e1, -0.63850e1, -0.17616e2, -0.37828e2, -0.23785e2, -0.13459e3],
        "exp": [0.388, 1.17, 3.0, 6.0, 8.0, 15.0]}
