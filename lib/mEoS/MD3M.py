#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class MD3M(MEoS):
    """Multiparameter equation of state for dodecamethylpentasiloxane"""
    name = "dodecamethylpentasiloxane"
    CASNumber = "141-63-9"
    formula = "C12H36Si5O4"
    synonym = "MD3M"
    rhoc = unidades.Density(263.9218791237794)
    Tc = unidades.Temperature(628.36)
    Pc = unidades.Pressure(945.0, "kPa")
    M = 384.839  # g/mol
    Tt = unidades.Temperature(192.0)
    Tb = unidades.Temperature(503.03)
    f_acent = 0.7218
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 39

    CP1 = {"ao": 463.2,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [609372332.2, 0, 4290277999.0, 0],
           "hyp": [908.5, 0, 2117.1, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD3M of Colonna et al. (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., and Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and [O-Si-(CH3)2]6", 
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi":  "10.1016/j.fluid.2007.10.001"}, 
            
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.54, 
        "Pmin": 0.4e-12, "rhomin": 2.54, 

        "nr1": [1.20540386, -2.42914797, 0.69016432, -0.69268041, 0.18506046,
                0.31161436e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.99862519, 0.74229034e-1, -0.80259136, -0.20865337,
                -0.36461791e-1, 0.19174051e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.92608e1, 0.15861e1, -0.32859e1, -0.75194e1, -0.34883e1],
        "exp": [1.0, 1.5, 2.46, 3.7, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.74156, 0.21723e1, 0.66412e2, -0.17125e3, 0.10848e3],
        "exp": [0.22, 0.51, 5.5, 6.0, 6.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.19054e1, -0.74526e1, -0.10520e3, 0.24548e3, -0.23783e3,
               -0.21226e3],
        "exp": [0.332, 0.88, 3.25, 4.0, 4.6, 12.0]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=DodecaC1_5Siloxane(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
