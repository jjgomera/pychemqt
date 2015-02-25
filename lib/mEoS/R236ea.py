#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades
from R134a import R134a


class R236ea(MEoS):
    """Multiparameter equation of state for R236ea

    >>> c4f10=R236ea(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.4f %0.4f %0.2f" % (c4f10.T, c4f10.rho, c4f10.h.kJkg, c4f10.s.kJkgK, c4f10.cv.kJkgK, c4f10.cp.kJkgK, c4f10.w)
    300.0 0.08077 27.469 0.14617 10.7187 14.8450 1309.82
    """
    name = "1,1,1,2,3,3-hexafluoropropane"
    CASNumber = "431-63-0"
    formula = "CF3CHFCHF2"
    synonym = "R236ea"
    rhoc = unidades.Density(563.0044946256)
    Tc = unidades.Temperature(412.44)
    Pc = unidades.Pressure(3501.98, "kPa")
    M = 152.03928  # g/mol
    Tt = unidades.Temperature(170.0)
    Tb = unidades.Temperature(279.34)
    f_acent = 0.3794
    momentoDipolar = unidades.DipoleMoment(1.129, "Debye")
    id = 693

    CP1 = {"ao": 5.30694,
           "an": [0.03973, -1.859e-5], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    ecs = {"__type__": "ECS",
           "__name__": "Thermodynamic Extended Corresponding States model w/ T- and rho-dependent shape factors.",
           "__doc__":  u"""Huber, M.L. and Ely, J.F., "A predictive extended corresponding states model for pure and mixed refrigerants including an equation of state for R134a," Int. J. Refrigeration, 17:18-31, 1994.""",
           "cp": CP1,
           "ref": R134a,
           "eq": "helmholtz1",
           # "eq": "MBWR",
           "R": 8.314471,
        
            "Tmin": 242.0, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 10.465, 
#            "Pmin": aaaaaaa, "rhomin": aaaaaaa, 

           "ft": [-0.67786992, -0.52182651],
           "ft_add": [], "ft_add_exp": [],
           "fd": [0.113833347e-1], "fd_exp": [1],
           "ht": [0.142369159e1, 0.870214752e-1],
           "ht_add": [0.195298641e-1], "ht_add_exp": [1],
           "hd": [], "hd_exp": []}

    eq = ecs,

    _surface = {"sigma": [0.049561, 0.055607, -0.067899],
                "exp": [1.26, 1.76, 2.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.22360e2, 0.60938e3, -0.15037e4, 0.10657e4, -0.16142e3],
        "exp": [1, 1.5, 1.65, 1.8, 2.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.20433e1, -0.66050e1, 0.10613e2, -0.38994e1, 0.88965],
        "exp": [0.11, 0.3, 0.5, 0.8, 2.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.83044, -.99128e1, .1279e3, -.2739e4, .88175e4, -.85578e4],
        "exp": [0.08, 1.0, 5.0, 6.0, 7.0, 8.0]}

    trnECS = {"eq": "ecs",
              "__name__": "Extended Corresponding States model",
              "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., Model for the viscosity and thermal conductivity of refrigerants, including a new correlation for the viscosity of R134a, Ind.Eng.Chem.Res. 42: 3163-3178 (2003).""",

              "ref": R134a,
              "ref_eq": "helmholtz1",
              "eq_visco": "visco0",
              "eq_thermo": "thermo0",

              "f_int": [1.32e-3],
              "psi": [1.0],
              "phi": [1.0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.5e-9, "Tcref": 579.49}

#    _viscosity=trnECS,
#    _thermal=trnECS,

if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=R236ea(T=300., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
