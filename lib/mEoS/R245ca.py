#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades
from R134a import R134a


class R245ca(MEoS):
    """Multiparameter equation of state for R245ca

    >>> c4f10=R245ca(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.4f %0.4f %0.2f" % (c4f10.T, c4f10.rho, c4f10.h.kJkg, c4f10.s.kJkgK, c4f10.cv.kJkgK, c4f10.cp.kJkgK, c4f10.w)
    300.0 0.08077 27.469 0.14617 10.7187 14.8450 1309.82
    """
    name = "1,1,2,2,3-pentafluoropropane"
    CASNumber = "679-86-7"
    formula = "CHF2CF2CH2F"
    synonym = "R245ca"
    rhoc = unidades.Density(523.6)
    Tc = unidades.Temperature(447.57)
    Pc = unidades.Pressure(3925.0, "kPa")
    M = 134.04882  # g/mol
    Tt = unidades.Temperature(191.5)
    Tb = unidades.Temperature(298.28)
    f_acent = 0.3536
    momentoDipolar = unidades.DipoleMoment(1.740, "Debye")
    id = 693

    CP1 = {"ao": -3.8444,
           "an": [5.24008e-1, -3.74976e-4], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    ecs = {"__type__": "ECS",
           "__name__": "Thermodynamic Extended Corresponding States model w/ T- and rho-dependent shape factors.",
           "__doc__":  u"""Huber, M.L. and Ely, J.F., "A predictive extended corresponding states model for pure and mixed refrigerants including an equation of state for R134a," Int. J. Refrigeration, 17:18-31, 1994.""",
           "cp": CP1,
           "ref": R134a,
           "eq": "helmholtz1",
           "R": 8.314471,
        
            "Tmin": 200.0, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 11.995, 
#            "Pmin": aaaaaaa, "rhomin": aaaaaaa, 

           "ft": [-0.241011472, -0.788477331],
           "ft_add": [], "ft_add_exp": [],
           "fd": [], "fd_exp": [],
           "ht": [0.160567866e1, -0.727455038],
           "ht_add": [], "ht_add_exp": [],
           "hd": [], "hd_exp": []}

    eq = ecs,

    _surface = {"sigma": [0.069297, -0.022419], "exp": [1.2795, 3.1368]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.77617e1, 0.15867e1, -0.48984e2, 0.12055e3, -0.77763e2],
        "exp": [1.0, 1.5, 3.5, 3.8, 4.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-.57615e3, .20016e4, -.30523e4, .33777e4, -.29692e4, .12236e4],
        "exp": [0.79, 0.91, 1.1, 1.4, 1.7, 1.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [.67777, -.53374e1, -.1124e2, -.13244e2, .38054e2, -.10678e3],
        "exp": [0.082, 0.41, 2.0, 3.0, 5.0, 6.0]}

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

    cyc5=R245ca(T=300., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
