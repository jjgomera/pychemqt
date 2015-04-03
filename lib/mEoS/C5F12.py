#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades
from R113 import R113


class C5F12(MEoS):
    """Multiparameter equation of state for perfluropentane

    >>> c4f10=C5F12(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.4f %0.4f %0.2f" % (c4f10.T, c4f10.rho, c4f10.h.kJkg, c4f10.s.kJkgK, c4f10.cv.kJkgK, c4f10.cp.kJkgK, c4f10.w)
    300.0 0.08077 27.469 0.14617 10.7187 14.8450 1309.82
    """
    name = "perfluoropentane"
    CASNumber = "678-26-2"
    formula = "678-26-2"
    synonym = ""
    rhoc = unidades.Density(609.47148)
    Tc = unidades.Temperature(420.555)
    Pc = unidades.Pressure(2045.0, "kPa")
    M = 288.03  # g/mol
    Tt = unidades.Temperature(148.363)
    Tb = unidades.Temperature(302.904)
    f_acent = 0.423
    momentoDipolar = unidades.DipoleMoment(0.14, "Debye")
    id = 693

    CP1 = {"ao": 0.24705743e1,
           "an": [0.11875895, -0.12235660e-3, 0.45790525e-7], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    ecs = {
        "__type__": "ECS",
        "__name__": "Thermodynamic Extended Corresponding States model w/ T- and rho-dependent shape factors.",
        "__doc__":  u"""Huber, M.L. and Ely, J.F., "A predictive extended corresponding states model for pure and mixed refrigerants including an equation of state for R134a," Int. J. Refrigeration, 17:18-31, 1994.""",
        "cp": CP1,
        "ref": R113,
        "eq": "helmholtz1",
        "R": 8.314471,

        "Tmin": Tt, "Tmax": 500., "Pmax": 30000.0, "rhomax": 6.7, 
#        "Pmin": 0.61166, "rhomin": 55.497, 

        "ft": [0.960871894e-2, -0.820122088],
        "ft_add": [], "ft_add_exp": [],
        "fd": [], "fd_exp": [],
        "ht": [-0.367946699e-1, 0.736529816e-1],
        "ht_add": [], "ht_add_exp": [],
        "hd": [], "hd_exp": []}

    eq = ecs,

    _surface = {"sigma": [0.04394], "exp": [1.254]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.86275e1, 0.13120e2, -0.12690e2, -0.63017e1, -0.65096e2],
        "exp": [1.0, 1.5, 1.6, 4.0, 15.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.44087, 0.10672e2, -0.95203e1, 0.18688e1, -0.31782e1],
        "exp": [0.09, 0.75, 0.9, 2.0, 8.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.57928, -0.60202e1, -0.16560e2, -0.67528e2, -0.55108e3],
        "exp": [0.13, 0.7, 2.3, 5.5, 15.0]}

    trnECS = {"eq": "ecs",
              "__name__": "Extended Corresponding States model",
              "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., Model for the viscosity and thermal conductivity of refrigerants, including a new correlation for the viscosity of R134a, Ind.Eng.Chem.Res. 42: 3163-3178 (2003).""",

              "ref": R113,
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

    cyc5=C5F12(T=400., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)

