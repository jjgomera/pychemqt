#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R1234ze(MEoS):
    """Multiparameter equation of state for R1234ze

    >>> r227ea=R1234ze(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r227ea.T, r227ea.rho, r227ea.u.kJkg, r227ea.h.kJkg, r227ea.s.kJkgK, r227ea.cv.kJkgK, r227ea.cp.kJkgK, r227ea.w)
    300.0 7.0020 71.96 86.24 0.4638 0.75948 0.81568 122.13
    """
    name = "trans-1,3,3,3-tetrafluoropropene"
    CASNumber = "1645-83-6"
    formula = "CHF=CHCF3"
    synonym = "R-1234ze"
    rhoc = unidades.Density(489.238433112)
    Tc = unidades.Temperature(382.52)
    Pc = unidades.Pressure(3636.25, "kPa")
    M = 114.0415928  # g/mol
    Tt = unidades.Temperature(168.62)
    Tb = unidades.Temperature(254.2)
    f_acent = 0.313
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 671

    CP1 = {"ao": 6.259,
           "an": [], "pow": [],
           "ao_exp": [7.303, 8.597, 2.333], "exp": [691, 1705, 4216],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze of McLinden et al. (2010).",
        "__doc__":  u"""McLinden, M.O., Thol, M., and Lemmon, E.W. "Thermodynamic Properties of trans-1,3,3,3-Tetrafluoropropene [R1234ze(E)]: Measurements of Density and Vapor Pressure and a Comprehensive Equation of State," International Refrigeration and Air Conditioning Conference at Purdue, July 12-15, 2010.""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 420.0, "Pmax": 20000.0, "rhomax": 13.20, 
        "Pmin": 0.23, "rhomin": 13.19, 

        "nr1": [0.4434245e-1, 0.1646369e1, -0.2437488e1, -0.517056, 0.1815626],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.31, 0.923, 1.06, 0.44],

        "nr2": [-0.1210104e1, -0.5944653, 0.7521992, -0.6747216, -0.2448183e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.08, 2.32, 1.25, 2.0, 1.0],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.1379434e1, -0.4697024, -0.2036158, -0.8407447e-1, 0.5109529e-3],
        "d3": [1, 1, 3, 3, 2],
        "t3": [0.93, 1.93, 2.69, 1.0, 2.0],
        "alfa3": [1.0, 1.4, 1.134, 7.68, 24.],
        "beta3": [1.64, 1.57, 1.49, 257.0, 45.0],
        "gamma3": [1.13, 0.61, 0.65, 1.13, 1.34],
        "epsilon3": [0.711, 0.856, 0.753, 0.772, 1.88],
        "nr4": []}

    eq = helmholtz1,

    _surface = {"sigma": [0.05681], "exp": [1.23]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.76813e1, 0.31759e1, -0.26397e1, -0.35234e1],
        "exp": [1.0, 1.5, 1.8, 3.9]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.16130e1, 0.46976e1, -0.68759e1, 0.34227e1],
        "exp": [0.31, 0.94, 1.2, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.24897e1, -0.63324e1, -0.20262e2, -0.62612e2],
        "exp": [0.36, 1.07, 3.0, 6.8]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2010)",
               "__doc__": """Perkins, R.A. and Huber, M.L., unpublished work, 2010.""",

               "Tref": 382.52, "kref": 1,
               "no": [0.285145e-2, -0.439091e-2, 0.232616e-1],
               "co": [0, 1, 2],

               "Trefb": 367.85, "rhorefb": 4.17, "krefb": 1.,
               "nb": [-0.10750600e-2, -0.11560800e-1, 0.76230400e-2, 0.0, 0.0,
                      0.0, 0.10181100e-1, -0.22450100e-2, 0.0, 0.0],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.285e-10, "Tcref": 573.78}

    _thermal = thermo0,


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=R1234ze(T=300., rho=5.0)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
