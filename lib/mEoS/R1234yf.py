#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R1234yf(MEoS):
    """Multiparameter equation of state for R1234yf

    >>> r227ea=R1234yf(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r227ea.T, r227ea.rho, r227ea.u.kJkg, r227ea.h.kJkg, r227ea.s.kJkgK, r227ea.cv.kJkgK, r227ea.cp.kJkgK, r227ea.w)
    300.0 7.0020 71.96 86.24 0.4638 0.75948 0.81568 122.13
    """
    name = "2,3,3,3-tetrafluoropropene"
    CASNumber = "754-12-1"
    formula = "CF3CF=CH2"
    synonym = "R-1234yf"
    rhoc = unidades.Density(475.553441976)
    Tc = unidades.Temperature(367.85)
    Pc = unidades.Pressure(3382.2, "kPa")
    M = 114.0415928  # g/mol
    Tt = unidades.Temperature(220.)
    Tb = unidades.Temperature(243.7)
    f_acent = 0.276
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 671

    CP1 = {"ao": 5.944,
           "an": [], "pow": [],
           "ao_exp": [7.549, 1.537, 2.03, 7.455],
           "exp": [718, 877, 4465, 1755],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234yf of Richter et al. (2010)",
        "__doc__":  u"""Richter, M., McLinden, M.O., and Lemmon, E.W. "Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf): p-rho-T Measurements and an Equation of State," submitted to J. Chem. Eng. Data, 2010.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [0.4592563e-1, 0.1546958e1, -0.2355237e1, -0.4827835, 0.1758022],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.32, 0.929, 0.94, 0.38],

        "nr2": [-0.1210006e1, -0.6177084, 0.6805262, -0.6968555, -0.2695779e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.28, 1.76, 0.97, 2.44, 1.05],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*7,

        "nr3": [0.1389966e1, -0.4777136, -0.1975184, -0.1147646e1, 0.3428541e-3],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.4, 3.0, 3.5, 1.0, 3.5],
        "alfa3": [1.02, 1.336, 1.055, 5.84, 16.2],
        "beta3": [1.42, 2.31, 0.89, 80., 108.],
        "gamma3": [1.13, 0.67, 0.46, 1.28, 1.2],
        "epsilon3": [0.712, 0.910, 0.677, 0.718, 1.64],
        "nr4": []}

    eq = helmholtz1,

    _surface = {"sigma": [0.05983], "exp": [1.367]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74697e1, 0.27915e1, -0.21312e1, -0.29531e1],
        "exp": [1.0, 1.5, 1.8, 3.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.19083e1, -0.21383e1, 0.93653e1, -0.98659e1, 0.35859e1],
        "exp": [0.32, 0.56, 0.8, 1.0, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.23511e1, -0.11515e2, -0.53984e1, -0.37937e2],
        "exp": [0.355, 2.45, 1.0, 5.1]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2010)",
               "__doc__": """Perkins, R.A. and Huber, M.L., unpublished work, 2010.""",

               "Tref": 367.85, "kref": 1,
               "no": [-0.237681e-2, 0.781106e-2, 0.147206e-1],
               "co": [0, 1, 2],

               "Trefb": 367.85, "rhorefb": 4.17, "krefb": 1.,
               "nb": [0.223071e-2, -0.278314e-1, 0.149657e-1, 0.0,
                      -0.252976e-3, 0.0, 0.188873e-1, -0.608795e-2, 0.0, 0.0],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.285e-10, "Tcref": 551.775}

#    _thermal = thermo0,


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=R1234yf(T=300., rho=5.0)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)


