#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class C1Cyclohexane(MEoS):
    """Multiparameter equation of state for methylcyclohexane

    >>> metilciclohexano=C1Cyclohexane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "methylcyclohexane"
    CASNumber = "108-87-2"
    formula = "C6H11-CH3"
    synonym = ""
    rhoc = unidades.Density(267.0660832)
    Tc = unidades.Temperature(572.2)
    Pc = unidades.Pressure(3470.0, "kPa")
    M = 98.18606  # g/mol
    Tt = unidades.Temperature(146.7)
    Tb = unidades.Temperature(374.)
    f_acent = 0.23
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 39

    CP1 = {"ao": 2.04122,
           "an": [0.016417, 0.000185315, -3.14826e-7, 1.65567e-10],
           "pow": [1, 2, 3, 4],
           "ao_exp": [],
           "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for methylcyclohexane of Lemmon (2007).",
        "__doc__":  """""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.3026, -2.6270, 0.68834, -0.16415, 0.092174, 0.0003842],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.38, 1.2, 2.14, 1.6, 0.3, 0.7],

        "nr2": [-0.29737, -0.078187, -0.049139, -0.30402, -0.074888],
        "d2": [1, 2, 5, 1, 4],
        "t2": [2.7, 3.25, 2.35, 3.7, 4.1],
        "c2": [1, 1, 1, 2, 2],
        "gamma2": [1]*17}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.65871e1, -0.56553e1, 0.68947e1, -0.41281e1, -0.25444e1],
        "exp": [1.0, 1.5, 1.6, 3.2, 10.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18273e-1, 0.15215e2, -0.21951e2, 0.94466e1, 0.16781],
        "exp": [0.1, 0.64, 0.8, 1.0, 4.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.52572e1, -0.13417e2, -0.24271e1, -0.54482e2, -0.15791e3],
        "exp": [0.544, 2.3, 2.5, 6.1, 15.0]}

    visco0 = {"eq": 4, "omega": 1,
              "__doc__": """regression as in S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "Tref": 572.2, "muref": 1.0,
              "ek": 454.3, "sigma": 0.5801,
              "n_ideal": [32.8082, -104.308, 98.4289, -13.7085],
              "t_ideal": [0, 0.25, 0.5, 0.75],

              "a": [-0.464134e-5, 0.0, 0.397245e-6],
              "b": [-0.381691e-4, 0.866218e-4, 0.414300e-6],
              "c": [0.389947e-3, -0.194159e-3, 0.0],
              "A": [-0.297679e-7, 0.223799e-9, 0.0],
              "B": [0.384063e-8, 0.0, 0.0],
              "C": [0.0, 0.0, 0.0],
              "D": [0.0, 0.0, 0.0]}

    visco1 = {"eq": 5, "omega": 3,
              "__doc__": """T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E. "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties" Ind. Eng. Chem. Res. 1998, 27, 671-679""",
              "__name__": "Chung (1988)",
              "w": 0.885, "mur": 0.0, "k": 0.0}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2008)",
               "__doc__": """Perkins, R.A. Hammerschmidt, U. and Huber, M.L., "Measurement and Correlation of the Thermal Conductivity of Methylcyclohexane and Propylcyclohexane from 300 K to 600 K at Pressures to 60 MPa," J. Chem. Eng. Data, 53:2120-2127(2008).""",

               "Tref": 572.2, "kref": 1,
               "no": [0.289968e-2, -0.180666e-1, 0.727576e-1, -0.129778e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 572.2, "rhorefb": 2.72, "krefb": 1.,
               "nb": [0.91914900e-1, -0.79040800e-1, -0.81708800e-1,
                      0.92391100e-1, 0.29644900e-1, -0.42849800e-1,
                      -0.29983400e-2, 0.72786000e-2, 0.0, 0.0],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 6.24e-10, "Tcref": 858.3}

    _thermal = thermo0,


if __name__ == "__main__":
    # import doctest
    # doctest.testmod()

    cyc5 = C1Cyclohexane(T=300., P=0.1, visco=1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
