#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class DMC(MEoS):
    """Multiparameter equation of state for dimethyl carbonate

    >>> metilciclohexano=DMC(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "dimethyl carbonate"
    CASNumber = "616-38-6"
    formula = "C3H6O3"
    synonym = ""
    rhoc = unidades.Density(358.050803706)
    Tc = unidades.Temperature(557.376)
    Pc = unidades.Pressure(4835.08, "kPa")
    M = 90.07794  # g/mol
    Tt = unidades.Temperature(277.06)
    Tb = unidades.Temperature(363.112)
    f_acent = 0.33327
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 184
    # id=1798

    CP1 = {"ao": 8.6169,
           "an": [],
           "pow": [],
           "ao_exp": [0.69884, 0.13132e2, 0.69241, 0.83174e1],
           "exp": [1150., 1339., 1590., 3111.],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for DMC of Zhou et al. (2010).",
        "__doc__":  u"""Zhou, Y., Wu, J., and Lemmon, E.W. "Equations for the Thermophysical Properties of Dimethyl Carbonate," submitted to J. Phys. Chem. Ref. Data, 2010.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [0.376822e-1, 0.128917e1, -0.254658e1, -0.208420, 0.187066],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.2147, 1.033, 0.9774, 0.6649],

        "nr2": [-0.212435, 0.275468, -0.337738e-1, -0.353955, -0.209746, -0.504864],
        "d2": [1, 2, 7, 1, 2, 3],
        "t2": [0.5999, 1.095, 1.088, 5.117, 6.633, 2.465],
        "c2": [1, 1, 1, 2, 2, 2],
        "gamma2": [1]*6,

        "nr3": [-0.111732, 0.140922e1, -0.402018, -0.756710e-1, -0.875139e-1],
        "d3": [1, 1, 1, 3, 3],
        "t3": [0.9935, 1.805, 3.0, 2.235, 2.948],
        "alfa3": [1.609, 1.014, 1.867, 1.307, 1.798],
        "beta3": [9.306, 0.3138, 0.6848, 0.4879, 0.8123],
        "gamma3": [0.6345, 1.87, 1.172, 1.55, 0.6256],
        "epsilon3": [0.8785, 1.051, 1.049, 1.351, 1.179],
        "nr4": []}

    eq = helmholtz1,

    #FIXME: Dan malas estimaciones que provocan bucle
#    _vapor_Pressure={ "eq": 5, "ao": [-0.75913e1, -0.14610e1, 0.85506e1, -0.90296e1, -0.35161e1], "exp": [1.0, 1.5, 1.9, 2.4, 6.0]}
#    _liquid_Density={ "eq": 1, "ao": [0.14791, 0.12391e2, -0.25090e2, 0.22517e2, -0.70191e1], "exp": [0.21, 0.7, 1.0, 1.3, 1.6]}
#    _vapor_Density={ "eq": 3, "ao": [-0.23324e1, -0.50936e1, -0.20195e2, -0.94319e2, 0.25801e3, -0.49161e3], "exp": [0.461, 0.81, 2.8, 6.8, 10.0, 12.0]}

    visco0 = {"eq": 1, "omega": 3,
              "__name__": "Zhou (2010)",
              "__doc__":  u"""Zhou, Y., Wu, J., and Lemmon, E.W. "Equations for the Thermophysical Properties of Dimethyl Carbonate," submitted to J. Phys. Chem. Ref. Data, 2010.""",
              "ek": 442.3, "sigma": 0.510747,
              "Tref": 557.376, "rhoref": 3.9749*M,
              "n_chapman": 0.20555,
              "n_ideal": [],
              "t_ideal": [],

              "Tref_res": 557.376, "rhoref_res": 3.9749*M, "etaref_res": 1,
              "n_poly": [5.07808, -0.056734, 0.00832177, 35.459838, 0.0513528],
              "t_poly": [-0.1, -3.0968, -2.8945, 0.0731, -3.9871],
              "d_poly": [4, 10, 12, 2, 0],
              "g_poly": [0, 0, 0, 0, 0, 0],
              "c_poly": [0, 1, 1, 2, 3],
#             "n_num": [],
#             "t_num": [],
#             "d_num": [],
#             "g_num": [],
#             "c_num": [],
#             "n_den": [],
#             "t_den": [],
#             "d_den": [],
#             "g_den": [],
#             "c_den": []
}

    _viscosity = visco0,


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=DMC(T=400., rho=5.)
    print "%0.1f %0.5f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
    print cyc5.k.mWmK, cyc5.mu.muPas
