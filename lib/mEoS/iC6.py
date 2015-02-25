#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class iC6(MEoS):
    """Multiparameter equation of state for isohexane

    >>> isohexano=iC6(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.4f %0.4f %0.7f %0.4f %0.4f %0.2f" % (isohexano.T, isohexano.rho, isohexano.u.kJkg, isohexano.h.kJkg, isohexano.s.kJkgK, isohexano.cv.kJkgK, isohexano.cp.kJkgK, isohexano.w)
    300.0 646.90 -77.7760 -77.6215 -0.2451750 1.7371 2.2442 1037.16
    """
    name = "isohexane"
    CASNumber = "107-83-5"
    formula = "(CH3)2-CH-(CH2)2-CH3"
    synonym = ""
    rhoc = unidades.Density(233.966)
    Tc = unidades.Temperature(497.7)
    Pc = unidades.Pressure(3040.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(119.6)
    Tb = unidades.Temperature(333.36)
    f_acent = 0.2797
    momentoDipolar = None
    id = 52

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [7.9127, 16.871, 19.257, 14.075],
           "exp": [325, 1150, 2397, 5893],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isohexane of Lemmon and Span (2006)",
        "__doc__":  u"""Lemmon, E.W., Span, R. Short fundamental equations of state for 20 industrial fluids. J. Chem. Eng. Data 51 (2006), 785 – 850.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 1000000.0, "rhomax": 9.38, 
        "Pmin": 7.34e-9, "rhomin": 9.37, 

        "nr1": [1.1027, -2.9699, 1.0295, -0.21238, 0.11897, 0.00027738],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40103, -0.034238, -0.43584, -0.11693, -0.019262, 0.0080783],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.053], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74130e1, 0.16267e1, -0.22311e1, -0.26040e1, -0.29490e1],
        "exp": [1.0, 1.5, 2.62, 4.56, 16.3]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18489e2, -0.43541e2, 0.43985e2, -0.16581e2, 0.64563],
        "exp": [0.59, 0.77, 0.96, 1.15, 3.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.41180e1, -0.61956e1, -0.21190e2, -0.58972e2, -0.15824e3],
        "exp": [0.4824, 1.418, 3.32, 7.1, 16.1]}

    visco0 = {"eq": 2, "omega": 3,
              "__name__": "NIST",
              "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
              "ek": 368.52, "sigma": 0.61222,
              "n_chapman": 0.2267237/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-13.294469653994, -466.41004563, 15.438316998,
                    -3363.2028894, -0.11398677788, 171.32077134, 2849.7100897],
              "rhoc": 2.727}

    _viscosity = visco0,

    thermo0 = {"eq": 1, "critical": 0,
               "__name__": "NIST14",
               "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

               "Tref": 368.52, "kref": 1e-3,
               "no": [1.35558587, -0.152808259573429, 1],
               "co": [0, -1, -96],

               "Trefb": 498.05, "rhorefb": 2.727, "krefb": 1e-3,
               "nb": [13.747515904, 10.1607102792, -7.75232868497,
                      0.627943006907, 1.9518640415, -0.293574041046],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6}

    _thermal = thermo0,
