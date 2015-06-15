#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R1234ze(MEoS):
    """Multiparameter equation of state for R1234ze"""
    name = "trans-1,3,3,3-tetrafluoropropene"
    CASNumber = "29118-24-9"
    formula = "CHF=CHCF3"
    synonym = "R-1234ze"
    rhoc = unidades.Density(489.238464)
    Tc = unidades.Temperature(382.513)
    Pc = unidades.Pressure(3634.9, "kPa")
    M = 114.0416  # g/mol
    Tt = unidades.Temperature(168.62)
    Tb = unidades.Temperature(254.177)
    f_acent = 0.313
    momentoDipolar = unidades.DipoleMoment(1.27, "Debye")
    id = 671

    CP1 = {"ao": 4.0,
           "an": [], "pow": [],
           "ao_exp": [9.3575, 10.717], 
           "exp": [513, 1972],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 6.259,
           "an": [], "pow": [],
           "ao_exp": [7.303, 8.597, 2.333], 
           "exp": [691, 1705, 4216],
           "ao_hyp": [], "hyp": []}
           
    CP3 = {"ao": 5.8887,
           "an": [], "pow": [],
           "ao_exp": [7.0804, 9.3371, 2.5577],
           "exp": [620, 1570, 3953],
           "ao_hyp": [], "hyp": []}

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-10.8724711, -30.1326538],
           "ao_exp": [6.07536, 9.95795],
           "titao": [289/Tc, 1303/Tc], 
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze of Thol and Lemmon (2013).",
        "__doi__": {"autor": "Thol, M. and Lemmon, E.W.",
                    "title": "to be published in Int. J. Thermophys., 2013.", 
                    "ref": "",
                    "doi": ""}, 
                    
        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR", 
        
        "Tmin": Tt, "Tmax": 420.0, "Pmax": 20000.0, "rhomax": 13.26, 
        "Pmin": 0.2187, "rhomin": 13.26, 

        "nr1": [0.03982797, 1.812227, -2.537512, -0.5333254, 0.1677031],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.223, 0.755, 1.24, 0.44],

        "nr2": [-1.323801, -0.6694654, 0.8072718, -0.7740229, -0.01843846],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2., 2.2, 1.2, 1.5, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.407916, -0.4237082, -0.2270068, -0.805213, 0.00994318, -0.008798793],
        "d3": [1, 1, 3, 3, 2, 1],
        "t3": [1.33, 1.75, 2.11, 1.0, 1.5, 1.0],
        "alfa3": [1.0, 1.61, 1.24, 9.34, 5.78, 3.08],
        "beta3": [1.21, 1.37, 0.98, 171, 47.4, 15.4],
        "gamma3": [0.943, 0.642, 0.59, 1.2, 1.33, 0.64],
        "epsilon3": [0.728, 0.87, 0.855, 0.79, 1.3, 0.71],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze of McLinden et al. (2010).",
        "__doi__": {"autor": "McLinden, M.O., Thol, M., and Lemmon, E.W.",
                    "title": "Thermodynamic Properties of trans-1,3,3,3-Tetrafluoropropene [R1234ze(E)]: Measurements of Density and Vapor Pressure and a Comprehensive Equation of State", 
                    "ref": "International Refrigeration and Air Conditioning Conference at Purdue, July 12-15, 2010.",
                    "doi": "10.0000_docs.lib.purdue.edu_generic-99DA7EA2C877"}, 
                    
        "R": 8.314472,
        "cp": CP2,
        "ref": "IIR", 
        
        "Tmin": Tt, "Tmax": 420.0, "Pmax": 20000.0, "rhomax": 13.20, 
        "Pmin": 0.23, "rhomin": 13.19, 

        "nr1": [0.055563, 1.66927, -2.53408, -0.475075, 0.190055],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.34, 0.91, 1.23, 0.46],

        "nr2": [-1.25154, -0.742195, 0.537902, -0.741246, -0.0355064],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.26, 2.50, 2.0, 2.24, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.58506, -0.502086, -0.19136, -0.975576],
        "d3": [1, 1, 3, 3],
        "t3": [1.06, 1.79, 3.75, 0.92],
        "alfa3": [1.02, 1.34, 1.08, 6.41],
        "beta3": [1.19, 2.29, 1.15, 131.8],
        "gamma3": [1.14, 0.667, 0.505, 1.22],
        "epsilon3": [0.711, 0.914, 0.694, 0.731],
        "nr4": []}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze of McLinden et al. (2010).",
        "__doi__": {"autor": "McLinden, M.O., Thol, M., and Lemmon, E.W.",
                    "title": "Thermodynamic Properties of trans-1,3,3,3-Tetrafluoropropene [R1234ze(E)]: Measurements of Density and Vapor Pressure and a Comprehensive Equation of State", 
                    "ref": "unpublished equation, similar to helmholtz2",
                    "doi": ""}, 
        "R": 8.314472,
        "cp": CP3,
        "ref": "IIR", 
        
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

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234yf of Akasaka (2011).",
        "__doi__": {"autor": "Akasaka, R.",
                    "title": "New Fundamental Equations of State with a Common Functional Form for 2,3,3,3-Tetrafluoropropene (R-1234yf) and trans-1,3,3,3-Tetrafluoropropene (R-1234ze(E))", 
                    "ref": "Int J Thermophys (2011) 32:1125–1147",
                    "doi": "10.1007/s10765-011-0992-0"}, 
        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR", 

        "Tmin": 240., "Tmax": 420.0, "Pmax": 15000.0, "rhomax": 13.20, 
        "Pmin": 0.23, "rhomin": 13.19, 

        "nr1": [0.85579765e1, -0.94701332e1, -0.25013623, 0.13789870, 0.12177113e-1],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.66886, 0.83392, 1.6982, 1.8030, 0.36657],
 
        "nr2": [-0.14227996, 0.10096648, 0.17504319e-1, -0.17627303e-1, -0.14705120e-1, 0.37202269, -0.30138266, -0.92927274e-1, 0.87051177e-1, 0.18113770e-1, -0.16018424e-1, 0.53809860e-2],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.8666, 1.0194, 0, 1.1655, 8.3101, 6.1459, 8.3495, 6.0422,
               7.444, 15.433, 21.543, 15.499],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    eq = helmholtz1, helmholtz2, helmholtz3, helmholtz4

    _surface = {"__doi__": 
                    {"autor": "Tanaka, K., Higashi, Y.",
                     "title": "Surface Tension of trans-1,3,3,3-Tetrafluoropropene and trans-1,3,3,3-Tetrafluoropropene + Difluoromethane Mixture", 
                     "ref": "J. Chem. Eng. Japan, 2013",
                     "doi": "10.1252/jcej.13we021"}, 
                "sigma": [0.05681], "exp": [1.23]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.5888, 1.9696, -2.0827, -4.1238],
        "exp": [1.0, 1.5, 2.2, 4.6]}
    _liquid_Density = {
        "eq": 1,
        "ao": [1.1913, 2.2456, -1.7747, 1.3096],
        "exp": [0.27, 0.7, 1.25, 1.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.0308, -5.0422, -11.5, -37.499, -77.945],
        "exp": [0.24, 0.72, 2.1, 4.8, 9.5]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2011)",
               "__doi__": {"autor": "Perkins, R.A. and Huber, M.L.",
                            "title": "Measurement and Correlation of the Thermal Conductivity of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf) and trans-1,3,3,3-Tetrafluoropropene (R1234ze(E))", 
                            "ref": "J. Chem. Eng. Data, 2011, 56 (12), pp 4868–4874",
                            "doi": "10.1021/je200811n"}, 
               "__test__": """
                    >>> st=R1234ze(T=250, P=5e4, eq=1)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    2.80451 0.0098503
                    >>> st=R1234ze(T=300, P=1e5, eq=1)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    4.67948 0.013933
                    >>> st=R1234ze(T=250, P=2e7, eq=1)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    1349.37 0.100066
                    >>> st=R1234ze(T=300, P=2e7, eq=1)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    1233.82 0.085389
                    """, # Table 2, Pag 4872

               "Tref": 382.52, "kref": 1,
               "no": [-0.0103589, 0.0308929, 0.000230348],
               "co": [0, 1, 2],

               "Trefb": 382.52, "rhorefb": 4.29, "krefb": 1.,
               "nb": [-0.428296e-1, 0.927099e-1, -0.702107e-1, 0.249708e-1,
                      -0.301838e-2, 0.434288e-1, -0.605844e-1, 0.440187e-1,
                      -0.155082e-1, 0.210190e-2],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.835e-10, "Tcref": 573.78}

    _thermal = thermo0,
