#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoSBlend
from lib import unidades


class R410a(MEoSBlend):
    """Multiparameter equation of state for R410A (50% R32, 50% R125)"""
    name = "R410A"
    CASNumber = ""
    formula = "R32+R125"
    synonym = "R410A"
    rhoc = unidades.Density(459.0300696)
    Tc = unidades.Temperature(344.494)
    Pc = unidades.Pressure(4901.2, "kPa")
    M = 72.5854  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(221.71)
    f_acent = 0.296
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 62
    # id = None

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.1],
           "ao_pow": [36.8871, 7.15807, -46.87575],
           "ao_exp": [2.0623, 5.9751, 1.5612],
           "titao": [697/Tc, 1723/Tc, 3875/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-410A of Lemmon (2003)",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "Pseudo-Pure Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C", 
                    "ref": "Int. J. Thermophys., 24(4):991-1006, 2003.",
                    "doi": "10.1023/A:1025048800563"}, 
        "__test__": """
            >>> st=R410a(T=300, rhom=0)
            >>> print "%0.3g %0.1f %0.1f %0.3f %0.3f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 0.0 0.0 50.400 58.714 200.08
            >>> st=R410a(T=300, P=R410a._bubbleP(300))
            >>> print "%0.3g %0.4f %0.5f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 1.7404 14.45917 67.147 125.50 418.60
            >>> st=R410a(T=300, P=R410a._dewP(300))
            >>> print "%0.3g %0.4f %0.5f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 1.7351 0.95997 67.335 107.60 160.94
            >>> st=R410a(T=250, rhom=18)
            >>> print "%0.3g %0.3f %0.1f %0.3f %0.3f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            250 17.651 18.0 62.521 98.401 800.83
            """, # Table V, Pag 998

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR", 
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 19.51, 
        "Pmin": 29.16, "rhomin": 19.51, 

        "Tj": 344.494, "Pj": 4.9012, 
        "dew": {"i": [1*2, 1.6*2, 2.4*2, 5*2], 
                "n": [-7.4411, 1.9883, -2.4925, -3.2633]}, 
        "bubble": {"i": [1*2, 1.8*2, 2.4*2, 4.9*2], 
                   "n": [-7.2818, 2.5093, -3.2695, -2.8022]}, 

        "nr1": [0.987252, -0.103017e1, 0.117666e1, -0.138991, 0.302373e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.44, 1.2, 2.97, 2.95, 0.2],

        "nr2": [-0.253639e1, -0.196680e1, -0.830480, 0.172477, -0.261116,
                -0.745473e-1, 0.679757, -0.652431, 0.553849e-1, -0.710970e-1,
                -0.875332e-3, 0.200760e-1, -0.139761e-1, -0.185110e-1,
                0.171939e-1, -0.482049e-2],
        "d2": [1, 2, 3, 5, 5, 5, 1, 1, 4, 4, 9, 2, 2, 4, 5, 6],
        "t2": [1.93, 1.78, 3., 0.2, 0.74, 3, 2.1, 4.3, 0.25, 7, 4.7, 13, 16,
               25, 17, 7.4],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3],
        "gamma2": [1]*16}

    eq = helmholtz1,

    _surface = {"sigma": [0.06443], "exp": [1.245]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.4411, 1.9883, -2.4924, -3.2633],
        "exp": [1, 1.6, 2.4, 5]}
    _liquid_Pressure = {
        "eq": 5,
        "ao": [-7.2818, 2.5093, -3.2695, -2.8022],
        "exp": [1, 1.8, 2.4, 4.9]}
