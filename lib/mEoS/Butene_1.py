#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Butene_1(MEoS):
    """Multiparameter equation of state for 1-butene"""
    name = "butene"
    CASNumber = "106-98-9"
    formula = "CH3-CH2-CH=CH2"
    synonym = ""
    rhoc = unidades.Density(237.8907968)
    Tc = unidades.Temperature(419.29)
    Pc = unidades.Pressure(4005.1, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(87.8)
    Tb = unidades.Temperature(266.84)
    f_acent = 0.192
    momentoDipolar = unidades.DipoleMoment(0.339, "Debye")
    id = 24
           
    Fi1 = {"ao_log": [1, 2.9197],
           "pow": [0, 1],
           "ao_pow": [-0.00101126, 2.3869174],
           "ao_exp": [2.9406, 6.5395, 14.5395, 5.8971],
           "titao": [274/Tc, 951/Tc, 2127/Tc, 5752/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for 1-butene of Lemmon and Ihmels (2005)",
        "__doi__": {"autor": "Lemmon, E.W., Ihmels, E.C.",
                    "title": "Thermodynamic properties of the butenes: Part II. Short fundamental equations of state", 
                    "ref": "Fluid Phase Equilibria 228 – 229 (2004), 173 – 187.",
                    "doi":  "10.1016/j.fluid.2004.09.004"}, 
        "__test__": """
            >>> st=Butene_1(T=350, rho=0)
            >>> print "%0.0f %0.1f %0.1f %0.5g %0.5g %0.5g %0.5g" % (st.T, st.rhoM, st.P.MPa, st.hM.kJkmol, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            350 0.0 0.0 29617 88.208 96.522 238.24
            >>> st=Butene_1(T=350, rho=0.3*56.10632)
            >>> print "%0.0f %0.1f %.5g %.5g %0.5g %0.5g %0.5g %0.5g" % (st.T, st.rhoM, st.P.MPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            350 0.3 0.75679 28321 87.626 92.719 108.45 211.38
            >>> st=Butene_1(T=350, rho=10*56.10632)
            >>> print "%0.0f %0.1f %.5g %.5g %0.5g %0.5g %0.5g %0.5g" % (st.T, st.rhoM, st.P.MPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            350 10.0 17.864 11377 31.563 97.760 135.36 843.31
            >>> st=Butene_1(T=440, rho=4*56.10632)
            >>> print "%0.0f %0.1f %.5g %.5g %0.5g %0.5g %0.5g %0.5g" % (st.T, st.rhoM, st.P.MPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            440 4.0 5.3245 29454 80.191 124.13 416.03 151.49
            """, # Table 9, Pag 186
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 525., "Pmax": 70000.0, "rhomax": 14.59, 
        "Pmin": 0.0000000008, "rhomin": 14.58, 

        "nr1": [0.78084, -2.8258, 0.99403, 0.017951, 0.088889, 0.00024673],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.22846, -0.074009, -0.22913, -0.062334, -0.025385, 0.011040],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.05644], "exp": [1.248]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.71727e1, 0.26360e1, -0.20781e1, -0.28860e1, -0.13041e1],
        "exp": [1, 1.5, 2, 4.35, 16.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.16857e2, -0.46280e2, 0.53727e2, -0.23314e2, 0.18889e1],
        "exp": [0.547, 0.73, 0.92, 1.14, 2.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31106e1, -0.63103e1, -0.19272e2, -0.48739e2, -0.99898e2, -0.19001e3],
        "exp": [0.415, 1.27, 3.34, 7.0, 14.5, 28.0]}
