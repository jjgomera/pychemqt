#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class SO2(MEoS):
    """Multiparameter equation of state for sulfur dioxide"""
    name = "sulfur dioxide"
    CASNumber = "7446-09-5"
    formula = "SO2"
    synonym = "R-764"
    rhoc = unidades.Density(525.)
    Tc = unidades.Temperature(430.64)
    Pc = unidades.Pressure(7884.0, "kPa")
    M = 64.0638  # g/mol
    Tt = unidades.Temperature(197.7)
    Tb = unidades.Temperature(263.13)
    f_acent = 0.2557
    momentoDipolar = unidades.DipoleMoment(1.6, "Debye")
    id = 51

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1],
           "ao_pow": [-4.5328346436, 4.4777967379, -0.01560057996],
           "ao_exp": [1.062, 1.9401],
           "titao": [775/Tc, 1851/Tc]}

    CP1 = {"ao": 4,
           "an": [0.72453e-4], "pow": [1],
           "ao_exp": [1.062, 1.9401], "exp": [775, 1851],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 0.4021066/8.3143*64.066,
           "an": [0.87348570e-3/8.3143*64.066, -0.45968820e-6/8.3143*64.066,
                  -0.13328400e-11/8.3143*64.066, 0.23785000e-13/8.3143*64.066],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for sulfur dioxide of Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=SO2(T=432, rho=8*64.0638)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            432 8 8052.256 20821.200 56.819 61.478 4877.456 171.538
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 525.0, "Pmax": 35000.0, "rhomax": 25.30, 
        "Pmin": 1.66, "rhomin": 25.29, 

        "nr1": [0.93061, -1.9528, -0.17467, 0.061524, 0.00017711],
        "d1": [1, 1, 1, 3, 7, ],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.21615, 0.51353, 0.010419, -0.25286, -0.054720, -0.059856,
                -0.016523],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur dioxide of Polt (1987).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe", 
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""}, 
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP", 
        
        "Tmin": 273.0, "Tmax": 523.0, "Pmax": 32000.0, "rhomax": 22.91, 
        "Pmin": 11.82, "rhomin": 23.0, 

        "nr1": [0.789407019882, -0.170449580056e1, 0.115984637964e1,
                -0.576307837294, 0.249237283833e1, -0.518115678632e1,
                0.320766081899e1, -0.123636065893e1, 0.144419600938e-1,
                -0.15380705504, 0.386324300525, 0.292550313202, -0.372445361392,
                -0.636924333910e-1, 0.986166451596e-1, -0.216993783055e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.789407019882, 0.170449580056e1, -0.115984637964e1,
                -0.480876182378, 0.164910076886e1, -0.133861069604e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1]*6}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.1016572, -0.020501, -0.0057962],
                "exp": [1.25, 2.25, 3.25]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73845e1, 0.22867e1, -0.24669e1, -0.32217e1, 0.23109],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.17156e2, -0.60441e2, 0.81407e2, -0.51871e2, 0.16754e2],
        "exp": [0.57, 0.8, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.33832e1, -0.76873e1, -0.23614e2, -0.13720e3, 0.18664e4, -0.24469e4],
        "exp": [0.424, 1.4, 3.6, 8.5, 13.0, 14.0]}
