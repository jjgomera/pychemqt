#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class iC5(MEoS):
    """Multiparameter equation of state for isopentane

#    >>> isopentano=iC5(T=300, P=0.1)
#    >>> print "%0.1f %0.2f %0.4f %0.4f %0.7f %0.4f %0.4f %0.2f" % (isopentano.T, isopentano.rho, isopentano.u.kJkg, isopentano.h.kJkg, isopentano.s.kJkgK, isopentano.cv.kJkgK, isopentano.cp.kJkgK, isopentano.w)
#    300.0 613.09 -2.3963 -2.2332 -0.0074246 1.7111 2.2836 951.69
#
    """
    name = "isopentane"
    CASNumber = "78-78-4"
    formula = "(CH3)2-CH-CH2-CH3"
    synonym = "R-601a"
    rhoc = unidades.Density(236.)
    Tc = unidades.Temperature(460.35)
    Pc = unidades.Pressure(3378.0, "kPa")
    M = 72.14878  # g/mol
    Tt = unidades.Temperature(112.65)
    Tb = unidades.Temperature(300.98)
    f_acent = 0.2274
    momentoDipolar = unidades.DipoleMoment(0.11, "Debye")
    id = 7

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [2.5822330405, 1.1609103419],
           "ao_exp": [7.4056, 9.5772, 15.765, 12.119],
           "titao": [442/Tc, 1109/Tc, 2069/Tc, 4193/Tc]}

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [7.4056, 9.5772, 15.765, 12.119],
           "exp": [442, 1109, 2069, 4193],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [11.7618, 20.1101, 33.1688, 0],
           "hyp": [0.635392636*Tc, 1.977271641*Tc, 4.169371131*Tc, 0]}

    CP3 = {"ao": 0.396504/8.3143*72.151,
           "an": [0.260678e-2/8.3143*72.151, 0.93677e-5/8.3143*72.151,
                  -0.158286e-7/8.3143*72.151,  0.76525e-11/8.3143*72.151],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 21.3861/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [2.1524504e8/8.3159524*4.184, 2.8330244e7/8.3159524*4.184, 0, 0],
           "hyp": [1.70158e3, 7.75899e2, 0, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isopentane of Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=iC5(T=462, rho=3*72.14878)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            462 3 3458.617 37318.534 94.634 190.631 4660.943 96.324
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 1000000.0, "rhomax": 13.3, 
        "Pmin": 0.83e-7, "rhomin": 10.925, 

        "nr1": [1.0963, -3.0402, 1.0317, -0.15410, 0.11535, 0.00029809],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.39571, -0.045881, -0.35804, -0.10107, -0.035484, 0.018156],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 1000000.0, "rhomax": 13.3, 
        "Pmin": 0.83e-7, "rhomin": 10.925, 

        "nr1": [0.11017531966644e1, -0.30082368531980e1, 0.99411904271336,
                -0.14008636562629, 0.11193995351286, 0.29548042541230e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.36370108598133, -0.48236083488293e-1, -0.35100280270615,
                -0.10185043812047, -0.35242601785454e-1, 0.19756797599888e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP3,

        "Tmin": 200.0, "Tmax": 553.0, "Pmax": 7500.0, "rhomax": 5.2252, 
        "Pmin": 51.964, "rhomin": 8.7248, 

        "nr1": [-0.143819012123e1, 0.138298276836e1, -0.203328695121,
                0.619304204378, -0.311353942178e1, 0.316914412369e1,
                -0.218812895934e1, 0.211230723299, 0.765790344231,
                -0.851773312153, 0.706192861166, -0.165802139239,
                0.781356542750e-1, 0.106516957202, -0.205642736936,
                0.360787537633e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.143819012123e1, -0.138298276836e1, 0.203328695121,
                -0.213463476736e1, 0.547491842897e1, -0.335666356499e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.002528]*6}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Starling (1973)",
        "__doc__":  u"""Starling, K.E., "Fluid Thermodynamic Properties for Light Petroleum Systems," Gulf Publishing Company, 1973.""",
        "R": 8.3159524,
        "cp": CP4,

        "Tmin": 199.82, "Tmax": 589.0, "Pmax": 55000.0, "rhomax": 9.9258626, 
        "Pmin": 0.34375, "rhomin": 9.9259, 

        "nr1": [0.179378842786e1, 0.258488286720, -0.812072482201,
                -0.753941018871, 0.565338153509e-1, -0.115706201242e-2,
                0.406090628523, -0.469700474204, -0.967480812300e-1,
                0.958936263943e-2, 0.197520012548e-2],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.179378842786e1, -0.431019031876],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.48056842]*2}

    eq = helmholtz1, #GERG, helmholtz3, helmholtz4

    _surface = {"sigma": [0.05106], "exp": [1.21]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.26977],  "expt0": [-1.], "expd0": [1.],
                   "a1": [25.31, 0.025], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [108.9, 63.68, -15447, -5449.3],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.83e-7,
                "Tmin": Tt, "Tmax": 2000.0,
                "a1": [-7127700000000, 7127700000001], "exp1": [0, 1.563],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.72392e1, 0.22635e1, -0.18237e1, -0.29997e1, -0.27752e1],
        "exp": [1., 1.5, 2.02, 4.24, 16.1]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18367e2, -0.30283e2, 0.13557e2, -0.90533, 0.20927e1],
        "exp": [1.21, 1.41, 1.65, 0.09, 0.164]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.38825e2, 0.79040e2, -0.48791e2, -0.21603e2, -0.57218e2, -0.15164e3],
        "exp": [0.565, 0.66, 0.77, 3.25, 7.3, 16.6]}

    visco0 = {"eq": 2, "omega": 3,
              "__name__": "NIST",
              "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
              "ek": 341.06, "sigma": 0.56232,
              "n_chapman": 0.2267237/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-4.57981980159405, -3393.5243856, 9.3806654324, 33641.3512,
                    0.15624235969, 122.90017543, -20914.795166],
              "rhoc": 3.24}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "NIST14",
               "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

               "Tref": 341.06, "kref": 1e-3,
               "no": [1.35558587, -0.152666315743857, 1],
               "co": [0, -1, -96],

               "Trefb": 460.51, "rhorefb": 3.24, "krefb": 1e-3,
               "nb": [18.6089331038, -5.836570612990, 3.489871005290,
                      0.704467355508, -0.206501417728, -0.223070394020],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.9316e-9, "Tcref": 690.525}

    _thermal = thermo0,
