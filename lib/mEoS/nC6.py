#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class nC6(MEoS):
    """Multiparameter equation of state for n-hexane

    >>> hexano=nC6(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (hexano.T, hexano.rho, hexano.u.kJkg, hexano.h.kJkg, hexano.s.kJkgK, hexano.cv.kJkgK, hexano.cp.kJkgK, hexano.w)
    300.0 653.08 14.32 14.47 0.69082 1.7504 2.2596 1048.3
    """
    name = "hexane"
    CASNumber = "110-54-3"
    formula = "CH3-(CH2)4-CH3"
    synonym = ""
    rhoc = unidades.Density(233.1819)
    Tc = unidades.Temperature(507.82)
    Pc = unidades.Pressure(3034.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(177.83)
    Tb = unidades.Temperature(341.86)
    f_acent = 0.299
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 10
    _Tr = unidades.Temperature(487.762087)
    _rhor = unidades.Density(235.700888)
    _w = 0.298052404

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [14.345969349, -96.165722367],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [11.6977, 26.8142, 38.6164, 0],
           "hyp": [0.359036667, 1.691951873, 3.596924107, 0]}

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.3888640e6, 0.1979523e8, 0.1288410e9, 0],
           "hyp": [0.182326e3, 0.8592070e3, 0.1826590e4, 0]}

    CP2 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [11.6977, 26.8142, 38.6164, 0],
           "hyp": [0.359036667*Tc, 1.691951873*Tc, 3.596924107*Tc, 0]}

    CP3 = {"ao": 2.5200507,
           "an": [0.52806530e-1, -0.57861557e-5, -0.10899040e-7, -0.18988742e-12],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 26.6225/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [2.3738446e8/8.3159524*4.184, 3.5806766e7/8.3159524*4.184, 0, 0],
           "hyp": [1.71849e3, 8.02069e2, 0, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for hexane of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.", 
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"}, 
        "__test__": """
            >>> st=nC6(T=700, rho=200)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            3.1802 10.221 3.6535
            >>> st2=nC6(T=750, rho=100)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            213.09 0.33219
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 8.85, 
        "Pmin": 0.001277, "rhomin": 8.8394, 

        "nr1": [0.10553238013661e1, -0.26120615890629e1, 0.76613882967260,
                -0.29770320622459, 0.11879907733358, 0.27922861062617e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.46347589844105, 0.11433196980297e-1, -0.48256968738131,
                -0.93750558924659e-1, -0.67273247155994e-2, -0.51141583585428e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 8.85, 
#        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.10553238013661e1, -0.26120615890629e1, 0.76613882967260,
                -0.29770320622459, 0.11879907733358, 0.27922861062617e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.46347589844105, 0.11433196980297e-1, -0.48256968738131,
                -0.93750558924659e-1, -0.67273247155994e-2, -0.51141583585428e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP3,

        "Tmin": 223.0, "Tmax": 623.0, "Pmax": 510000.0, "rhomax": 8.726125, 
        "Pmin": 0.001277, "rhomin": 8.8394, 

        "nr1": [-0.157654494847e1, 0.178731485778e1, -0.341262936801,
                0.114919468260e1, -0.381451065649e1, 0.356688884337e1,
                -0.274863278063e1, 0.391987699726, 0.346062554746,
                -0.139140552239, 0.489013943543, -0.529751545354e-1,
                -0.149303737787, 0.455990262306e-1, -0.564866336099e-1,
                0.152437539639e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.157654494847e1, -0.178731485778e1, 0.341262936801,
                0.139479099785, 0.5076238131, -0.655600474113],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.00773692]*6}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexane of Starling (1973)",
        "__doc__":  u"""Starling, K.E., "Fluid Thermodynamic Properties for Light Petroleum Systems," Gulf Publishing Company, 1973.""",
        "R": 8.3159524,
        "cp": CP4,

        "Tmin": 222.04, "Tmax": 644.0, "Pmax": 55000.0, "rhomax": 8.6724844, 
        "Pmin": 0.001277, "rhomin": 8.8394, 

        "nr1": [0.261128818398e1, 0.451396780770, -0.783362300734,
                -0.108785843809e1, 0.124906986929, -0.155020819852e-1,
                0.42399441457, -0.636532521368, -0.524764104726e-1,
                0.120405133154e-1, 0.992632580157e-3],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.261128818398e1, -0.558196781075],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.42752599]*2}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.31451,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [2.43433265, 1.18137185, -4.24411947, 1.08655334e-1,
                2.87828538e-4, -2.51781047e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],
 
        "nr2": [2.16096570e-2, -4.58052979e-1, 1.63940974e-1, -2.55034034e-2,
                -2.47418231e-1, -8.05544799e-3, -7.78926202e-2, -2.69044742e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, GERG, helmholtz3, helmholtz4, helmholtz5

    _surface = {"sigma": [0.0522937, 0.0061685, -0.0035869],
                "exp": [1.25, 2.25, 3.25]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                   "a1": [30.18, 0.03], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [222.31, 232.62, -36872, -25733],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74172e1, 0.12897e1, -0.32544e1, -0.14609e1, 0.81765e-1],
        "exp": [1., 1.5, 3.1, 5.3, 5.6]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.14686e3, -0.26585e3, 0.12200e3],
        "exp": [0.75, 0.81, 0.88]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.13309, -0.50653e1, -0.11602e2, -0.28530e2, -0.51731e2, -0.13482e3],
        "exp": [0.107, 0.553, 2.006, 4.46, 8.0, 16.]}

    visco0 = {"eq": 2, "omega": 3,
              "__name__": "NIST",
              "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
              "ek": 399.3, "sigma": 0.5949,
              "n_chapman": 0.247780666/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-18.180936383994, 3263.1590998, 17.765176425, -53270.220915,
                    -0.32352381766, 195.54170454, 38519.153073],
              "rhoc": 2.704}

    visco1 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quinones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quinones-Cisneros (2006)",
              "Tref": 507.82, "muref": 1.0,
              "ek": 393.1, "sigma": 0.5949, "n_chapman": 0,
              "n_ideal": [16.9975, -54.2985, 48.0065],
              "t_ideal": [0, 0.25, 0.5],

              "a": [-6.63500718148775e-5, -2.14251735181008e-5, 7.74647275349291e-14],
              "b": [1.64280427908191e-4, -1.34908441238411e-4, -2.17284146069693e-14],
              "c": [7.25570985000000e-5, -3.12153040000000e-6, 0.0],
              "A": [1.45983786505096e-9, -8.15150058452202e-10, 0.0],
              "B": [2.59524353609885e-8, 1.69361972245028e-9, 0.0],
              "C": [-2.29226420147789e-6, 1.18011366260701e-6, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "NIST14",
               "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

               "Tref": 399.3, "kref": 1e-3,
               "no": [1.35558587, -0.143662461021788, 1],
               "co": [0, -1, -96],

               "Trefb": 507.35, "rhorefb": 2.704, "krefb": 1e-3,
               "nb": [15.275017704, 11.2896277792, -8.61369853497,
                      0.697714450907, 2.1687378215, -0.326193379046],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.0327e-9, "Tcref": 761.73}

    _thermal = thermo0,

if __name__ == "__main__":
    for eq in (0, 1, 4):
        st=nC6(T=300, P=1e5, eq=eq)
        print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
            st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
