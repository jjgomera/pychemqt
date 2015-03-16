#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class nC4(MEoS):
    """Multiparameter equation of state for n-butane

    >>> butano=nC4(h=50.470, P=25)
    >>> print "%0.1f %0.2f %0.4f %0.3f %0.5f %0.4f %0.4f %0.2f" % (butano.T, butano.rho, butano.u.kJkg, butano.h.kJkg, butano.s.kJkgK, butano.cv.kJkgK, butano.cp.kJkgK, butano.w)
    450.0 463.97 -3.4126 50.470 -0.35405 2.3235 2.9759 632.80
    """
    name = "n-butane"
    CASNumber = "106-97-8"
    formula = "CH3-(CH2)2-CH3"
    synonym = "R-600"
    rhoc = unidades.Density(228.)
    Tc = unidades.Temperature(425.125)
    Pc = unidades.Pressure(3796.0, "kPa")
    M = 58.1222  # g/mol
    Tt = unidades.Temperature(134.895)
    Tb = unidades.Temperature(272.660)
    f_acent = 0.201
    momentoDipolar = unidades.DipoleMoment(0.05, "Debye")
    id = 6
    _Tr = unidades.Temperature(406.785141)
    _rhor = unidades.Density(230.384826)
    _w = 0.194240287

    CP1 = {"ao": 4.24680487,
           "an": [], "pow": [],
           "ao_exp": [5.54913289, 11.4648996, 7.59987584, 9.66033239],
           "exp": [329.40404, 1420.17366, 2113.08938, 4240.8573],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 4.33944,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [9.44893, 6.89406, 24.4618, 14.7824],
           "hyp": [1.101487798*Tc, 0.43195766*Tc, 4.502440459*Tc, 2.124516319*Tc]}

    CP3 = {"ao": 4.240207,
           "an": [], "pow": [],
           "ao_exp": [5.513671, 7.388450, 10.250630, 11.061010],
           "exp": [327.55988, 1319.06935, 4138.63184, 1864.36783],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": -1.3491511376e1,
           "an": [3.8802310194e5, -1.5444296890e5, 2.8455082239e3,
                  6.6142595353e-2, -2.4307965028e-5, 1.5044248429e-10],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-8.3933423467], "exp": [3000],
           "ao_hyp": [], "hyp": []}

    CP5 = {"ao": 4.33944,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.2071931e7, 0.2324827e6, 0.8962262e8, 0.1205864e8],
           "hyp": [0.46827e3, 0.183636e3, 0.19141e4, 0.903185e3]}

    CP6 = {"ao": 0.801601/8.3143*58.124,
           "an": [0.655936e-3/8.3143*58.124, 0.12277e-4/8.3143*58.124,
                  -0.165626e-7/8.3143*58.124, 0.67736e-11/8.3143*58.124],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Buecker and Wagner (2006)",
        "__doc__":  u"""Bücker, D., Wagner, W. Reference equations of state for the thermodynamic properties of fluid phase n-butane and isobutane. J. Phys. Chem. Ref. Data 35 (2006), 929 – 1020.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 575., "Pmax": 200000.0, "rhomax": 13.86, 
        "Pmin": 0.000653, "rhomin": 12.645, 

        "nr1": [0.25536998241635e1, -0.44585951806696e1, 0.82425886369063,
                0.11215007011442, -0.35910933680333e-1, 0.16790508518103e-1,
                0.32734072508724e-1],
        "d1": [1, 1, 1, 2, 3, 4, 4],
        "t1": [0.50, 1.00, 1.50, 0.00, 0.50, 0.50, 0.75],
        "nr2": [0.95571232982005, -0.10003385753419e1, 0.85581548803855e-1,
                -0.25147918369616e-1, -0.15202958578918e-2, 0.47060682326420e-2,
                -0.97845414174006e-1, -0.48317904158760e-1, 0.17841271865468,
                0.18173836739334e-1, -0.11399068074953, 0.19329896666669e-1,
                0.11575877401010e-2, 0.15253808698116e-3, -0.43688558458471e-1,
                -0.82403190629989e-2],
        "d2": [1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6],
        "t2": [2.00, 2.50, 2.50, 1.50, 1.00, 1.50, 4.00, 7.00, 3.00, 7.00,
               3.00, 1.00, 6.00, 0.00, 6.00, 13.00],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*16,

        "nr3": [-0.28390056949441e-1, 0.14904666224681e-2],
        "d3": [1, 2],
        "t3": [2., 0.],
        "alfa3": [10, 10],
        "beta3": [150, 200],
        "gamma3": [1.16, 1.13],
        "epsilon3": [0.85, 1.]}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for butane of Younglove and Ely (1987)",
        "__doc__":  u"""Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
        "R": 8.31434,
        "cp": CP4,

        "Tmin": 134.86, "Tmax": 500., "Pmax": 70000.0, "rhomax": 13.2, 
        "Pmin": 6.736e-4, "rhomin": 12.65, 

        "b": [None, 0.153740104603e-1, -0.160980034611, -0.979782459010e1,
              0.499660674504e3, -0.102115607687e7, 0.236032147756e-2,
              -0.137475757093e1, -0.907038733865e3, 0.385421748213e6,
              -0.349453710700e-4, 0.157361122714, 0.102301474068e3,
              0.182335737331e-1, -0.404114307787e1, 0.187979855783e1,
              0.362088795040, -0.738762248266e-2, -0.218618590563e1,
              0.118802729027, 0.706854198713e6, -0.219469885796e9,
              -0.182454361268e5, 0.206790377277e10, 0.111757550145e3,
              0.558779925986e5, -0.159579054026e2, -0.148034214622e7,
              -0.245206328201, 0.218305259309e3, -0.923990627338e-4,
              -0.205267776639e1, 0.387639044820e2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 575., "Pmax": 69000.0, "rhomax": 13.2, 
        "Pmin": 0.000653, "rhomin": 12.645, 

        "nr1": [0.10626277411455e1, -0.28620951828350e1, 0.88738233403777,
                -0.12570581155345, 0.10286308708106, 0.25358040602654e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.32325200233982, -0.37950761057432e-1, -0.32534802014452,
                -0.79050969051011e-1, -0.20636720547775e-1, 0.57053809334750e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Miyamoto and Watanabe (2001)",
        "__doc__":  u"""Miyamoto, H. and Watanabe, K. "A Thermodynamic Property Model for Fluid-Phase n-Butane," Int. J. Thermophys., 22(2):459-475, 2001.""",
        "R": 8.314472,
        "cp": CP3,

        "Tmin": 134.87, "Tmax": 589., "Pmax": 69000.0, "rhomax": 13.15, 
        "Pmin": 0.000688, "rhomin": 12.652, 

        "nr1": [2.952054e-1, -1.32636, -2.031317e-3, 2.240301e-1,
                -3.635425e-2, 1.905841e-3, 7.409154e-5, -1.401175e-6],
        "d1": [1, 1, 2, 2, 3, 5, 8, 8],
        "t1": [-0.25, 1.5, -0.75, 0, 1.25, 1.5, 0.5, 2.5],

        "nr2": [-2.492172, 2.386920, 1.424009e-3, -9.393388e-3, 2.616590e-3,
                -1.977323e-1, -3.809534e-2, 1.523948e-3, -2.391345e-2,
                -9.535229e-3, 3.928384e-5],
        "d2": [3, 3, 8, 5, 6, 1, 5, 7, 2, 3, 15],
        "t2": [1.5, 1.75, -0.25, 3, 3, 4, 2, -1, 2, 19, 5],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*11}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for butane of Span and Wagner (2003)",
        "__doc__":  u"""Span, R., Wagner, W. Equations of state for technical applications. II. Results for nonpolar fluids. Int. J. Thermophys. 24 (2003), 41 – 109.""",
        "R": 8.31451,
        "cp": CP5,

        "Tmin": 134.86, "Tmax": 600., "Pmax": 100000.0, "rhomax": 13.20, 
        "Pmin": 0.00064578, "rhomin": 12.671, 

        "nr1": [0.10626277e1, -0.28620952e1, 0.88738233, -0.12570581,
                0.10286309, 0.25358041e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.323252, -0.37950761e-1, -0.32534802, -0.79050969e-1,
                -0.20636721e-1, 0.57053809e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP6,

        "Tmin": 140.0, "Tmax": 589., "Pmax": 30000.0, "rhomax": 12.81, 
        "Pmin": 0.00161, "rhomin": 12.573, 

        "nr1": [-0.504188295325, 0.541067401063, -0.760421383062e-1,
                0.846035653528, -0.191317317203e1, 0.521441860186,
                -0.783511318207, 0.689697797175e-1, 0.947825461055e-1,
                -0.141401831669, 0.382675021672, -0.423893176684e-1,
                0.677591792029e-1, 0.567943363340e-1, -0.131517698401,
                0.221136942526e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.504188295325, -0.541067401063, 0.760421383062e-1,
                -0.619109535460e-1, 0.423035373804, -0.390505508895],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.08974964]*6}

    helmholtz6 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [1.18936994, 1.05407451, -3.24964532, 8.25263908e-2,
                2.76467405e-4, -8.09869214e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-9.38097492e-2, 1.46213532e-1, 4.01168502e-1, -1.28716120e-2,
                -2.75191070e-1, -1.62708971e-2, -7.04082962e-2, -2.32871995e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4, helmholtz5, helmholtz6

    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.0557549],  "expt0": [-1.], "expd0": [1.],
                   "a1": [20.611, 0.02], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [66.64, 24.44, -7461.2, -1983.6],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.00066566,
                "Tmin": 134.895, "Tmax": 575.0,
                "a1": [-558558235.4, 558558236.4], "exp1": [0, 2.206],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _surface = {"sigma": [0.05418], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.71897e1, 0.26122e1, -0.21729e1, -0.27230e1],
        "exp": [1, 1.5, 2., 4.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.52341e1, -0.62011e1, 0.36063e1, 0.22137],
        "exp": [0.44, 0.6, 0.76, 5.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.27390e1, -0.57347e1, -0.16408e2, -0.46986e2, -0.10090e3],
        "exp": [0.39, 1.14, 3.0, 6.5, 14.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.17067154, -0.48879666, 0.039038856],
              "__name__": "Vogel (1999)",
              "__doc__": """Vogel, E., Kuechenmeister, C., and Bich, E., "Viscosity for n-Butane in the Fluid Region," High Temp. - High Pressures, 31(2):173-186, 1999.""",
              "ek": 280.51, "sigma": 0.57335,
              "Tref": 1, "rhoref": 1.*M, "etaref": 1.,
              "n_chapman": 0.1628213/M**0.5,

              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158],
              "t_virial": [0.0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 280.51, "etaref_virial": 0.1135034,

              "Tref_res": 425.125, "rhoref_res": 3.92*M, "etaref_res": 1,
              "n_packed": [2.30873963359, 2.03404037254],
              "t_packed": [0, 0.5],
              "n_poly": [-54.7737770846, 58.0898623034, 0, 35.2658446259,
                         -39.6682203832, 0, -1.83729542151, 0, 0,
                         -0.833262985358, 1.93837020663, 0, -188.075903903],
              "t_poly": [0, -1, -2, 0, -1, -2, 0, -1, -2, 0, -1, -2, 0],
              "d_poly": [2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              "n_num": [188.075903903],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [1, 0]}

    visco1 = {"eq": 2, "omega": 2,
              "__name__": "Younglove (1987)",
              "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
              "ek": 440., "sigma": 0.503103,
              "n_chapman": 0.20352457/M**0.5,
              "F": [0.1630521851e1, 0.0, 1.40, 425.16],
              "E": [-0.2724386845e2, 0.8012766611e3, 0.2503978646e2,
                    -0.1309704275e5, -0.8313305258e-1, 0.6636975027e2,
                    0.9849317662e4],
              "rhoc": 3.920}

    visco2 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "Tref": 425.125, "muref": 1.0,
              "ek": 440., "sigma": 0.503103, "n_chapman": 0,
              "n_ideal": [18.3983, -57.1255, 49.3197],
              "t_ideal": [0, 0.25, 0.5],
              "a": [-1.34110938674421e-5, -8.56587924603951e-5, -6.45720639242339e-13],
              "b": [1.49859653515567e-4, -1.71133855507542e-4, 7.37953726544736e-13],
              "c": [3.53018109777015e-7, -1.93040375218067e-5, -1.26469933968355e-14],
              "A": [-3.63389393526204e-9, -7.73717469888952e-10, 0.0],
              "B": [3.70980259815724e-8, 2.07658634467549e-9, 0.0],
              "C": [-1.12495594619911e-7, 7.66906137372152e-8, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2002)",
               "__doc__": """Perkins, R.A, Ramires, M.L.V., Nieto de Castro, C.A. and Cusco, L., "Measurement and Correlation of the Thermal Conductivity of Butane from 135 K to 600 K at Pressures to 70 MPa," J. Chem. Eng. Data, 47(5):1263-1271, 2002.""",

               "Tref": 425.16, "kref": 1.,
               "no": [1.62676e-3, 9.75703e-4, 2.89887e-2],
               "co": [0, 1, 2],

               "Trefb": 425.16, "rhorefb": 3.92, "krefb": 1.,
               "nb": [-3.04337e-2, 4.18357e-2, 1.65820e-1, -1.47163e-1,
                      -1.48144e-1, 1.33542e-1, 5.25500e-2, -4.85489e-2,
                      -6.29367e-3, 6.44307e-3],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.875350e-9, "Tcref": 637.68}

    thermo1 = {"eq": 2, "omega": 2,
               "__name__": "Younglove (1987)",
               "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
               "visco": visco1,
               "n_chapman": 2.0352526600e-1,
               "G": [0.1530992335e1, -0.2114511021],
               "E": [0.4024170074e-2, 0.1561435847e1, -0.6004381127e3,
                     -0.7547260841e-3, -0.2069676662e-1, 0.9382534978e2,
                     -0.1711371457, 0.3647724935e2],

               "critical": 2,
               "X": [0.000769608, 13.2533, 0.485554, 1.01021],
               "Z": 9.10218e-10}

    _thermal = thermo0, thermo1

if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=nC4(T=300., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
