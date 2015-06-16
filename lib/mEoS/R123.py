#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R123(MEoS):
    """Multiparameter equation of state for R123"""
    name = "2,2-dichloro-1,1,1-trifluoroethane"
    CASNumber = "306-83-2"
    formula = "CHCl2CF3"
    synonym = "R123"
    rhoc = unidades.Density(550.)
    Tc = unidades.Temperature(456.831)
    Pc = unidades.Pressure(3661.8, "kPa")
    M = 152.931  # g/mol
    Tt = unidades.Temperature(166.0)
    Tb = unidades.Temperature(300.973)
    f_acent = 0.28192
    momentoDipolar = unidades.DipoleMoment(1.356, "Debye")
    id = 236
    # id = 1631

    CP1 = {"ao": 17.01154/8.31451,
           "an": [0.4046308/8.31451, -4.644803e-4/8.31451, 2.347418e-7/8.31451],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-123 of Younglove and McLinden (1994).",
        "__doi__": {"autor": "Younglove, B.A. and McLinden, M.O.",
                    "title": "An International Standard Equation of State for the Thermodynamic Properties of Refrigerant 123 (2,2-Dichloro-1,1,1-trifluoroethane)", 
                    "ref": "J. Phys. Chem. Ref. Data, 23:731-779, 1994.",
                    "doi":  "10.1063/1.555950"}, 
        #TODO: Add test from file
        #FIXME: The file include derived heltmholtz expresion for MBWR equations
        "R": 8.31451,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 40000.0, "rhomax": 11.60, 
        "Pmin": 0.0042, "rhomin": 11.60, 

        "b": [None, -0.657453133659e-2, 0.293479845842e1, -0.989140469845e2,
              0.201029776013e5, -0.383566527886e7, 0.227587641969e-2,
              -0.908726819450e1, 0.434181417995e4, 0.354116464954e7,
              -0.635394849670e-3, 0.320786715274e1, -0.131276484299e4,
              -0.116360713718, -0.113354409016e2, -0.537543457327e4,
              0.258112416120e1, -0.106148632128, 0.500026133667e2,
              -0.204326706346e1, -0.249438345685e7, -0.463962781113e9,
              -0.284903429588e6, 0.974392239902e10, -0.637314379308e4,
              0.314121189813e6, -0.145747968225e3, -0.843830261449e7,
              -0.241138441593e1, 0.108508031257e4, -0.106653193965e-1,
              -0.121343571084e2, -0.257510383240e3]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz transform of MBWR EOS for R-123 of Younglove & McLinden (1994)",
        "__doi__": {"autor": "Younglove, B.A. and McLinden, M.O.",
                    "title": "An International Standard Equation of State for the Thermodynamic Properties of Refrigerant 123 (2,2-Dichloro-1,1,1-trifluoroethane)", 
                    "ref": "J. Phys. Chem. Ref. Data, 23:731-779, 1994.",
                    "doi":  "10.1063/1.555950"}, 
                    
        "R": 8.31451,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 40000.0, "rhomax": 11.60, 
        "Pmin": 0.0042, "rhomin": 11.60, 

        "nr1": [-0.100242647494e2, -0.280607656419, 0.206814471606e-1,
                -0.284379431451, 0.593928110321e1, -0.936560389528e1,
                0.416660793675e1, -0.174023292951e1, 0.177019905365,
                -0.154721692260e1, 0.161820495590e1, 0.288903529383e1,
                -0.118493874757, 0.130952266209e1, -0.117308103711e1,
                -0.128125131950, -0.786087387513e-1, -0.816000499305e-1,
                0.536451054311e-1, -0.680078211929e-2, 0.701264082191e-2,
                -0.901762397311e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8],
        "t1": [3, 4, 5, 0, 0.5, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 1, 2, 3, 2, 2, 3, 3],

        "nr2": [0.100242647494e2, 0.280607656419, -0.206814471606e-1,
                0.798923878145e1, -0.547972072476, -0.206814470584e-1,
                0.249142724365e1, -0.273986034884, 0.236001863614,
                0.540528251211, -0.600457561959e-1, 0.786672874826e-1,
                0.708085874508e-1, -0.150114389748e-1, 0.182205199477e-2,
                0.314978575163e-2, 0.784455573794e-2, 0.364410397155e-3],
        "d2": [0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8, 8, 10, 10, 10],
        "t2": [3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5],
        "c2": [2]*18,
        "gamma2": [1]*18}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-123 of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids", 
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"}, 
        "__test__": """
            >>> st=R123(T=700, rho=200, eq=2)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            0.8667 6.018 1.9509
            >>> st2=R123(T=750, rho=100, eq=2)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            144.33 0.29582
            """, # Table III, Pag 117
            
        "R": 8.31451,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 11.62, 
        "Pmin": 0.0041534, "rhomin": 11.613, 

        "nr1": [0.1116973e1, -0.3074593e1, 0.51063873, 0.94478812e-1,
                0.29532752e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.66974438, 0.96438575, -0.14865424e-1, -0.49221959,
                -0.22831038e-1, -0.1407486, -0.25117301e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = MBWR, helmholtz1, helmholtz2

    _surface = {"sigma": [0.056151], "exp": [1.2367]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74610e1, 0.20293e1, -0.21897e1, -0.34945e1],
        "exp": [1.0, 1.5, 2.25, 4.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.19996e1, 0.41823, 0.24849, 0.18831, 0.13737],
        "exp": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.30205e1, -0.74537e1, -0.21880e2, -0.57430e2, 0.11239e2, -0.16640e3],
        "exp": [0.3905, 1.29, 3.4, 7.0, 12.0, 15.0]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Tanaka (1996)",
              "__doc__": """Tanaka, Y. and Sotani, T., "Transport Properties (Thermal Conductivity and Viscosity), Int. J. Thermophys., 17(2):293-328, 1996""",
              "ek": 275.16, "sigma": 0.5909,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0,
              "n_ideal": [5.099859e-2, -2.402786e-5],
              "t_ideal": [1, 2],

              "Tref_res": 1., "rhoref_res": 6.538897e-3*M, "etaref_res": 1,
              "n_packed": [1.828263e3],
              "t_packed": [0],
              "n_poly": [-1.762849e2, -2.226484e-2, 5.550623e-5, -1.009812e-1,
                         6.161902e-5, -8.840480e-8],
              "t_poly": [0, 0, 1, 0, 0, 0],
              "d_poly": [0, 1, 1, 1, 2, 3],
              "g_poly": [0, 0, 0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0, 0, 0],
              "n_num": [-3.222951e5],
              "t_num": [0],
              "d_num": [0],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [1, 0],
              "g_den": [0, 1],
              "c_den": [0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Laesecke (1996)",
               "__doc__": """Laesecke, A., Perkins, R.A., and Howley, J.B., "An improved correlation for the thermal conductivity of HCFC123 (2,2-dichloro-1,1,1-trifluoroethane)," Int. J. Refrigeration, 19:231-238, 1996""",

               "Tref": 1., "kref": 1,
               "no": [-0.00778, 5.695e-5],
               "co": [0, 1],

               "Trefb": 456.831, "rhorefb": 3.596417, "krefb": 1,
               "nb": [0.642894e-1, -0.530474e-1, 0.453522e-4, -0.139928,
                      0.16654, -0.162656e-1, 0.136819, -0.183291, 0.357146e-1,
                      -0.231210e-1, 0.341945e-1, -0.757341e-2],
               "tb": [-1.5, -2, -6, 0, -0.5, -1.5, 0, -0.5, -1.5, 0, -0.5, -1.5],
               "db": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
               "cb": [0]*12,

               "critical": 1,
               "crit_num_Tref": -456.831, "crit_num_rhoref": 3.596417,
               "crit_num_k": 1.,
               "crit_num_n": [0.486742e-2],
               "crit_num_alfa": [0],
               "crit_num_t": [0],
               "crit_num_beta": [0],
               "crit_num_d": [0],
               "crit_num_c": [0],
               "crit_den_n": [-100.0, -7.08535],
               "crit_den_alfa": [-1, 0],
               "crit_den_t": [4, 0],
               "crit_den_beta": [0, -1],
               "crit_den_d": [0, 2],
               "crit_den_c": [0, 0]}

    _thermal = thermo0,
