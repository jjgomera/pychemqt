#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Ethanol(MEoS):
    """Multiparameter equation of state for ethanol"""
    name = "ethanol"
    CASNumber = "64-17-5"
    formula = "C2H6O"
    synonym = "R-170"
    rhoc = unidades.Density(275.9906)
    Tc = unidades.Temperature(513.9)
    Pc = unidades.Pressure(6148., "kPa")
    M = 46.06844  # g/mol
    Tt = unidades.Temperature(159)
    Tb = unidades.Temperature(351.39)
    f_acent = 0.644
    momentoDipolar = unidades.DipoleMoment(1.6909, "Debye")
    id = 134

    CP1 = {"ao": 6.4112,
           "an": [], "pow": [],
           "ao_exp": [1.95988750679, 7.60084166080, 3.89583440622, 4.23238091363],
           "exp": [694, 1549, 2911, 4659],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Dillon and Penoncello (2004)",
        "__doi__": {"autor": "Dillon, H.E. and Penoncello, S.G.",
                    "title": "A Fundamental Equation for Calculation of the Thermodynamic Properties of Ethanol", 
                    "ref": "Int. J. Thermophys., 25(2):321-335, 2004.",
                    "doi": "10.1023/B:IJOT.0000028470.49774.14"}, 
        "__test__": """
            >>> st=Ethanol(T=350, P=1e5)
            >>> print "%0.1f %0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0.1 350 737.85 259.95 1.0363 2.6542 3.1707 964.32
            >>> st=Ethanol(T=650, P=1e5)
            >>> print "%0.1f %0.0f %0.4g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0.1 650 0.8534 1771.7 4.8011 2.3679 2.5512 355.1
            
            >>> st=Ethanol(T=450, P=1e6)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            1 450 13.736 1270.3 3.4711 2.0217 2.3954 276.5
            
            >>> st=Ethanol(T=250, P=1e7)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            10 250 831.84 9.4714 0.1625 1.6724 2.0325 1372.8
            >>> st=Ethanol(T=600, P=1e7)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            10 600 129 1474.5 3.5258 2.7347 4.0688 273.66
            
            >>> st=Ethanol(T=400, P=1e8)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            100 400 784.52 500.06 1.3279 2.8005 3.2736 1348.2
            """, # Table IV, Pag 329
            
        "R": 8.314472,
        "cp": CP1,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 45800, "so": 180}, 

        "Tmin": 250.0, "Tmax": 650.0, "Pmax": 280000.0, "rhomax": 19.4, 
        "Pmin": 0.00000088, "rhomin": 19.4, 

        "nr1": [0.114008942201e2, -0.395227128302e2, 0.413063408370e2,
                -0.188892923721e2, 0.472310314140e1, -0.778322827052e-2,
                0.171707850032, -0.153758307602e1, 0.142405508571e1,
                0.132732097050, -0.114231649761, 0.327686088736e-5,
                0.495699527725e-3, -0.701090149558e-4, -0.225019381648e-5],
        "d1": [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 6, 7, 8, 8],
        "t1": [-0.5, 0, 0.5, 1.5, 2, 5, -0.5, 1, 2, 0, 2.5, 6, 2, 2, 4],

        "nr2": [-0.255406026981, -0.632036870646e-1, -0.314882729522e-1,
                0.256187828185e-1, -0.308694499382e-1, 0.722046283076e-2,
                0.299286406225e-2, 0.972795913095e-3],
        "d2": [1, 3, 3, 6, 7, 8, 2, 7],
        "t2": [5, 3, 7, 5.5, 4, 1, 22, 23],
        "c2": [2, 2, 2, 2, 2, 2, 4, 4],
        "gamma2": [1]*8}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Sun and Ely (2004).",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
            
        "R": 8.314472,
        "cp": CP1,
        "ref": {"name": "CUSTOM",
            "Tref": 273.15, "Pref": 1., "ho": 45800, "so": 180}, 

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 280000.0, "rhomax": 19.6, 
        "Pmin": 0.00000064, "rhomin": 19.55, 

        "nr1":  [-2.95455387, 1.95055493, -1.31612955, -1.47547651e-2,
                 1.39251945e-4, 5.04178939e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],
        
        "nr2": [2.52310166e-1, 1.97074652, 8.73146115e-1, 4.27767205e-2,
                9.68966545e-2, -8.39632113e-1, -7.71828521e-2, 1.63430744e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, helmholtz2
    _PR = 0.0043733

    _surface = {"sigma": [0.05], "exp": [0.952]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.81829e1, -0.62767, -0.33289e1, -0.35278e1, 0.93103e1],
        "exp": [1, 1.5, 3, 5.6, 7.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.11818e1, -0.36120e1, 0.54325e1, -0.47789, -0.17766e-1],
        "exp": [0.098, 0.22, 0.35, 0.7, 2.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.93315, -0.40761e2, 0.6325e2, -0.45195e2, 0.15114e1, -0.56666e2],
        "exp": [0.09, 1.07, 1.3, 1.7, 4.0, 5.0]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Kiselev (2005)",
              "__doc__": """Kiselev, S. B., Ely, J. F., Abdulagatov, I. M., Huber, M. L.,"Generalized SAFT-DFT/DMT Model for the Thermodynamic, Interfacial, and Transport Properties of Associating Fluids: Application for n-Alkanols", Ind. Eng. Chem. Res., 2005, 44, 6916-6927""",
              "ek": 362.6, "sigma": 0.453,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0,
              "n_ideal": [-1.03116, 3.48379e-2, -6.50264e-6],
              "t_ideal": [0, 1, 2],

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 362.6, "etaref_virial": 0.0559816,

              "Tref_res": 513.9, "rhoref_res": 5.991*M, "etaref_res": 1000,
              "n_packed": [-3.38264465, 1.27568864e1],
              "t_packed": [0, 0.5],
              "n_poly": [1.31194057e-1, -8.05700894e-2, -3.82240694e-1,
                         1.53811778e-1, 0.0, -1.10578307e-1, -2.37222995e1],
              "t_poly": [0, 0, -1, -1, -2, -2, 0],
              "d_poly": [2, 3, 2, 3, 2, 3, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 1],
              "n_num": [2.37222995e1],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Marsh (2002)",
               "__doc__": """Marsh, K., Perkins, R., and Ramires, M.L.V., "Measurement and Correlation of the Thermal Conductivity of Propane from 86 to 600 K at Pressures to 70 MPa," J. Chem. Eng. Data, 47(4):932-940, 2002""",

               "Tref": 513.9, "kref": 1,
               "no": [0.123120e-1, -0.153612e-1, 0.426611e-1],
               "co": [0, 1, 2],

               "Trefb": 513.9, "rhorefb": 5.991, "krefb": 1.,
               "nb": [0.266894e-1, 0.0, -0.482953e-1, 0.414022e-1, 0.172939e-1,
                      -0.977825e-2, 0.0, 0.0, 0.0, 0.0],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 2.307981e-10, "Tcref": 770.85}

    _thermal = thermo0,

# TODO: Add discard eq and thermal in refprop
