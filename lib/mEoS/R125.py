#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R125(MEoS):
    """Multiparameter equation of state for R125

    >>> r125=R125(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r125.T, r125.rho, r125.u.kJkg, r125.h.kJkg, r125.s.kJkgK, r125.cv.kJkgK, r125.cp.kJkgK, r125.w)
    300.0 4.8836 342.90 363.38 1.7177 0.72455 0.79910 149.16
    """
    name = "pentafluoroethane"
    CASNumber = "354-33-6"
    formula = "CHF2CF3"
    synonym = "R125"
    rhoc = unidades.Density(573.5822706)
    Tc = unidades.Temperature(339.173)
    Pc = unidades.Pressure(3617.7, "kPa")
    M = 120.0214  # g/mol
    Tt = unidades.Temperature(172.52)
    Tb = unidades.Temperature(225.06)
    f_acent = 0.3052
    momentoDipolar = unidades.DipoleMoment(1.563, "Debye")
    id = 236
    # id = 1231

    CP1 = {"ao": 0,
           "an": [3.063], "pow": [0.1],
           "ao_exp": [2.303, 5.086, 7.3], "exp": [314, 756, 1707],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 5.911212,
           "an": [], "pow": [],
           "ao_exp": [6.856764, 4.885985, 3.292859],
           "exp": [670.10271, 1626.81842, 1863.11546],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 4.3987,
           "an": [0.0242728, -0.4099e-5], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 25.87069,
           "an": [0.2690914, -1.331388e-4, 4.10133e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP5 = {"ao": 3.111514,
           "an": [10.982115/339.33, -1.843797/339.33**2, 0.019273/339.33**3],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP6 = {"ao": 12.990267052,
           "an": [-5.2262199/339.165**-0.25], "pow": [-0.25],
           "ao_exp": [7.028445731, 4.586635360],
           "exp": [1664.325535, 705.7406967],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Lemmon and Jacobsen (2005).",
        "__doc__":  u"""Lemmon, E.W. and Jacobsen, R.T, "A New Functional Form and New Fitting Techniques for Equations of State with Application to Pentafluoroethane (HFC-125)," J. Phys. Chem. Ref. Data, 34(1):69-108, 2005.""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.09, 
        "Pmin": 2.914, "rhomin": 14.086, 

        "nr1": [0.5280760e1, -0.8676580e1, 0.7501127, 0.7590023, 0.1451899e-1],
        "d1": [1, 1, 1, 2, 4],
        "t1": [0.669, 1.05, 2.75, 0.956, 1],

        "nr2": [0.4777189e1, -0.3330988e1, 0.3775673e1, -0.2290919e1,
                0.8888268, -0.6234864, -0.4127263e-1, -0.8455389e-1,
                -0.1308752, 0.8344962e-2],
        "d2": [1, 1, 2, 2, 3, 4, 5, 1, 5, 1],
        "t2": [2, 2.75, 2.38, 3.37, 3.47, 2.63, 3.45, 0.72, 4.23, 0.2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 3],
        "gamma2": [1]*10,

        "nr3": [-0.1532005e1, -0.5883649e-1, 0.2296658e-1],
        "t3": [4.5, 29, 24],
        "d3": [2, 3, 5],
        "alfa3": [1]*3,
        "beta3": [1]*3,
        "gamma3": [0]*3,
        "epsilon3": [0]*3}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-125 of Outcalt and McLinden (1995)",
        "__doc__":  u"""Outcalt, S.L. and McLinden, M.O., "Equations of state for the thermodynamic properties of R32 (difluoromethane) and R125 (pentafluoroethane)," Int. J. Thermophysics, 16:79-89, 1995.""",
        "R": 8.314471,
        "cp": CP4,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.10, 
        "Pmin": 2.921, "rhomin": 14.095, 

        "b": [None, -0.523369607050e-1, 0.378761878904e1, -0.807152818990e2,
              0.115654605248e5, -0.152175619161e7, 0.597541484451e-2,
              -0.145990589966e1, -0.992338995652e3, -0.399180535687e6,
              -0.722591037504e-3, 0.358108080969, -0.108627994573e3,
              0.229821626570e-1, 0.149537670449e1, 0.911199833952e3,
              -0.254479949722, 0.102433894096e-1, -0.645583164735e1,
              0.218649963191, 0.114748721552e7, -0.118389825386e9,
              0.306539775027e5, 0.542870289406e9, 0.903502635609e3,
              -0.153646507435e6, 0.314617903718e1, 0.429297546671e6,
              0.109652021582, -0.329350271819e2, -0.338796950505e-3,
              0.384533651902, -0.491511706857e2]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Sunaga et al. (1998)",
        "__doc__":  u"""Sunaga, H., Tillner-Roth, R., Sato, H., and Watanabe, K., "A Thermodynamic Equation of State for Pentafluoroethane (R-125)," Int. J. Thermophys., 19(6):1623-1635, 1998.""",
        "R": 8.314471,
        "cp": CP2,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.09, 
        "Pmin": 2.943, "rhomin": 14.088, 

        "nr1": [0.12439220, 0.27922179, -0.11822597e1, 0.23616512, -0.11571810e-1],
        "d1": [1, 2, 2, 3, 2],
        "t1": [-0.5, 0, 1.5, 1.5, 3],

        "nr2": [0.1225177e1, -0.2147964e1, -0.298138, 0.3391211, -0.6322995e-3,
                0.1271747e-3, 0.5026962e-5, -0.1667058, -0.733275e-1,
                -0.637878e-1, 0.683311e-5, -0.1995426e-1, 0.1260026e-1],
        "d2": [1, 1, 1, 3, 8, 10, 12, 1, 2, 4, 15, 3, 4],
        "t2": [0.5, 1, 3, 2.75, 2, -1, 1.25, 4, 4, 3, 0.25, 23, 14],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*13}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Piao and Noguchi (1998)",
        "__doc__":  u"""Piao, C.-C. and Noguchi, M., "An international standard equation of state for the thermodynamic properties of HFC-125 (pentafluoroethane)," J. Phys. Chem. Ref. Data, 27(4):775-806, 1998.""",
        "R": 8.314471,
        "cp": CP3,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.11, 
        "Pmin": 2.9562, "rhomin": 14.1, 

        "nr1": [0.85393382372e-1, -0.133260499658, 0.257817782488,
                -0.735018179542, -0.787454743426, -0.190320468891e-1,
                0.388329449013, -0.631901774641, 0.623842653447, 0.109925047828,
                -0.993099630896e-1, -0.104601585904e-1, -0.769998709731e-1,
                0.149829594347e-1, 0.166640927925e-1, -0.181492321758e-2],
        "d1": [0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6],
        "t1": [1, 2, 0, 1, 3, 4, 0, 1, 3, 1, 3, 1, 2, 1, 2, 1],

        "nr2": [-0.85393382372e-1, 0.133260499658, 0.410983574575, -0.45298892633],
        "d2": [0, 0, 2, 2],
        "t2": [1, 2, 1, 2],
        "c2": [2, 2, 2, 2],
        "gamma2": [1]*4}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-125 of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids", 
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"}, 
        "__test__": """
            >>> st=R125(T=700, rho=200, eq=4)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            1.0745 14.620 1.3219
            >>> st2=R125(T=750, rho=100, eq=4)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            151.30 0.35860
            """, # Table III, Pag 117
            
        "R": 8.31451,
        "cp": CP5,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 14.1, 
        "Pmin": 2.9213, "rhomin": 14.096, 

        "nr1": [0.11290996e1, -0.28349269e1, 0.29968733, 0.87282204e-1,
                0.26347747e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.61056963, 0.90073581, -0.68788457e-2, -0.44211186,
                -0.35041493e-1, -0.1269863, -0.25185874e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Astina and Sato (2004)",
        "__doc__":  u"""Astina, I.M. and Sato, H. "A Rational Fundamental Equation of State for Pentafluoroethane with Theoretical and Experimental Bases," Int. J. Thermophys., 25(1):113-131, 2004.""",
        "R": 8.314472,
        "cp": CP6,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 70000.0, "rhomax": 14.1, 
        "Pmin": 2.94, "rhomin": 14.1, 

        "nr1": [1.51628822, -1.4959805, -1.2893965, 1.47295195, -2.22976436,
                1.02082011, -9.61695881e-3, 4.14142522e-2],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4],
        "t1": [0.5, 0.75, 2.25, 0.5, 0.875, 2, 3, 0.5],

        "nr2": [1.46217490e-1, -6.56486371e-2, -9.18319727e-2, -2.90343386e-2,
                -1.74343357e-2, -8.77406498e-4, -5.10648362e-3, 3.52425947e-3,
                4.98022850e-4],
        "d2": [3, 6, 4, 2, 4, 4, 3, 5, 7],
        "t2": [4, 2, 3.25, 9.5, 4.5, 10.5, 25, 5, 28],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    helmholtz6 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [7.41057508e-1, 1.13555445, -3.12563760, 9.32031442e-2,
                2.76844975e-4, -5.64403707e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [9.63969526e-3, 4.30480259e-1, 7.65668079e-1, -1.13913859e-2,
                -4.41468178e-1, -2.00943884e-2, -1.26041587e-1, -2.32331768e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, helmholtz2, helmholtz3, helmholtz4, helmholtz5, helmholtz6
    _PR = -0.00247

    _surface = {"sigma": [0.05260], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.5295, 1.9026, -2.2966, -3.448],
        "exp": [1, 1.5, 2.3, 4.6]}
    _liquid_Density = {
        "eq": 2,
        "ao": [1.6684, 0.88415, 0.44383],
        "exp": [1, 1.8, 8.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-2.8403, -7.2738, -21.89, -58.825],
        "exp": [0.38, 1.22, 3.3, 6.9]}

    visco0 = {"eq": 1, "omega": 3,
              "__name__": "Huber (2006)",
              "__doc__": """Huber, M.L., and Laesecke, A. , "Correlation for the Viscosity of Pentafluoroethane (R125) from the Triple Point to 500 K at Pressures up to 60 MPa", Ind. Eng. Chem. Res,2006, 45, 4447-4453""",
              "ek": 237.077, "sigma": 0.5235,
              "collision": [0.355404, -0.464337, 0.257353e-1],
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.2924206/M**0.5,

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 237.077, "etaref_virial": 0.0863974,

              "Tref_res": 339.173, "rhoref_res": 4.779*M, "etaref_res": 1000,
              "n_packed": [3.03379692, 0.299246403],
              "t_packed": [0, 0.5],
              "n_poly": [0.0, -5.09666198e-3, 5.67744840e-3, 0.0, -0.141256365],
              "t_poly": [-1, -1, -2, -2, 0],
              "d_poly": [2, 3, 2, 3, 1],
              "g_poly": [0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0],
              "n_num": [0.141256365],
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
               "__name__": "Perkins (2006)",
               "__doc__": """Perkins, R.A. and Huber, M.L.,"Measurement and Correlation of the Thermal Conductivity of Pentafluoroethane (R125)from 190 K to 512 K at pressures to 70 MPa", Journal of Chemical and Engineering Data, 2006, 51, 898-904""",

               "Tref": 339.173, "kref": 1.,
               "no": [-0.460820e-2, 0.168688e-1, 0.488345e-2],
               "co": [0, 1, 2],

               "Trefb": 339.173, "rhorefb": 4.779, "krefb": 1.,
               "nb": [-0.729410e-2, 0.110497e-1, 0.416339e-1, -0.289236e-1,
                      -0.311487e-1, 0.278399e-1, 0.112682e-1, -0.121100e-1,
                      -0.138322e-2, 0.211196e-2],
               "tb": [0, 1]*5,
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.834646e-10, "Tcref": 508.7475}

    _thermal = thermo0,
