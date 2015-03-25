#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R113(MEoS):
    """Multiparameter equation of state for R113

    >>> r113=R113(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r113.T, r113.rho, r113.u.kJkg, r113.h.kJkg, r113.s.kJkgK, r113.cv.kJkgK, r113.cp.kJkgK, r113.w)
    300.0 1558.7841 -33.86 -33.79 0.1152 0.67117 0.91951 693.14
    """
    name = "1,1,2-trichloro-1,2,2-trifluoroethane"
    CASNumber = "76-13-1"
    formula = "CCl2FCClF2"
    synonym = "R113"
    rhoc = unidades.Density(560.)
    Tc = unidades.Temperature(487.21)
    Pc = unidades.Pressure(3392.2, "kPa")
    M = 187.375  # g/mol
    Tt = unidades.Temperature(236.93)
    Tb = unidades.Temperature(320.735)
    f_acent = 0.25253
    momentoDipolar = unidades.DipoleMoment(0.803, "Debye")
    id = 232

    CP1 = {"ao": 3.9999966,
           "an": [], "pow": [],
           "ao_exp": [12.4464495, 2.72181845, 0.692712415, 3.32248298],
           "exp": [5.1143280e2, 1.60676324e3, 4.20292102e3, 1.60618738e3],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-113 of Marx et al. (1992).",
        "__doc__":  u"""Marx, V., Pruss, A., and Wagner, W.,"Neue Zustandsgleichungen fuer R 12, R 22, R 11 und R 113.  Beschreibung des thermodynamishchen Zustandsverhaltens bei Temperaturen bis 525 K und Druecken bis 200 MPa," Duesseldorf: VDI Verlag, Series 19 (Waermetechnik/Kaeltetechnik), No. 57, 1992.""",
        "R": 8.314471,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 525.0, "Pmax": 200000.0, "rhomax": 9.10, 
        "Pmin": 1.87, "rhomin": 9.099, 

        "nr1": [0.8432092286, -0.2019185967e1, 0.2920612996, 0.5323107661e-1,
                0.3214971931e-2, 0.4667858574e-4, -0.1227522799e-5],
        "d1": [1, 1, 2, 3, 4, 8, 8],
        "t1": [0.5, 1.5, 1.5, -0.5, 2, 0, 3],

        "nr2": [0.8167288718, -0.1340790803e1, 0.4065752705, -0.1534754634,
                -0.2414435149e-1, -0.2113056197e-1, -0.3565436205e-1,
                0.1364654968e-2, -0.1251838755e-1, -0.1385761351e-2,
                0.7206335486e-3],
        "d2": [3, 3, 3, 5, 1, 2, 2, 9, 3, 7, 8],
        "t2": [-0.5, 0, 2, 1.5, 6, 2, 10, 6, 18, 15, 33],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4],
        "gamma2": [1]*11}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-113 of Span and Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids", 
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"}, 
        "__test__": """
            >>> st=R113(T=700, rho=200, eq=1)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            0.8055 3.962 4.0344
            >>> st2=R113(T=750, rho=100, eq=1)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            131.00 0.26004
            """, # Table III, Pag 117
            
        "R": 8.31451,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 9.09, 
        "Pmin": 1.869, "rhomin": 9.0893, 

        "nr1": [0.10519071e1, -0.28724742e1, 0.41983153, 0.87107788e-1,
                0.24105194e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.70738262, 0.93513411, -0.96713512e-2, -0.52595315,
                0.22691984e-1, -0.14556325, -0.2741995e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.05566], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73838e1, 0.32594e1, -0.27761e1, -0.37758e1, -0.19921],
        "exp": [1.0, 1.5, 1.8, 4.3, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.15785e1, 0.12404e1, -0.66933, 0.49775e1, -0.55253e1],
        "exp": [0.3, 0.7, 2.0, 4.0, 5.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.26225e1, -0.60753e1, -0.15768e2, -0.42361e2, -0.79071e1, -0.31966e3],
        "exp": [0.379, 1.13, 2.9, 6.0, 7.0, 15.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.355404, -0.464337, 0.257353e-1],
              "__name__": "Huber (2003)",
              "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., "Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a", Ind. Eng. Chem. Res., 42:3163-3178, 2003.""",
              "ek": 376.035, "sigma": 0.6019,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.2509/M**0.5,

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 376.035, "etaref_virial": 0.13132,

              "Tref_res": 487.21, "rhoref_res": 2.988659*M, "etaref_res": 1310,
              "n_packed": [3.16369563558749, -0.8901733752064137e-1,
                           0.1000352946668359],
              "t_packed": [0, 1, 2],
              "n_poly": [-0.2069007192080741e-1, 0.3560295489828222e-3,
                         0.2111018162451597e-2, 0.1396014148308975e-1,
                         -0.4564350196734897e-2, -0.3515932745836890e-2,
                         -0.2147633195397038],
              "t_poly": [0, -6, -2, -0.5, 2, 0, 0],
              "d_poly": [1, 2, 2, 2, 2, 3, 0],
              "g_poly": [0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0],
              "n_num": [0.2147633195397038],
              "t_num": [0],
              "d_num": [0],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2000)",
               "__doc__": """Perkins, R.A., Laesecke, A., Howley, J., Ramires, M.L.V., Gurova, A.N., and Cusco, L., "Experimental thermal conductivity values for the IUPAC round-robin sample of 1,1,1,2-tetrafluoroethane (R134a)," NISTIR, 2000.""",

               "Tref": 487.21, "kref": 1.1,
               "no": [-0.460820e-2, 0.168688e-1, 0.488345e-2],
               "co": [0, 1, 2],

               "Trefb": 487.21, "rhorefb": 2.988659, "krefb": 0.66,
               "nb": [-0.729410e-2, 0.110497e-1, 0.416339e-1, -0.289236e-1,
                      -0.311487e-1, 0.278399e-1, 0.112682e-1, -0.121100e-1,
                      -0.138322e-2, 0.211196e-2],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.5e-9, "Tcref": 730.8}

    _thermal = thermo0,


if __name__ == "__main__":
    r113=R113(T=300, P=0.1)
    print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r113.T, r113.rho, r113.u.kJkg, r113.h.kJkg, r113.s.kJkgK, r113.cv.kJkgK, r113.cp.kJkgK, r113.w)
    print r113.k, r113.mu
