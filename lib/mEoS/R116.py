#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R116(MEoS):
    """Multiparameter equation of state for R116"""
    name = "hexafluoroethane"
    CASNumber = "76-16-4"
    formula = "CF3CF3"
    synonym = "R116"
    rhoc = unidades.Density(613.3245)
    Tc = unidades.Temperature(293.03)
    Pc = unidades.Pressure(3048.0, "kPa")
    M = 138.01182  # g/mol
    Tt = unidades.Temperature(173.1)
    Tb = unidades.Temperature(195.06)
    f_acent = 0.2566
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 236

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-10.7088650331, 8.9148979056],
           "ao_exp": [2.4818, 7.0622, 7.9951],
           "titao": [190/Tc, 622/Tc, 1470/Tc]}

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [2.4818, 7.0622, 7.9951],
           "exp": [190, 655, 1470],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 27.4009901,
           "an": [-2.6057376855e-6, 9.7501305219e-10, -6559.250418,
                  787904.9649, -34166787.86],
           "pow": [2, 3, -1.001, -2, -3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-116 of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=R116(T=295, rho=4*138.01182)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            295 4 3180.336 34509.528 161.389 120.218 2189.730 73.317
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 425.0, "Pmax": 50000.0, "rhomax": 12.31, 
        "Pmin": 26.1, "rhomin": 12.3, 

        "nr1": [1.1632, -2.8123, 0.77202, -0.14331, 0.10227, 0.00024629],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30893, -0.028499, -0.30343, -0.068793, -0.027218, 0.010665],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-116 of Kozlov (1996).",
        "__doc__":  u"""private communication with Dr. Alexander  D. Kozlov, Director, VNITs SMV Russian Research Center for Standartization Information and Certification of Materials, Nahimovsky prospect, 31, bld. 2 Moscow 117418, Russia. aldrkozlov@mail.ru""",
        "R": 8.31451,
        "cp": CP2,
        
        "Tmin": Tt, "Tmax": 425.0, "Pmax": 50000.0, "rhomax": 12.23, 
        "Pmin": 32.09, "rhomin": 12.231, 
        "Pmin": 32.09, "rhomin": 12.231, 

        "nr1": [2.1775273, -5.5052198, -1.3675742, -8.1284229e-1,
                -4.0207525e-1, 2.5890073, 1.4500537, -1.0445036, 9.8965288e-1,
                -8.6794888e-1, 2.8240917e-1, 4.5154220e-2, -3.0294024e-2,
                -1.7668398e-2, 2.0592774e-3],
        "d1": [1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 6, 7, 8],
        "t1": [0.25, 1, 3, 4, 0.25, 1, 3.5, 1.5, 2.5, 3, 3, 1, 3, 1, 1],

        "nr2": [4.2059839, 2.1500380e-1, -1.6449561e-1, -1.2396086e-1,
                1.5814552e-1, -1.4362345e-1, 1.8637877e-2, 1.6342835e-2],
        "d2": [1, 1, 4, 4, 5, 5, 8, 4],
        "t2": [2, 5, 2, 4, 8, 10, 10, 18],
        "c2": [1, 2, 2, 2, 3, 3, 3, 4],
        "gamma2": [1]*8}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.047593, -0.0073402], "exp": [1.2666, 1.9892]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73997e1, 0.22554e1, -0.23385e1, -0.35244e1, 0.40350],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.68490e2, -0.24772e3, 0.35824e3, -0.25290e3, 0.76880e2],
        "exp": [0.64, 0.79, 0.95, 1.14, 1.33]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.34135e1, -0.14529e3, 0.23651e3, -0.22276e3, 0.23103e3, -0.17433e3],
        "exp": [0.428, 2.0, 2.24, 3.0, 4.0, 5.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.355404, -0.464337, 0.257353e-1],
              "__name__": "Huber (2003)",
              "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., "Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a", Ind. Eng. Chem. Res., 42:3163-3178, 2003.""",
              "ek": 226.16, "sigma": 0.5249,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.2509/M**0.5,

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 226.16, "etaref_virial": 0.08709,

              "Tref_res": 513.9, "rhoref_res": 5.991*M, "etaref_res": 1000,
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

               "Tref": 1, "kref": 1.05,
               "no": [-1.05248e-2, 8.00982e-5],
               "co": [0, 1],

               "Trefb": 1.0, "rhorefb": 4.444, "krefb": 1.64e-3,
               "nb": [1.836526, 5.126143, -1.436883, 6.261441e-1],
               "tb": [0, 0, 0, 0],
               "db": [1, 2, 3, 4],
               "cb": [0, 0, 0, 0],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496e-1, "qd": 0.5e-9, "Tcref": 439.545}

    _thermal = thermo0,


if __name__ == "__main__":
    r116=R116(T=300, P=0.1)
    print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r116.T, r116.rho, r116.u.kJkg, r116.h.kJkg, r116.s.kJkgK, r116.cv.kJkgK, r116.cp.kJkgK, r116.w)
    print r116.k, r116.mu
