#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp, log

from lib.meos import MEoS
from lib import unidades


class H2(MEoS):
    """Multiparameter equation of state for hydrogen (normal)"""
    name = "hydrogen"
    CASNumber = "1333-74-0"
    formula = "H2"
    synonym = "R-702"
    rhoc = unidades.Density(31.26226704)
    Tc = unidades.Temperature(33.145)
    Pc = unidades.Pressure(1296.4, "kPa")
    M = 2.01588  # g/mol
    Tt = unidades.Temperature(13.957)
    Tb = unidades.Temperature(20.369)
    f_acent = -0.219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4579856475, 1.888076782],
           "ao_exp": [1.616, -0.4117, -0.792, 0.758, 1.217],
           "titao": [16.0205159149, 22.6580178006, 60.0090511389,
                     74.9434303817, 206.9392065168]}
                     
    Fi2 = {"ao_log": [1, 1.47906],
           "pow": [0, 1],
           "ao_pow": [13.796443393, -175.864487294],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [0.95806, 0.45444, 1.56039, -1.3756],
           "hyp": [6.891654113, 9.84763483, 49.76529075, 50.367279301]}

    CP1 = {"ao": 0.72480209e3,
           "an": [0.12155215e11, -0.36396763e10, 0.43375265e9, -0.23085817e8,
                  -0.38680927e4, 0.88240136e5, -0.78587085e4, -0.18426806e3,
                  0.21801550e2, -0.13051820e1, 0.21003175e-1, 0.23911604e-2,
                  -0.18240547e-3, 0.56149561e-5, -0.73803310e-7, 0.66357755e-11],
           "pow": [-7, -6, -5, -4, -3, -2, -1.001, 0.5, 1, 1.5, 2, 2.5, 3, 3.5,
                   4, 5],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for normal hydrogen of Leachman et al. (2009).",
        "__doi__": {"autor": "Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.",
                    "title": "Fundamental equations of state for parahydrogen, normal hydrogen, and orthohydrogen", 
                    "ref": "J. Phys. Chem. Ref. Data, 38 (2009), 721 – 748",
                    "doi": "10.1063/1.3160306"}, 
        "__test__": """
            >>> st=H2(T=13.957, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            13.957 7.3580 77.004 0.12985 −53.926 399.83 −3.0723 29.438 5.1616 6.2433 7.0212 10.564 1269.2 307.14
            """, # Table 14, Pag 746
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 2000000.0, "rhomax": 102.0, 
        "Pmin": 7.36, "rhomin": 38.2, 

        "nr1": [-6.93643, 0.01, 2.1101, 4.52059, 0.732564, -1.34086, 0.130985],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.6844, 1., 0.989, 0.489, 0.803, 1.1444, 1.409],

        "nr2": [-0.777414, 0.351944],
        "d2": [1, 3],
        "t2": [1.754, 1.311],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.0211716, 0.0226312, 0.032187, -0.0231752, 0.0557346],
        "d3": [2, 1, 3, 1, 1],
        "t3": [4.187, 5.646, 0.791, 7.249, 2.986],
        "alfa3": [1.685, 0.489, 0.103, 2.506, 1.607],
        "beta3": [0.171, 0.2245, 0.1304, 0.2785, 0.3967],
        "gamma3": [0.7164, 1.3444, 1.4517, 0.7204, 1.5445],
        "epsilon3": [1.506, 0.156, 1.736, 0.67, 1.662],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for hydrogen of Younglove (1982)",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen", 
                    "ref": "J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.",
                    "doi": ""}, 
        
        "R": 8.31434,
        "cp": CP1,
        "ref": "IIR", 

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 38.148, 
        "Pmin": 7.70, "rhomin": 38.3, 

        "b": [None, 0.4675528393416e-3, 0.4289274251454e-1, -0.5164085596504,
              0.2961790279801e1, -0.3027194968412e2, 0.1908100320379e-4,
              -0.1339776859288e-2, 0.3056473115421, 0.5161197159532e2,
              0.1999981550224e-6, 0.2896367059356e-3, -0.2257803939041e-1,
              -0.2287392761826e-5, 0.2446261478645e-4, -0.1718181601119e-2,
              -0.5465142603459e-6, 0.4051941401315e-8, 0.1157595123961e-5,
              -0.1269162728389e-7, -0.4983023605519e2, -0.1606676092098e3,
              -0.1926799185310, 0.9319894638928e1, -0.3222596554434e-3,
              0.1206839307669e-2, -0.3841588197470e-6, -0.4036157453608e-4,
              -0.1250868123513e-9, 0.1976107321888e-8, -0.2411883474011e-12,
              -0.4127551498251e-12, 0.8917972883610e-11]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032-3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 38.148, 
#        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.53579928451252e1, -0.62050252530595e1,  0.13830241327086,
                -0.71397954896129e-1,  0.15474053959733e-1],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.5, 0.625, 0.384, 0.625, 1.125],

        "nr2": [-0.14976806405771, -0.26368723988451e-1,  0.56681303156066e-1,
                -0.60063958030436e-1, -0.45043942027132,  0.42478840244500,
                -0.21997640827139e-1, -0.1049952137453e-1, -0.28955902866816e-2],
        "d2": [1, 5, 5, 5, 1, 1, 2, 5, 1],
        "t2": [2.625, 0, 0.25, 1.375, 4, 4.25, 5, 8, 8],
        "c2": [1, 1, 1, 1, 2, 2, 3, 3, 5],
        "gamma2": [1]*9,

        "nr3": [],
        "nr4": []}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen of Bender (1982).",
        "__doi__": {"autor": "Bender, E.",
                    "title": "Equation of state of normal hydrogen in the range 18 to 700 K and 1 to 500 bar", 
                    "ref": "VDI-Forschungsheft, no. 609, 1982, p. 15-20",
                    "doi": ""}, 
                    
        "R": 8.3143,
        "cp": CP1,
        "ref": "IIR", 

        "Tmin": 18.0, "Tmax": 700.0, "Pmax": 50000.0, "rhomax": 38.74, 
        "Pmin": 8.736, "rhomin": 38.7, 

        "nr1": [0.133442326203e1, -0.104116843433e1, 0.227202245707,
                0.300374270906, -0.463984214813, -0.178010492282e1,
                0.100460103605e1, -0.187200622541, 0.980276957749e-2,
                0.543224866339e-1, -0.263496312610e-1, 0.315432315759e-1,
                -0.525788294155e-1, -0.685380627808e-2, 0.344540276656e-1,
                -0.555747275982e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.133442326203e1, 0.104116843433e1, -0.227202245707,
                -0.378598758038, 0.249888797892, -0.498847982876e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.711139834571]*6,

        "nr3": [],
        "nr4": []}

    eq = helmholtz1, MBWR, GERG, helmholtz3
    _PR = -0.004803

    _surface = {"sigma": [-1.4165, 0.746383, 0.675625],
                "exp": [0.63882, 0.659804, 0.619149]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [2.0306, 0.0056], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [0.181, 0.021, -7.4],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 7.3578,
                "Tmin": Tt, "Tmax": 400.0,
                "a1": [1], "exp1": [0],
                "a2": [5626.3, 2717.2], "exp2": [1, 1.83],
                "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 7.7,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-8.065], "exp2": [0.93],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.489789e1, 0.988558, 0.349689, 0.499356],
        "exp": [1.0, 1.5, 2.0, 2.85]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.15456e2, -0.41720e2, 0.50276e2, -0.27947e2, 0.56718e1],
        "exp": [0.62, 0.83, 1.05, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.29962e1, -0.16724e2, 0.15819e2, -0.16852e2, 0.34586e2, -0.53754e2],
        "exp": [0.466, 2, 2.4, 4., 7., 8.]}

    visco0 = {"eq": 0,
             "method": "_visco0",
             "__name__": "Muzny (2013)", 
             "__doi__": {"autor": "Muzny, C.D., Huber, M.L., and Kazakov, A.F.",
                         "title": "Correlation for the Viscosity of Normal Hydrogen Obtained from Symbolic Regression", 
                         "ref": "J. Chem. Eng. Data, 2013, 58 (4), pp 969–979",
                         "doi": "10.1021/je301273j"}}

    def _visco0(self, rho, T, fase):
        sigma = 0.297
        ek = 30.41
        
        # Zero-Density Limit, Eq. 3-4
        T_ = T/ek
        ai = [2.0963e-1, -4.55274e-1, 1.43602e-1, -3.35325e-2, 2.76981e-3]
        suma = 0
        for i, a in enumerate(ai):
            suma += a*log(T_)**i
        S = exp(suma)
        no=0.021357*(self.M*T)**0.5/sigma/S
        
        # Excess Contribution, Eq. 5-7
        bi = [-0.187, 2.4871, 3.7151, -11.0972, 9.0965, -3.8292, 0.5166]
        B_ = 0
        for i, b in enumerate(bi):
            B_ += b/T**i
        B = B_*sigma**3
        n1=B*no
        
        # Simbolic Regression, Eq. 9
        rhor = rho/90.5
        Tr = T/self.Tc
        c = [6.43449673, 4.56334068e-2, 2.32797868e-1, 9.5832612e-1,
             1.27941189e-1, 3.63576595e-1]
        nc=c[0]*rhor**2*exp(c[1]*Tr+c[2]/Tr+c[3]*rhor**2/(c[4]+Tr)+c[5]*rhor**6)
        
        return unidades.Viscosity(no+n1+nc, "muPas")
        
    visco1 = {"eq": 0,
             "method": "_visco1",
             "__name__": "McCarty (1972)", 
             "__doi__": {"autor": "McCarty, R.D. and Weber, L.A.",
                         "title": "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 Psia", 
                         "ref": "NBS Technical Note 617",
                         "doi": ""}}

    def _visco1(self, rho, T, fase):
        DELV = lambda rho1, T1, rho2, T2: DILV(T1) + EXCESV(rho1, T2) \
            - DILV(T2)-EXCESV(rho2, T2)

        def EXVDIL(rho, T):
            A = exp(5.7694+log(rho.gcc)+0.65e2*rho.gcc**1.5-6e-6*exp(127.2*rho.gcc))
            B = 10+7.2*((rho.gcc/0.07)**6-(rho.gcc/0.07)**1.5)-17.63*exp(-58.75*(rho.gcc/0.07)**3)
            return A*exp(B/T)*0.1

        def DILV(T):
            b = [-0.1841091042788e2, 0.3185762039455e2, -0.2308233586574e2,
                 0.9129812714730e1, -0.2163626387630e1, 0.3175128582601,
                 -0.2773173035271e-1, 0.1347359367871e-2, -0.2775671778154e-4]
            suma = 0
            for i, b in enumerate(b):
                suma += b*T**((-3.+i)/3)
            return suma*100

        def EXCESV(rho, T):
            c = [-0.1324266117873e2, 0.1895048470537e2, 0.2184151514282e2,
                 0.9771827164811e5, -0.1157010275059e4, 0.1911147702539e3,
                 -0.3186427506942e4, 0.0705565000000]
            R2 = rho.gcc**0.5*(rho.gcc-c[7])/c[7]
            A = c[0]+c[1]*R2+c[2]*rho.gcc**0.1+c[3]*R2/T**2+c[4]*rho.gcc**0.1/T**1.5+c[5]/T+c[6]*R2/T
            B = c[0]+c[5]/T
            return 0.1*(exp(A)-exp(B))

        if T > 100:
            n = DILV(100) + EXVDIL(rho, 100) \
                + DELV(rho, T, rho, 100)
        else:
            n = DILV(T)+EXVDIL(rho, T)
        return unidades.Viscosity(n, "muPas")

    visco2 = {"eq": 4, "omega": 1,
              "__name__": "Quiñones-Cisneros (2011)",
              "__doi__": {"autor": "S.E.Quinones-Cisneros, M.L. Huber and U.K. Deiters",
                  "title": "model of 1-march-2011", 
                  "ref": "unpublished",
                  "doi": ""}, 

              "Tref": 33.145, "muref": 1.0,
              "ek": 59.7, "sigma": 0.2827, "n_chapman": 0,
              "n_ideal": [72.46400680522131e-1, -352.3929484813708e-1,
                          664.5332385860778e-1, -566.74979475607415e-1,
                          265.66570031561248e-1, -54.81307488054635e-1,
                          4.595978383724549e-1],
              "t_ideal": [0, 0.5, 0.75, 1, 1.25, 1.5],
              "n_poly": [1],
              "t_poly": [0.75],
              "n_polyden": [1, 1],
              "t_polyden": [0, 1],
              
              "nb": [1.0, -0.187, 75.6327, 3435.61, -312078, 7.77929e6,
                    -9.95841e7, 4.08557e8], 
              "tb": [0.0157768, 0, -1, -2, -3, -4, -5, -6], 

              "a": [-0.00002348389676311179e3, 0.00002197232806029717e3,
                    2.4547322430816313e-3, 3.9791170684039065e-8,
                    4.581319859008102e-3],
              "b": [0.000026869839733943842e3, 0.000027387647542474032e3,
                    0.000013065230652860072e3, 3.0723581102227345e-7,
                    -0.00007033089468735152e3],
              "c": [0, 0, 0, 0, 0],
              "A": [-3.912305916140789e-5, -2.1198288980972056e-6,
                    4.690087618888682e-6, 1.6938783854559677e-11,
                    9.39021777998824e-5],
              "B": [-6.381148168720446e-5, 5.178086941554603e-4, 
                    -4.5508093750991845e-5 -1.3780811004280076e-9
                    -3.7679840470735697e-4],
              "C": [0, 0, 0, 0, 0],
              "D": [4.3699367404316626e-7, 0.0, -1.1321685281996792e-8, 0, 0]}

    visco3 = {"eq": 1, "omega": 1,
              "__name__": "Vargaftik (1996)",
              "__doi__": {"autor": "Vargaftik, N.B., Vinogradov, Y.K. and Yargin, V.S.",
                          "title": "Handbook of Physical Properties of Liquids and Gases", 
                          "ref": "Hemisphere Publishing Corporation,New York, NY",
                          "doi": ""}, 

              "ek": 59.7, "sigma": 0.2827,
              "Tref": 32.938, "rhoref": 1.*M,
              "n_virial": [-2.1505e-1, 10.727e-1, -16.935e-1, 0.0, 22.702e-1,
                           2.2123e-1, 0.34163e-1, -0.043206e-1],
              "t_virial": [-1.5, -1, -0.5, 0, 0.5, 1.5, 2.],
              "Tref_virial": 32.938, "etaref_virial": 1.*M,

              "Tref_res": 32.938, "rhoref_res": 15.556*M, "etaref_res": 1.,
              "n_packed": [], "t_packed": [], 
              "n_poly": [-9.22703e-1, 6.41602, -5.98018, 2.89715e-1, 2.36429,
                         -2.78870e-1, -1.10595e1, 1.11582e1, 7.18928,
                         -7.76971, -1.21827,  1.47193],
              "t_poly": [0, -1, -2, -3, 0, 0, -1, -2, -1, -2, -1, -2],
              "d_poly": [1, 1, 1, 1, 2, 3, 3, 3, 4, 4, 5, 5],
              "g_poly": [0]*12,
              "c_poly": [0]*12,
              "n_num": [], "t_num": [], "d_num": [], "g_num": [], "c_num": [],
              "n_den": [], "t_den": [], "d_den": [], "g_den": [], "c_den": []}

    _viscosity = visco0, visco1, visco2, visco3

    thermo0 = {"eq": 1,
               "__name__": "Assael (2011)",
               "__doi__": {"autor": " Assael, M.J., Assael. J.-A.M., Huber, M.L., Perkins, R.A. and Takata, Y.",
                           "title": "Correlation of the Thermal Conductivity of Normal and Parahydrogen from the Triple Point to 1000 K and up to 100 MPa", 
                           "ref": "J. Phys. Chem. Ref. Data 40, 033101 (2011)",
                           "doi": "10.1063/1.3606499"}, 
               "__test__": """
                   >>> st=H2(T=298.15, rho=0)
                   >>> print "%0.5g" % st.k.mWmK
                   185.67
                   >>> st=H2(T=298.15, rho=0.80844)
                   >>> print "%0.5g" % st.k.mWmK
                   186.97
                   >>> st=H2(T=298.15, rho=14.4813)
                   >>> print "%0.5g" % st.k.mWmK
                   201.35
                   >>> st=H2(T=35, rho=0)
                   >>> print "%0.5g" % st.k.mWmK
                   26.988
                   >>> st=H2(T=35, rho=30)
                   >>> print "%0.5g" % st.k.mWmK
                   75.594
                   >>> st=H2(T=35, rho=30)
                   >>> print "%0.5g" % st.k.mWmK
                   71.854
                   >>> st=H2(T=18, rho=0)
                   >>> print "%0.5g" % st.k.mWmK
                   13.875
                   >>> st=H2(T=18, rho=75)
                   >>> print "%0.5g" % st.k.mWmK
                   104.48
                   """, # Table 4, Pag 8

               "Tref": 1.0, "kref": 1e-3,
               "no": [-1.24159e7, 5.04056e6, -4.80868e4, 3.26394e2,
                      9.56218e-2, 1.73488e-4, -3.12802e-8],
               "co": [0, 1, 2, 3, 4, 5, 6],
               "noden": [5.04305e6, -2.43753e4, 1.51523e2, 1.0], 
               "coden": [0, 1, 2, 3], 
 
               "Trefb": 33.145, "rhorefb": 15.508, "krefb": 1.,
               "nb": [.363081e-1, -.207629e-1, .31481e-1, -.143097e-1,
                      .17498e-2, .18337e-2, -.886716e-2, .15826e-1,
                      -.106283e-1, .280673e-2],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 0.4e-9, "Tcref": 49.7175}
               
    _thermal = thermo0,
