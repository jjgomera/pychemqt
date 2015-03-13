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

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [1.616, -0.4117, -0.792, 0.758, 1.217],
           "exp": [531., 751., 1989., 2484., 6859.],
           "ao_hyp": [], "hyp": []}

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4579856475, 1.888076782],
           "ao_exp": [1.616, -0.4117, -0.792, 0.758, 1.217],
           "titao": [16.0205159149, 22.6580178006, 60.0090511389,
                     74.9434303817, 206.9392065168]}

    CP2 = {"ao": 0.72480209e3,
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
        "__name__": "Helmholtz equation of state for normal hydrogen of Leachman et al. (2007).",
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
        "__doc__": u"""Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
        "R": 8.31434,
        "cp": CP2,

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
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP1,

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
        "__doc__":  u"""Bender, E. "Equation of state of normal hydrogen in the range 18 to 700 K and 1 to 500 bar," VDI Forschungsheft N 609, pp. 15-20, 1982.""",
        "R": 8.3143,
        "cp": CP1,

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

    _surface = {"sigma": [0.005369], "exp": [1.065]}
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
             "__name__": "McCarty (1972)"}

    _viscosity = visco0,

    def _visco0(self, rho, T, fase):
        """McCarty, R.D. and Weber, L.A., "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia," Natl. Bur. Stand., Tech. Note 617, 1972."""

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

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "McCarty (1972)"}

    _thermal = thermo0,

    def _thermo0(self, rho, T, fase):
        """McCarty, R.D. and Weber, L.A., "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia," Natl. Bur. Stand., Tech. Note 617, 1972."""
        # TODO:
        return unidades.ThermalConductivity(0)


if __name__ == "__main__":
    hidrogeno=H2(T=300, P=0.1)
    print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (hidrogeno.T, hidrogeno.rho, hidrogeno.u.kJkg, hidrogeno.h.kJkg, hidrogeno.s.kJkgK, hidrogeno.cv.kJkgK, hidrogeno.cp.kJkgK, hidrogeno.w)
    hidrogeno=H2(T=300, P=0.1, eq=2)
    print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (hidrogeno.T, hidrogeno.rho, hidrogeno.u.kJkg, hidrogeno.h.kJkg, hidrogeno.s.kJkgK, hidrogeno.cv.kJkgK, hidrogeno.cp.kJkgK, hidrogeno.w)
    hidrogeno=H2(T=300, P=0.1, eq=3)
    print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (hidrogeno.T, hidrogeno.rho, hidrogeno.u.kJkg, hidrogeno.h.kJkg, hidrogeno.s.kJkgK, hidrogeno.cv.kJkgK, hidrogeno.cp.kJkgK, hidrogeno.w)
