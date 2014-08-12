#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp

from lib.meos import MEoS
from lib import unidades


class pH2(MEoS):
    """Multiparamenter equation of state for hydrogen (para)

    >>> hidrogeno=pH2(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.5f %0.4f %0.4f %0.2f" % (hidrogeno.T, hidrogeno.rho, hidrogeno.h.kJkg, hidrogeno.s.kJkgK, hidrogeno.cv.kJkgK, hidrogeno.cp.kJkgK, hidrogeno.w)
    300.0 0.08077 27.469 0.14617 10.7187 14.8450 1309.82
    """
    name = "parahydrogen"
    CASNumber = "1333-74-0p"
    formula = "H2"
    synonym = "R-702p"
    rhoc = unidades.Density(31.32274344)
    Tc = unidades.Temperature(32.938)
    Pc = unidades.Pressure(1285.8, "kPa")
    M = 2.01588  # g/mol
    Tt = unidades.Temperature(13.8033)
    Tb = unidades.Temperature(20.271)
    f_acent = -0.219
    momentoipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973,
                      0.536078],
           "exp": [499, 826.5, 970.8, 1166.2, 1341.4, 5395, 10185],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 2.4995169,
           "an": [-0.11125185e-2, 0.27491461e-3, -0.10005269e-4, 0.22695404e-8,
                  -0.21031029e-12],
           "pow": [1, 1.5, 2, 3, 4],
           "ao_exp": [0.12353388e2, -0.17777676e2, 0.64309174e1, 0.73347521e1],
           "exp": [598, 778, 1101, 6207],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for parahydrogen of Leachman et al. (2007)",
        "__doc__":  u"""Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W. Fundamental equations of state for parahydrogen, normal hydrogen, and orthohydrogen. J. Phys. Chem. Ref. Data, 38 (2009), 721 – 748.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [-7.33375, 0.01, 2.60375, 4.66279, 0.682390, -1.47078, 0.135801],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.6855, 1, 1, 0.489, 0.774, 1.133, 1.386],

        "nr2": [-1.05327, 0.328239],
        "d2": [1, 3],
        "t2": [1.619, 1.162],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.119510],
        "d3": [2, 1, 3, 1, 1],
        "t3": [3.96, 5.276, 0.99, 6.791, 3.19],
        "alfa3": [1.7437, 0.5516, 0.0634, 2.1341, 1.777],
        "beta3": [0.194, 0.2019, 0.0301, 0.2383, 0.3253],
        "gamma3": [0.8048, 1.5248, 0.6648, 0.6832, 1.493],
        "epsilon3": [1.5487, 0.1785, 1.28, 0.6319, 1.7104],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for parahydrogen of Younglove (1982).",
        "__doc__": u"""Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
        "R": 8.31434,
        "cp": CP2,

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

    eq = helmholtz1, MBWR

    _surface = {"sigma": [0.005328], "exp": [1.065]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [2.0297, 0.0069], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [0.181, 0.021, -7.4],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 3]}
    _sublimation = {"eq": 2, "Tref": 1, "Pref": 0.133332237,
                    "a1": [4.009857354, -90.77568949], "exp1": [0, -1],
                    "a2": [], "exp2": [], "a3": [2.48983094], "exp3": [1]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.487767e1, 0.103359e1, 0.826680, -0.129412],
        "exp": [1.0, 1.5, 2.65, 7.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.13509, 0.40739e1, -0.53985e1, 0.55230e1, -0.23643e1],
        "exp": [0.15, 0.44, 0.7, 0.99, 1.31]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.57545e1, 0.38153e1, -0.12293e2, 0.15095e2, -0.17295e2, -0.34190e2],
        "exp": [0.53, 0.7, 1.7, 2.4, 3.3, 10]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "McCarty (1972)"}

    _viscosity = visco0,

    def _visco0(self):
        """McCarty, R.D. and Weber, L.A., "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia," Natl. Bur. Stand., Tech. Note 617, 1972."""

        DELV = lambda rho1, T1, rho2, T2: DILV(T1)+EXCESV(rho1, T2)-DILV(T2)-EXCESV(rho2, T2)

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
            c = [-0.1099981128e2, 0.1895876508e2, -0.3813005056e3,
                 0.5950473265e2, 0.1099399458e1, 0.8987269839e1,
                 0.1231422148e4, 0.311]
            R2 = rho.gcc**0.5*(rho.gcc-c[7])/c[7]
            A = c[0]+c[1]*R2+c[2]*rho.gcc**0.1+c[3]*R2/T**2+c[4]*rho.gcc**0.1/T**1.5+c[5]/T+c[6]*R2/T
            B = c[0]+c[5]/T
            return 0.1*(exp(A)-exp(B))

        if self.T > 100:
            n = DILV(100)+EXVDIL(self.rho, 100)+DELV(self.rho, self.T, self.rho, 100)
        else:
            n = DILV(self.T)+EXVDIL(self.rho, self.T)
        return unidades.Viscosity(n, "muPas")

    @classmethod
    def _Melting_Pressure(cls, T=None):
        Tref = 1
        Pref = 1000
        Tita = T/Tref
        a = [-0.265289115e2, 0.248578596, -0.21272389e2, 0.125746643]
        expo = [0, 0.1764739e1, 0, 0.1955e1]
        if T > 22:
            suma = -0.265289115e2+0.248578596*Tita**0.1764739e1
        else:
            suma = -0.21272389e2+0.125746643*Tita**0.1955e1

        return unidades.Pressure(suma*Pref, "kPa")

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "McCarty (1972)"}

    _thermal = thermo0,

    def _thermo0(self):
        """McCarty, R.D. and Weber, L.A., "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia," Natl. Bur. Stand., Tech. Note 617, 1972."""
        # TODO:
        return unidades.ThermalConductivity(0)
