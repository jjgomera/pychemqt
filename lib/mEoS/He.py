#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp, log

from lib.meos import MEoS
from lib import unidades


class He(MEoS):
    """Multiparameter equation of state for helium"""
    name = "helium"
    CASNumber = "7440-59-7"
    formula = "He"
    synonym = "R-704"
    rhoc = unidades.Density(69.5800323874)
    Tc = unidades.Temperature(5.1953)
    Pc = unidades.Pressure(227.61, "kPa")
    M = 4.002602  # g/mol
    Tt = unidades.Temperature(2.1768)
    Tb = unidades.Temperature(4.2226)
    f_acent = -0.385
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 212

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [], "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [13.628409737, -143.470759602],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for helium of Ortiz-Vega et al. (2013).",
        "__doi__": {"autor": "Ortiz-Vega, D.O., Hall, K.R., Holste, J.C., Arp, V.D., and Lemmon, E.W.",
                    "title": "Interim equation", 
                    "ref": "final equation of state to be published in J. Phys. Chem. Ref. Data, 2013.",
                    "doi":  ""}, 
                    
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 1000000.0, "rhomax": 141.22, 
        "Pmin": 5.0335, "rhomin": 36.48, 

        "nr1": [0.014799269, 3.06281562, -4.25338698, 0.05192797,
                -0.165087335, 0.087236897],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1.0, 0.426, 0.631, 0.596, 1.705, 0.568],

        "nr2": [2.10653786, -0.62835030, -0.28200301, 1.04234019, -0.07620555,
                 -1.35006365],
        "d2": [1, 1, 3, 2, 2, 1],
        "t2": [0.9524, 1.471, 1.48, 1.393, 3.863, 0.803],
        "c2": [1, 2, 2, 1, 2, 1],
        "gamma2": [1]*6,

        "nr3": [0.11997252, 0.10724500, -0.35374839, 0.75348862, 0.00701871,
                0.226283167, -0.22464733, 0.12413584, 0.00901399],
        "d3": [1, 1, 1, 2, 2, 2, 3, 2, 2],
        "t3": [3.273, 0.66, 2.629, 1.4379, 3.317, 2.3676, 0.7545, 1.353,
                1.982],
        "alfa3": [8.674, 4.006, 8.1099, 0.1449, 0.1784, 2.432, 0.0414, 0.421,
                   5.8575],
        "beta3": [8.005, 1.15, 2.143, 0.147, 0.154, 0.701, 0.21, 0.134,
                   19.256],
        "gamma3": [1.1475, 1.7036, 1.6795, 0.9512, 4.475, 2.7284, 1.7167,
                    1.5237, 0.7649],
        "epsilon3": [0.912, 0.79, 0.90567, 5.1136, 3.6022, 0.6488, 4.2753,
                      2.744, 0.8736],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for helium of Ortiz-Vega et al. (2010).",
        "__doi__": {"autor": "Ortiz-Vega, D.O., Hall, K.R., Holste, J.C., Arp, V.D., and Lemmon, E.W.",
                    "title": "Interim equation", 
                    "ref": "final equation of state to be published in J. Phys. Chem. Ref. Data, 2013.",
                    "doi":  ""}, 
                    
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 1000000.0, "rhomax": 141.22, 
        "Pmin": 5.043, "rhomin": 36.46, 

        "nr1":  [0.9288766e-2, 0.9258069, -0.1718156e1, 0.7606137e0,
                 -0.1024864e1, 0.1052455],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1.0, 0.28, 0.735, 0.64, 0.82, 1.16],

        "nr2": [-0.1875722, -0.1287812, -0.2227619e-2, 0.1823465,
                -0.4450014e-1, -0.8729033e-4],
        "d2": [1, 1, 3, 2, 2, 8],
        "t2": [1.28, 2.0, 0.41, 1.33, 4.2, 0.6],
        "c2": [1, 2, 2, 1, 2, 1],
        "gamma2": [1]*6,

        "nr3": [0.3854320e-1, -0.9585106, -0.5454010e-1, -0.3687260e-1,
                -0.1021851e-2, 0.6166348e-1, 0.2493437e-1, -0.8127424e-2,
                -0.8233032e-2],
        "d3": [1, 1, 1, 2, 2, 2, 3, 3, 2],
        "t3": [3, 1, 8.2, 1, 2.71, 1, 1, 2, 1],
        "alfa3": [1.0833, 18.3824, 5.0573, 0.2832, 6.0582, 0.2444, 0.0539,
                  0.1850, 0.5941],
        "beta3": [0.0385, 19.8246, 9.3799, 0.8073, 0.0310, 0.0061, 0.3581,
                  0.7518, 7.4629],
        "gamma3": [1.9776, 1.6178, 0.4371, 0.5355, 0.7777, 0.4832, 0.8162,
                   1.2896, 0.3577],
        "epsilon3": [0.6914, 0.8590, 0.8787, 2.7182, 2.0301, 0.8900, 1.1790,
                     0.5680, 1.6412],
        "nr4": []}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Fundamental equation of state for helium of McCarty and Arp (1990).",
        "__doi__": {"autor": "McCarty, R.D. and Arp, V.D.",
                    "title": "A New Wide Range Equation of State for Helium", 
                    "ref": "Adv. Cryo. Eng., 35:1465-1475, 1990",
                    "doi":  ""}, 
                    
        "R": 8.31431,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 1500.0, "Pmax": 100000.0, "rhomax": 88.73, 
        "Pmin": 4.8565, "rhomin": 36.537, 

        "nr1": [-0.208984171567e1, 0.381792817549, -0.441393943069e-1,
                0.954038242224e-1, 0.115744872054e1, -0.287584069992e1,
                0.754294125268, -0.237177092854, 0.264743463330e-1,
                -0.164983375328, 0.764132237117, 0.312978947837,
                -0.107558759761e-2, 0.109732330796, -0.309354837550,
                0.823154284944e-3, 0.149309620852e-1, -0.150469153718e-1,
                -0.213800009686e-2, 0.198095303505e-3, 0.195115121471e-2,
                -0.320152846941e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8],
        "t1": [3, 4, 5, 0, 0.5, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 1, 2, 3, 2, 2, 3, 3],

        "nr2": [0.208984171567e1, -0.381792817549, 0.441393943069e-1,
                0.119594006419e1, -0.152740402594, 0.441393941765e-1,
                0.469369369369, -0.763702010715e-1, -0.206787489008e-2,
                0.744548107827e-1, 0.393354771579e-2, -0.689291627989e-3,
                0.601642971226e-2, 0.983386926042e-3, -0.235321870328e-3,
                0.201249794359e-3, 0.485142401906e-3, -0.470643739266e-4],
        "d2": [0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8, 8, 10, 10, 10],
        "t2": [3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5],
        "c2": [2]*18,
        "gamma2": [1]*18}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for helium of McCarty and Arp (1990).",
        "__doi__": {"autor": "McCarty, R.D. and Arp, V.D.",
                    "title": "A New Wide Range Equation of State for Helium", 
                    "ref": "Adv. Cryo. Eng., 35:1465-1475, 1990",
                    "doi":  ""}, 
                    
        "R": 8.31431,
        "cp": CP1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 1500.0, "Pmax": 100000.0, "rhomax": 88.73, 
        "Pmin": 4.8565, "rhomin": 36.537, 

        "b": [None, 0.4558980227431e-3, 0.1260692007853e-1, -0.7139657549318e-1,
              0.9728903861441e-1, -0.1589302471562, 0.1454229259623e-4,
              -0.4708238429298e-3, 0.1132915223587e-1, 0.2410763742104e-1,
              -0.5093547838381e-7, 0.2699726927900e-4, -0.3954146691114e-3,
              0.1551961438127e-7, 0.1050712335785e-6, -0.5501158366750e-6,
              -0.1037673478521e-8, 0.6446881346448e-11, 0.3298960057071e-9,
              -0.3555585738784e-11, -0.6885401367690e-1, 0.9166109232806e-1,
              -0.6544314242937e-4, -0.3315398880031e-3, -0.2067693644676e-6,
              0.3850153114958e-6, -0.1399040626999e-9, -0.1888462892389e-10,
              -0.4595138561035e-13, 0.6872567403738e-13, -0.6097223119177e-17,
              -0.7636186157005e-16, 0.3848665703556e-16]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for helium of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032-3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 1500.0, "Pmax": 100000.0, "rhomax": 88.73, 
        "Pmin": 4.8565, "rhomin": 36.537, 

        "nr1": [-0.45579024006737, 0.12516390754925e1, -0.15438231650621e1,
                0.20467489707221e-1],
        "d1": [1, 1, 1, 4],
        "t1": [0, 0.125, 0.75, 1.],

        "nr2": [-0.34476212380781, -0.20858459512787e-1, 0.16227414711778e-1,
                -0.57471818200892e-1, 0.19462416430715e-1, -0.33295680123020e-1 ,
                -0.10863577372367e-1, -0.22173365245954e-1],
        "d2": [1, 3, 5, 5, 5, 2, 1, 2],
        "t2": [0.75, 2.625, 0.125, 1.25, 2., 1., 4.5, 5.],
        "c2": [1, 1, 1, 1, 1, 2, 3, 3],
        "gamma2": [1]*8,

        "nr3": [],
        "nr4": []}

    eq = helmholtz1, helmholtz2, helmholtz3, MBWR, GERG
    _PR = -0.005886

    _surface = {"sigma": [0.0004656, 0.001889, -0.002006],
                "exp": [1.04, 2.468, 2.661]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [0.517254], "expt1": [0], "expd1": [1],
                   "a2": [-0.203, 0.039, 7.47],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 3]}
    _melting = {"eq": 1, "Tref": 1, "Pref": 1000,
                "Tmin": Tt, "Tmax": 1500.0,
                "a1": [-1.7455837, 1.6979793], "exp1": [0, 1.555414],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-3.8357, 1.7062, -0.71231, 1.0862],
        "exp": [1, 1.5, 1.25, 2.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-1.5789, -10.749, 17.711, -15.413, -14.352],
        "exp": [0.333, 1.5, 2.1, 2.7, 9.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.5789, -10.749, 17.711, -15.413, -14.352],
        "exp": [0.333, 1.5, 2.1, 2.7, 9.0]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "Arp,V.D (1998)", 
              "__doi__": {"autor": "Arp, V.D. McCarty, R.D., and Friend, D.G.",
                          "title": "Thermophysical properties of Helium-4 from 0.8 to 1500 K with pressures to 2000 Mpa", 
                          "ref": "NIST Technical Note 1334",
                          "doi": ""}}
    
    def _visco0(self, rho, T, fase=None, coef=False):
            
        Visco0 = lambda T: -0.135311743/log(T) + 1.00347841 + \
            1.20654649*log(T) - 0.149564551*log(T)**2 + 0.0125208416*log(T)**3

        def ViscoE(T, rho):
            x = log(T)
            B = -47.5295259/x+87.6799309-42.0741589*x+8.33128289*x**2-0.589252385*x**3
            C = 547.309267/x-904.870586+431.404928*x-81.4504854*x**2+5.37005733*x**3
            D = -1684.39324/x+3331.08630-1632.19172*x+308.804413*x**2-20.2936367*x**3
            return rho.gcc*B+rho.gcc**2*C+rho.gcc**3*D


        if T < 100:
            # Section 4.2.1 for 3.5 < T < 100
            no = Visco0(T)
            ne = ViscoE(T, rho)
            n = exp(no+ne)
        else:
            # Section 4.2.1 for T > 100
            no = 196*T**0.71938*exp(12.451/T-295.67/T**2-4.1249)
            ne = exp(Visco0(T)+ViscoE(T, rho))-exp(Visco0(T)+ViscoE(T, unidades.Density(0)))
            n = no+ne

        if coef:
            return ne
        else:
            return unidades.Viscosity(n*1e-6, "P")

    _viscosity = visco0,

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "Hands (1981)",
               "__doi__": {"autor": "Hands, B.A. and Arp, V.D.",
                           "title": "A Correlation of Thermal Conductivity Data for Helium", 
                           "ref": "Cryogenics, 21(12):697-703, 1981",
                           "doi": "10.1016/0011-2275(81)90211-3"}}

    def _thermo0(self, rho, T, fase=None):
        bkt = 0.0
        n0 = [3.739232544, -26.20316969, 59.82252246, -49.26397634]
        suma = sum([ni*self.T**(-i-1) for i, ni in enumerate(n0)])
        lg = 2.7870034e-3*self.T**0.7034007057*exp(suma)

        tau = self.T**(1.0/3.0)
        delta = self.rho/self.M/0.2498376
        t = [0, 3, 1, 2, 0, 1, 2, 0, 1, 2, -3, 0, 1, 2, -3]
        d = [1, 1, 1, 1, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2]
        c = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]
        n = [0.186297053e-3, -0.7275964435e-6, -0.1427549651e-3,
             0.3290833592e-4, -0.5213335363e-7, 0.4492659933e-7,
             -0.5924416513e-8, 0.7087321137e-5, -0.6013335678e-5,
             0.8067145814e-6, 0.3995125013e-6, -2.99050061466e-5,
             2.53733162271e-5, -3.40393839209e-6, -1.68574607754e-6]
        lk = 0
        for ti, di, ni, ci in zip(t, d, n, c):
            if ci == 0 or delta == 0:
                lk += ni*tau**ti*delta**di
            else:
                lk += ni*tau**ti*delta**di*log(delta**ci)

        lc = 0
        if 3.5 <= self.T <= 12:
            e1 = 2.8461
            e2 = .27156
            delta = 4.304
#            pcc = 227460.

            if self.rho > 0.0:
                bkt = 1.0/self.derivative("P", "rho", "T")/self.rho/1000.
            deld = abs((self.rho-69.158)/69.158)
            delt = abs((self.T-5.18992)/5.18992)
            r2 = (delt/0.2)**2+(deld/0.25)**2
            if r2 < 1.0 and self.rho > 0.0:
                xx = delt/deld**(1.0/.3554)
                x1 = (xx+.392)/.392
                x2b = x1**(2.0*.3554)
                x2be = (1.0+e2*x2b)**((1.1743-1.0)/2.0/.3554)
                hh = e1*x1*x2be
                dhdx = e1*x2be/.392+e1*e2/.392*x2b*x2be/(1.0+e2*x2b)*(1.1743-1.0)
                d2kt = (delta*hh-xx*dhdx/.3554)*deld**(delta-1.0)
                bkt1 = (69.158/self.rho)**2/d2kt/227460.
                bkt = r2*bkt+(1.0-r2)*bkt1

            eta = self._visco0(coef=True)
            bkcrit = 0
            if bkt >= 0 and self.rho > 0.0:
                bkcrit = 3.4685233e-17*self.T**2*bkt**0.5/self.rho/eta*self.dpdT**2*exp(-18.66*delt**2-4.25*deld**4)
            lc = 3.726229668*bkcrit

        return unidades.ThermalConductivity(lg+lk+lc)

    thermo1 = {"eq": 0,
               "method": "_thermo1",
               "__name__": "Peterser (1970)", 
               "__doi__": {"autor": "Peterser, H.",
                           "title": "The Properties of Helium: Density, Specific Heats. Viscosity, and Thermal Conductivity at Pressures from 1 to 100 bar and from Room Temperature to about 1800 K", 
                           "ref": "Denmark. Forskningscenter Risoe. Risoe-R; No. 224",
                           "doi": ""}}

    def _thermo1(self, rho, T, fase=None):
        k = 2.682e-3*(1+1.123e-3*self.P.bar)*T**(0.71*(1-2e-4*self.P.bar))
        return unidades.ThermalConductivity(k)

    _thermal = thermo0, thermo1
