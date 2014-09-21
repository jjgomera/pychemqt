#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp

from lib.meos import MEoS
from lib import unidades


class H2O(MEoS):
    """Multiparameter equation of state for water (including IAPWS95)

    >>> water=H2O(T=300, rho=996.5560)
    >>> print "%0.10f %0.8f %0.5f %0.9f" % (water.P.MPa, water.cv.kJkgK, water.w, water.s.kJkgK)
    0.0992418352 4.13018112 1501.51914 0.393062643

    >>> water=H2O(T=500, rho=0.435)
    >>> print "%0.10f %0.8f %0.5f %0.9f" % (water.P.MPa, water.cv.kJkgK, water.w, water.s.kJkgK)
    0.0999679423 1.50817541 548.31425 7.944882714

    >>> water=H2O(T=900., P=700)
    >>> print "%0.4f %0.8f %0.5f %0.8f" % (water.rho, water.cv.kJkgK, water.w, water.s.kJkgK)
    870.7690 2.66422350 2019.33608 4.17223802
    """
    name = "water"
    CASNumber = "7732-18-5"
    formula = "H2O"
    synonym = "R-718"
    Tc = unidades.Temperature(647.096)
    rhoc = unidades.Density(322.)
    Pc = unidades.Pressure(22064.0, "kPa")
    M = 18.015268  # g/mol
    Tt = unidades.Temperature(273.16)
    Tb = unidades.Temperature(373.1243)
    f_acent = 0.3443
    momentoDipolar = unidades.DipoleMoment(1.855, "Debye")
    id = 62

    CP1 = {"ao": 4.00632,
           "an": [], "pow": [],
           "ao_exp": [0.012436, 0.97315, 1.27950, 0.96956, 0.24873],
           "exp": [833, 2289, 5009, 5982, 17800],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 4.00392,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.01059, -0.98763, 3.06904, 0],
           "hyp": [0.415386589*Tc, 1.763895929*Tc, 3.874803739*Tc, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for water of Wagner and Pruß (2002).",
        "__doc__":  u"""Wagner, W., Pruß, A. The IAPWS formulation 1995 for the thermodyamic properties of ordinary water substance for general and scientific use. J. Phys. Chem. Ref. Data 31 (2002), 387 – 535.""",
        "R": 8.314371357587,
        "cp": CP1,

        "nr1": [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
                0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
                0.88089493102134e-2],
        "d1": [1, 1, 1, 2, 2, 3, 4],
        "t1": [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1],

        "nr2": [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
                -0.19232721156002, -0.25709043003438, 0.16074868486251,
                -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5,
                0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
                0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9,
                -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
                -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
                -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
                0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
                0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
                -0.20393486513704e-1, -0.16554050063734e-2, 0.19955571979541e-2,
                0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
                0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
                -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
                0.31777497330738, -0.11841182425981],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
               2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,
               6, 6],
        "d2": [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
               4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14,
               3, 6, 6, 6],
        "t2": [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10,
               10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,
               23, 10, 50, 44, 46, 50],
        "gamma2": [1]*44,

        "nr3": [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4],
        "d3": [3]*3,
        "t3": [0, 1, 4],
        "alfa3": [20]*3,
        "beta3": [150, 150, 250],
        "gamma3": [1.21, 1.21, 1.25],
        "epsilon3": [1.]*3,

        "nr4": [-0.14874640856724, 0.31806110878444],
        "a4": [3.5, 3.5],
        "b": [0.85, 0.95],
        "B": [0.2, 0.2],
        "C": [28, 32],
        "D": [700, 800],
        "A": [0.32, .32],
        "beta4": [0.3, 0.3]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for oxygen of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP2,

        "nr1": [0.82728408749586, -0.18602220416584e1, -0.11199009613744e1,
                0.15635753976056, 0.87375844859025, -0.36674403715731,
                0.53987893432436e-1],
        "d1": [1, 1, 1, 2, 2, 3, 4],
        "t1": [0.5, 1.25, 1.875, 0.125, 1.5, 1, 0.75],

        "nr2": [0.10957690214499e1, 0.53213037828563e-1, 0.13050533930825e-1,
                -0.41079520434476, 0.14637443344120, -0.55726838623719e-1,
                -0.11201774143800e-1, -0.66062758068099e-2, 0.46918522004538e-2],
        "c2": [1, 1, 1, 2, 2, 2, 3, 5, 5],
        "d2": [1, 5, 5, 1, 2, 4, 4, 1, 1],
        "t2": [1.5, 0.625, 2.625, 5, 4, 4.5, 3, 4, 6],
        "gamma2": [1]*9}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Saul and Wagner (1989).",
        "__doc__":  u"""Saul, A. and Wagner, W., "A Fundamental Equation for Water Covering the Range From the Melting Line to 1273 K at Pressures up to 25000 MPa," J. Phys. Chem. Ref. Data, 18(4):1537-1564, 1989.""",
        "R": 8.31434,
        "cp": CP1,

        "nr1": [0.2330009013, -0.1402091128e-1, 0.1172248041, -0.1850749499,
                0.1770110422, 0.5525151794e-1, -0.341325738e-3, 0.8557274367e-3,
                0.3716900685e-3, -0.1308871233e-3, 0.3216895199e-4,
                0.2785881034e-6],
        "d1": [1, 1, 2, 2, 2, 2, 3, 5, 5, 6, 7, 8],
        "t1": [0, 2, 0, 1, 2, 3, 5, 0, 1, 3, 2, 5],

        "nr2": [-0.352151113, 0.7881914536e-1, -0.151966661e-1, -0.1068458586,
                -0.2055046288, 0.9146198012, 0.3213343569e-3, -0.1133591391e1,
                -0.3107520749, 0.1217901527e1, -0.4481710831, 0.5494218772e-1,
                -0.8665222096e-4, 0.3844084088e-1, 0.9853044884e-2,
                -0.1767598472e-1, 0.1488549222e-2, -0.3070719069e-2,
                0.388080328e-2, -0.2627505215e-2, 0.5258371388e-3, -0.1716396901,
                0.7188823624e-1, 0.5881268357e-1, -0.145593888e-1, -0.12161394e-1],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
               2, 3, 3, 3, 3, 3],
        "d2": [1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 11, 11,
               11, 2, 2, 3, 3, 5],
        "t2": [5, 7, 9, 5, 4, 6, 13, 5, 2, 3, 2, 0, 11, 1, 4, 0, 0, 3, 5, 6,
               7, 13, 14, 13, 24, 15],
        "gamma2": [1]*26}

    eq = helmholtz1, GERG, helmholtz2
    _PR = 0.0043451

    _surface = {"sigma": [0.2358, -0.147375], "exp": [1.256, 2.256]}
#    _dielectric={"eq": 2, "Tref": 273.16, "rhoref": 1000.,
#                            "a0": [],  "expt0": [], "expd0": [],
#                            "a1": [], "expt1": [], "expd1": [],
#                            "a2": [], "expt2": [], "expd2": []}
#    _sublimation={"eq": 2, "Tref": 1, "Pref": 0.133332237, "a1": [-0.212144006e2, 0.273203819e2, -0.61059813e1], "exp1": [-0.9933333333, 0.206667, 0.703333], "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 6,
        "ao": [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719,
               1.80122502],
        "exp": [2, 3, 6, 7, 8, 15]}
    _liquid_Density = {
        "eq": 2,
        "ao": [1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352,
               -6.74694450e5],
        "exp": [1, 2, 5, 16, 43, 110]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581,
               -63.9201063],
        "exp": [1, 2, 4, 9, 18.5, 35.5]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (1997)"}

    visco1 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "Tref": 647.096, "muref": 1.0,
              "ek": 809.1, "sigma": 0.2641, "n_chapman": 0,
              "n_ideal": [151.138, -444.318, 398.262, -81.7008],
              "t_ideal": [0, 0.25, 0.5, 0.75],

              "a": [-1.17407105202836e-5, -3.78854818708520e-7, 3.56742875797909e-8],
              "b": [1.62216397984014e-6, -8.36595322447571e-6, 9.10862531286788e-8],
              "c": [1.92706925578893e-5, -1.28679815491711e-5, 0.0],
              "A": [-3.30144899918610e-10, 0.0, 1.02931444103415e-11],
              "B": [5.03139997945133e-10, 1.82304182380560e-10, 0.0],
              "C": [8.01449084635477e-10, 5.65613687804585e-9, 1.10163426018591e-10],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    def _visco0(self, coef=False):
        """International Association for the Properties of Water and Steam, "Revised Release on the IAPS Formulation 1985 for the Viscosity of Ordinary Water Substance," IAPWS, 1997
Kestin, J., Sengers, J.V., Kamgar-Parsi, B. and Levelt Sengers, J.M.H. "Thermophysical Properties of Fluid H2O," J. Phys. Chem. Ref. Data, 13(1):175-183, 1984."""
        Tr = self.T/647.226
        rhor = self.rho/self.M/17.6385386
        n0 = [1.0, 0.978197, 0.579829, -0.202354]
        fi0 = Tr**0.5/sum([n/Tr**i for i, n in enumerate(n0)])

        Li = [0, 1, 4, 5, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 3, 0, 3, 1, 3]
        Lj = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 6]
        Lij = [0.5132047, 0.3205656, -0.7782567, 0.1885447, 0.2151778,
               0.7317883, 1.241044, 1.476783, -0.2818107, -1.070786,
               -1.263184, 0.1778064, 0.4605040, 0.2340379, -0.4924179,
               -0.0417661, 0.1600435, -0.01578386, -0.003629481]
        lst = [lij*(1./Tr-1)**i*(rhor-1)**j for i, j, lij in zip(Li, Lj, Lij)]
        fi1 = exp(rhor*sum(lst))

        fi2 = 1.
        if 645.91 < self.T < 650.77 and 245.8 < self.rho < 405.3:
            x = rhor/self.derivative("P", "rho", "T")/self.M/17.6385386*22115
            if x >= 21.93:
                fi2 = 0.922*x**0.0263

        if coef:
            return fi0, fi1
        else:
            return unidades.Viscosity(55.071*fi0*fi1*fi2, "muPas")

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (1997)"}

    _thermal = thermo0,

    def _thermo0(self):
        """International Association for the Properties of Water and Steam, "Revised Release on the IAPS Formulation 1985 for the Thermal Conductivity of Ordinary Water Substance," IAPWS, 1998.
Kestin, J., Sengers, J.V., Kamgar-Parsi, B. and Levelt Sengers, J.M.H. "Thermophysical Properties of Fluid H2O," J. Phys. Chem. Ref. Data, 13(1):175-183, 1984."""
        rhor = self.rho/self.M/17.6385386
        Tr = self.T/647.226

        Lo = [1.0, 6.978267, 2.599096, -0.998254]
        L0 = Tr**0.5/sum([Li/Tr**i for i, Li in enumerate(Lo)])

        L1i = [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2]
        L1j = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5]
        L1ij = [1.3293046, 1.7018363, 5.2246158, 8.7127675, -1.8525999,
                -0.40452437, -2.2156845, -10.124111, -9.5000611, 0.9340469,
                0.2440949, 1.6511057, 4.9874687, 4.3786606, 0.018660751,
                -0.76736002, -0.27297694, -0.91783782, -0.12961068, 0.37283344,
                -0.43083393, 0.044809953, -0.11203160, 0.13333849]
        lst = [Lij*(1./Tr-1)**i*(rhor-1)**j for i, j, Lij in zip(L1i, L1j, L1ij)]
        L1 = exp(rhor*sum(lst))

        eta0, eta1 = self._visco0(True)
        x = rhor/self.derivative("P", "rho", "T")/self.M/17.6385386*22115
        L2 = 0.0013848*(self.T/rhor)**2/eta0/eta1*(self.dpdT/22115)**2*x**0.4678*rhor**0.5*exp(-18.66*(Tr-1)**2-(rhor-1)**4)
        return unidades.ThermalConductivity(.4945*(L0*L1+L2), "mWmK")

    def _Dielectric(self):
        """Equation for the Dielectric constant

        >>> "%.7f" % _Dielectric(999.242866, 298.15)
        '78.5907250'
        >>> "%.8f" % _Dielectric(26.0569558, 873.15)
        '1.12620970'
        """
        k = 1.380658e-23
        Na = 6.0221367e23
        alfa = 1.636e-40
        epsilon0 = 8.854187817e-12
        mu = 6.138e-30
        M = 0.018015268

        d = self.rho/self.rhoc
        Tr = self.Tc/self.T
        I = [1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10, None]
        J = [0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10, None]
        n = [0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396,
             -0.298217036956, -0.108863472196, 0.949327488264e-1,
             -0.980469816509e-2, 0.165167634970e-4, 0.937359795772e-4,
             -0.123179218720e-9, 0.196096504426e-2]

        g = 1+n[11]*d/(self.Tc/228/Tr-1)**1.2
        for i in range(11):
            g += n[i]*d**I[i]*Tr**J[i]
        A = Na*mu**2*self.rho*g/M/epsilon0/k/self.T
        B = Na*alfa*self.rho/3/M/epsilon0

        return unidades.Dimensionless((1+A+5*B+(9+2*A+18*B+A**2+10*A*B+9*B**2)**0.5)/4/(1-B))

    @classmethod
    def _Melting_Pressure(cls, T=None):
#        if not T:
#            T=self.T
#        if not P:
#            P=self.P
        if 251.165 <= T <= 273.16:
            Tref = cls.Tt
            Pref = 0.611657
            Tita = T/Tref
            a = [0.119539337e7, 0.808183159e5, 0.33382686e4]
            expo = [3., 0.2575e2, 0.10375e3]
            suma = 1
            for ai, expi in zip(a, expo):
                suma += ai*(1-Tita**expi)
            P = suma*Pref
        elif 251.165 < T <= 256.164:
            Tref = 251.165
            Pref = 208566.
            Tita = T/Tref
            P = Pref*(1-0.299948*(1-Tita**60.))
        elif 256.164 < T <= 273.31:
            Tref = 256.164
            Pref = 350100.
            Tita = T/Tref
            P = Pref*(1-1.18721*(1-Tita**8.))
        elif 273.31 < T <= 355:
            Tref = 273.31
            Pref = 632400.
            Tita = T/Tref
            P = Pref*(1-1.07476*(1-Tita**4.6))
        elif 355. < T <= 715:
            Tref = 355
            Pref = 2216000.
            Tita = T/Tref
            P = Pref*exp(1.73683*(1-1./Tita)-0.544606e-1*(1-Tita**5)+0.806106e-7*(1-Tita**22))

        return unidades.Pressure(P, "kPa")

    @classmethod
    def _Sublimation_Pressure(cls, T=None):
#        if not T:
#            T=self.T
        Tref = cls.Tt
        Pref = 0.611657
        Tita = T/Tref
        suma = 0
        a = [-0.212144006e2, 0.273203819e2, -0.61059813e1]
        expo = [0.333333333e-2, 1.20666667, 1.70333333]
        for ai, expi in zip(a, expo):
            suma += ai*Tita**expi
        return unidades.Pressure(exp(suma/Tita)*Pref, "kPa")


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=H2O(T=500., rho=838.025, recursion=False)
#    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)

#    print cyc5.cp.kJkgK, cyc5.cp0.kJkgK
