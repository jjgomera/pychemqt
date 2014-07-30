#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp

from lib.meos import MEoS
from lib import unidades

class D2O(MEoS):
    """Ecuación de estado de multiparametros para el agua pesada

    >>> water=D2O(T=300, rho=996.5560)
    >>> print "%0.10f %0.8f %0.5f %0.9f" % (water.P.MPa, water.cv.kJkgK, water.w, water.s.kJkgK)
    0.0992418352 4.13018112 1501.51914 0.393062643
    """
    name="heavy water"
    CASNumber="7789-20-0"
    formula="D2O"
    synonym="deuterium oxide"
    Tc=unidades.Temperature(643.89)
    rhoc=unidades.Density(357.992)
    Pc=unidades.Pressure(21671.0, "kPa")
    M=20.0275      #g/mol
    Tt=unidades.Temperature(276.97)
    Tb=unidades.Temperature(374.563)
    f_acent=0.364
    momentoDipolar=unidades.DipoleMoment(1.9, "Debye")
    id=62

    CP1={  "ao": 0.39176485e1,
                "an": [-0.31123915e-3, 0.41173363e-5, -0.28943955e-8, 0.63278791e-12, 0.78728740],
                "pow": [1.00, 2.00, 3.00, 4.00, -0.99],
                "ao_exp": [],
                "exp": [],
                "ao_hyp": [], "hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for heavy water of Hill et al. (1982).",
        "__doc__":  u"""Hill, P.G., MacMillan, R.D.C., and Lee, V., "A Fundamental Equation of State for Heavy Water," J. Phys. Chem. Ref. Data, 11(1):1-14, 1982.""",
        "R": 8.3143565,
        "cp": CP1,

        "nr1": [-0.384820628204e3, 0.108213047259e4, -0.110768260635e4, 0.164668954246e4 , -0.137959852228e4, 0.598964185629e3, -0.100451752702e3, 0.419192736351e3, -0.107279987867e4, 0.653852283544e3, -0.984305985655e3, 0.845444459339e3 , -0.376799930490e3, 0.644512590492e2, -0.214911115714e3, 0.531113962967e3, -0.135454224420e3, 0.202814416558e3, -0.178293865031e3, 0.818739394970e2, -0.143312594493e2, 0.651202383207e2, -0.171227351208e3, 0.100859921516e2, -0.144684680657e2, 0.128871134847e2, -0.610605957134e1, 0.109663804408e1, -0.115734899702e2, 0.374970075409e2, 0.897967147669, -0.527005883203e1, 0.438084681795e-1, 0.406772082680, -0.965258571044e-2, -0.119044600379e-1],
        "d1": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8],
        "t1": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 0, 1, 0, 1, 0, 1],

        "nr2": [0.382589102341e3, -0.106406466204e4, 0.105544952919e4, -0.157579942855e4, 0.132703387531e4, -0.579348879870e3, 0.974163902526e2, 0.286799294226e3, -0.127543020847e4, 0.275802674911e4, -0.381284331492e4, 0.293755152012e4, -0.117858249946e4, 0.186261198012e3],
        "c2": [1]*14,
        "d2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2],
        "t2": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6],
        "gamma2": [1.5394]*14}

    eq=helmholtz1,

    _surface={"sigma": [0.238, -0.152082], "exp": [1.25, 2.25]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.80236e1, 0.23957e1, -0.42639e2, 0.99569e2, -0.62135e2], "exp": [1.0, 1.5, 2.75, 3.0, 3.2]}
    _liquid_Density={ "eq": 1, "ao": [0.26406e1, 0.97090e1, -0.18058e2, 0.87202e1, -0.74487e1], "exp": [0.3678, 1.9, 2.2, 2.63, 7.3]}
    _vapor_Density={ "eq": 3, "ao": [-0.37651e1, -0.38673e2, 0.73024e2, -0.13251e3, 0.75235e2, -0.70412e2], "exp": [0.409, 1.766, 2.24, 3.04, 3.42, 6.9]}

    visco0={"eq": 0,
                    "method": "_visco0",
                    "__name__": "IAPWS (1994)"}

    _viscosity=visco0,

    def _visco0(self):
        """International Association for the Properties of Water and Steam, "Viscosity and Thermal Conductivity of Heavy Water Substance," Physical Chemistry of Aqueous Systems:  Proceedings of the 12th International Conference on the Properties of Water and Steam, Orlando, Florida, September 11-16, A107-A138, 1994"""
        Tr=self.T/643.89
        rhor=self.rho/self.M/17.87542

        fi0=Tr**0.5/sum([n/Tr**i for i, n in enumerate([1.0, 0.940695, 0.578377, -0.202044])])

        Li=[0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2, 3, 0, 1, 3, 5, 0, 1, 5, 3]
        Lj=[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6]
        Lij=[0.4864192, -0.2448372, -0.8702035, 0.8716056, -1.051126, 0.3458395, 0.3509007, 1.315436, 1.297752, 1.353448, -0.2847572, -1.037026, -1.287846, -0.02148229, 0.07013759, 0.4660127, 0.2292075, -0.4857462, 0.01641220, -0.02884911, 0.1607171, -0.009603846, -0.01163815, -0.008239587, 0.004559914, -0.003886659]
        fi1=exp(rhor*sum([Lij*(1./Tr-1)**i*(rhor-1)**j for i, j, Lij in zip(Li, Lj, Lij)]))

        return unidades.Viscosity(55.2651*fi0*fi1, "muPas")

    thermo0={"eq": 0,
                    "method": "_thermo0",
                    "__name__": "IAPWS (1994)"}

    _thermal=thermo0,

    def _thermo0(self):
        """International Association for the Properties of Water and Steam, "Viscosity and Thermal Conductivity of Heavy Water Substance," Physical Chemistry of Aqueous Systems:  Proceedings of the 12th International Conference on the Properties of Water and Steam, Orlando, Florida, September 11-16, A107-A138, 1994."""
        rhor=self.rho/self.M/17.87542
        Tr=self.T/643.89
        tau=Tr/(abs(Tr-1.1)+1.1)

        Lo=sum([Li*Tr**i for i, Li in enumerate([1.0, 37.3223, 22.5485, 13.0465, 0.0, -2.60735])])
        Lr=sum([Li*rhor**i for i, Li in enumerate([483.656, -191.039, 73.0358, -7.57467])])

        b=[-2.506, -167.310, 0.354296e5, 0.5e10, 0.144847, -5.64493, -2.8, -0.080738543, -17.943, 0.125698, -741.112]
        f1=exp(b[4]*Tr+b[5]*Tr**2)
        f2=exp(b[6]*(rhor-1)**2)+b[7]*exp(b[8]*(rhor-b[9])**2)
        f3=1+exp(60*(tau-1.)+20)
        f4=1+exp(100*(tau-1)+15)
        Lr+=b[1]*(1-exp(b[0]*rhor))
        Lc=b[2]*f1*f2*(1+f2**2*(b[3]*f1**4/f3+3.5*f2/f4))
        Ll=b[10]*f1**1.2*(1-exp(-(rhor/2.5)**10))

        return unidades.ThermalConductivity(0.742128e-3*(Lo*Lr+Lc+Ll))

if __name__ == "__main__":
    import doctest
    doctest.testmod()

