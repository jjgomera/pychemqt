#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from math import exp, log

from scipy import roots
from scipy.constants import R

from lib.EoS.cubic import Cubic


# Table 2 from [2]_
dat = {
    98: (0.3536, 0.1054),
    105: (0.4058, 0.1803),
    # : (0.3593, 0.1441),  # phosphine
    208: (0.4303, 0.1196),
    # 951 : (0.4714, 0.1407),
    104: (0.4257, 0.1694),
    1: (0.1015, 0),
    62: (0.6042, 0.1832),
    50: (0.3731, 0.3096),
    63: (0.4907, 0.3435),
    46: (0.3844, 0.1683),
    110: (0.4892, 0.2851),
    107: (0.4478, 0),
    47: (0.3770, 0.1049),
    51: (0.5516, 0.3468),
    111: (0.6180, 1.2403),
    215: (0.5058, 0.1767),
    216: (0.4468, 0.1919),
    101: (0.5311, 0.1746),
    217: (0.4913, 0.2128),
    100: (0.4733, 0.2386),
    218: (0.5037, 0.1866),
    48: (0.3790, 0.2281),
    49: (0.5663, 0.2308),
    220: (0.5155, 0.2268),
    642: (0.5267, 0.2223),
    115: (0.4456, 0.1881),
    2: (0.3322, 0.1363),
    117: (0.7840, 0.1380),
    # : (0.4621, 0.2130),  # C1Mercaptan
    65: (0.5229, 0.2262),
    125: (0.2703, 0.4115),
    22: (0.4029, 0.1705),
    130: (0.5733, 0.1836),
    131: (0.5080, 0.2809),
    132: (0.4443, 0.2221),
    3: (0.4060, 0.1493),
    133: (0.4587, 0.2315),
    134: (0.9385, 0.2262),
    137: (0.4619, 0.2930),
    136: (0.4587, 0.2315),
    66: (0.5010, 0.1896),
    258: (0.5445, 0.1375),
    23: (0.4577, 0.1849),
    140: (0.4988, 0.2844),
    141: (0.5283, 0.4757),
    142: (0.5010, 0.1896),
    4: (0.4529, 0.1701),
    146: (0.9354, 0.3127),
    24: (0.4862, 0.2060),
    155: (0.5931, 0.3423),
    156: (0.5919, 0.3268),
    157: (0.5699, 0.5655),
    6: (0.5036, 0.1963),
    5: (0.4824, 0.1851),
    162: (0.5640, 0.2772),
    290: (0.4834, 0.7262),
    294: (0.5439, 0.3124),
    166: (0.6147, 0.3639),
    309: (0.6179, 0.3521),
    310: (0.6211, 0.3113),
    311: (0.6048, 0.3160),
    8: (0.5517, 0.1961),
    7: (0.5346, 0.2093),
    318: (0.6387, 0.1969),
    171: (0.4996, 0.2522),
    172: (0.5080, 0.2571),
    319: (0.5247, 0.2829),
    173: (0.5101, 0.2821),
    40: (0.4763, 0.2796),
    38: (0.4994, 0.3004),
    10: (0.5619, 0.2437),
    11: (0.5954, 0.2617),
    12: (0.6290, 0.2874),
    13: (0.6810, 0.3291),
    14: (0.7252, 0.3411),
    15: (0.7618, 0.3312),
    16: (0.7970, 0.3552),
    17: (0.8359, 0.3645),
    18: (0.8006, 0.5691),
    19: (0.9180, 0.3409),
    20: (0.9087, 0.3700),
    21: (0.7926, 0.5491),
    90: (0.7297, 0.6887),
    91: (0.9768, 0.5428),
    92: (1.0644, 0.4451)}


class TB(Cubic):
    r"""Trebble-Bishnoi cubic equation of state

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V^2+\left(b+c\right)V+bc+d^2}\\
        a = a_c\exp{q_1\left(1-T_r\right)}\\
        b = b_c\left(1+q_2\left(1-T_r+\ln{T_r}\right)\right)\\
        c = \frac{RT_c}{P_c}\left(1-3\zeta\right)\\
        q_1 = f\left(\omega, Z_c\right)\\
        q_2 = f\left(\omega\right)\\
        \zeta = 1.075Z_c\\
        d(cc/mol) = 0.241V_c(cc/mol)-5\\
        \end{array}

    This four parameter cubic equation include a b parameter with temperature
    dependence. The model need two compound specific parameters but include too
    a generalized correlation for both parameters defined in [1]_. The
    parameters are the updated values given in [2]_

    Examples
    --------
    Example 4.3 from [3]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = TB(300, 9.9742e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '90.8'
    >>> eq = TB(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vg.ccmol)
    '89.1'

    There are a tiny desviation, 2025 89.4 for two phases state and 88.3 for
    single phase state
    """

    __title__ = "Trebble-Bishnoi (1987)"
    __status__ = "TB"
    __doi__ = (
      {
        "autor": "Trebble, M.A., Bishnoi, P.R.",
        "title": "Development of a New Four-Parameter Cubic Equation of State",
        "ref": "Fluid Phase Equilibria 35 (1987) 1-8",
        "doi": "10.1016/0378-3812(87)80001-8"},
      {
        "autor": "Trebble, M.A.",
        "title": "Calculation of Constants in the Trebble-Bishnoi Equation of "
                 "State with an Extended Corresponding States Approach",
        "ref": "Fluid Phase Equilibria 45 (1989) 165-172",
        "doi": "10.1016/0378-3812(89)80255-9"},
      {
         "autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""})

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""
        ai = []
        bi = []
        ci = []
        di = []
        for cmp in self.componente:
            a, b, c, d = self.__lib(cmp, T)

            ai.append(a)
            bi.append(b)
            ci.append(c)
            di.append(d)

        self.ai = ai
        self.bi = bi
        self.ci = ci
        self.di = di

    def _GEOS(self, xi):
        coef = [self.ai, self.bi, self.ci, self.di]
        am, bm, cm, dm = self._mixture(None, xi, coef)

        delta = bm+cm
        epsilon = -bm*cm - dm**2
        return am, bm, delta, epsilon

    def __lib(self, cmp, T):
        Tr = T/cmp.Tc

        TcPci = cmp.Tc*cmp.Pc.MPa
        if cmp.f_acent < 0.225:
            # Eq 36
            TcPcH = 775.9 + 12003*cmp.f_acent - 57335*cmp.f_acent**2 + \
                91393*cmp.f_acent**3
        elif cmp.f_acent <= 1:
            TcPcH = 1876 - 1160*cmp.f_acent                         # Eq 37
        else:
            TcPcH = TcPci

        if cmp.f_acent >= -0.14:
            # Eq 34
            Zc = 0.29 - 0.0885*cmp.f_acent - 0.0005/(TcPci**0.5-TcPcH**0.5)
        else:
            Zc = 0.3024

        Xc = 1.075*Zc                                               # Eq 28
        d = 0.341*cmp.Vc.ccg*cmp.M - 0.005                          # Eq 29

        # Convert d from cc/mol to m3/mol
        d /= 1e6

        Dc = d*cmp.Pc/R/cmp.Tc                                      # Eq 12
        Cc = 1-3*Xc                                                 # Eq 6

        # Bc calculated ad the smallest positive root of Eq 8
        B = roots([1, 2-3*Xc, 3*Xc**2, -Dc**2-Xc**3])
        Bpositivos = []
        for Bi in B:
            if Bi > 0:
                Bpositivos.append(Bi)
        Bc = min(Bpositivos)

        Ac = 3*Xc**2 + 2*Bc*Cc + Bc + Cc + Bc**2 + Dc**2            # Eq 7

        if cmp.id in dat:
            q1, q2 = dat[cmp.id]
        else:
            if cmp.id == 212 and Tr <= 1:
                q1 = -0.31913
            elif cmp.f_acent < -0.1:
                # Eq 30
                q1 = 0.66208 + 4.63961*cmp.f_acent + 7.45183*cmp.f_acent**2
            elif cmp.f_acent <= 0.4:
                # Eq 32
                q1 = 0.35 + 0.7924*cmp.f_acent + 0.1875*cmp.f_acent**2 - \
                    28.93*(0.3-Zc)**2
            elif cmp.f_acent > 0.4:
                # Eq 33
                q1 = 0.32 + 0.9424*cmp.f_acent - 28.93*(0.3-Zc)**2

            if cmp.f_acent < -0.0423:
                q2 = 0
            elif cmp.f_acent <= 0.3:
                # Eq 26
                q2 = 0.05246 + 1.15058*cmp.f_acent - \
                    1.99348*cmp.f_acent**2 + 1.5949*cmp.f_acent**3 - \
                    1.39267*cmp.f_acent**4
            else:
                q2 = 0.17959 + 0.23471*cmp.f_acent                  # Eq 27

        alfa = exp(q1*(1-Tr))
        ac = Ac*R**2*cmp.Tc**2/cmp.Pc
        a = ac*alfa                                                 # Eq 23

        if T <= cmp.Tc:
            beta = 1 + q2*(1 - Tr + log(Tr))                        # Eq 19
        else:
            beta = 1
        bc = Bc*R*cmp.Tc/cmp.Pc
        b = bc*beta                                                 # Eq 18

        c = Cc*R*cmp.Tc/cmp.Pc

        return a, b, c, d


    def TB_Fugacidad(self, T, P):
        """Método de cálculo de la fugacidad mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        z=self.TB_Z(T, P)
        A=a*P/R_atml**2/T**2
        B=b*P/R_atml/T
        u=1+c/b
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        tita=abs(t)**0.5
        if t>=0:
            lamda=log((2*z+B*(u-tita))/(2*z+B*(u+tita)))
        else:
            lamda=2*arctan((2*z+u*B)/B/tita)-pi
        fi=z-1-log(z-B)+A/B/tita*lamda
        return unidades.Pressure(P*exp(fi), "atm")

    def TB_U_exceso(self, T, P):
        """Método de cálculo de la energía interna de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        v=self.TB_V(T, P)
        z=P*v/R_atml/T
        A=a*P/R_atml**2/T**2
        B=b*P/R_atml/T
        u=1+c/b
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        tita=abs(t)**0.5
        if t>=0:
            lamda=log((2*z+B*(u-tita))/(2*z+B*(u+tita)))
        else:
            lamda=2*arctan((2*z+u*B)/B/tita)-pi

        delta=v**2+(b+c)*v-b*c-d**2
        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        da=-q1*a/self.Tc
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
        else:
            db=0
        U=lamda/b/tita*(a-da*T)+db*T*(-R_atml*T/(v-b)+a/b**2/t*((v*(3*c+b)-b*c+c**2-2*d**2)/delta+(3*c+b)*lamda/b/tita)) #atm*l/mol
        return unidades.Enthalpy(U*101325/1000/self.peso_molecular, "Jkg")

    def TB_H_exceso(self, T, P):
        """Método de cálculo de la entalpía de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        p=unidades.Pressure(P, "atm")
        U=self.TB_U_exceso(T, P)
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        if t>=0:
            v=self.TB_V(T, P).m3g*self.peso_molecular
            return unidades.Enthalpy(p*v-R/T/self.peso_molecular+U.Jg, "Jg")
        else:
            Z=self.TB_Z(T, P)
            return unidades.Enthalpy(R/self.peso_molecular*T*(Z-1)+U.Jg, "Jg")

    def TB_Entalpia(self, T, P):
        """Método de cálculo de la entalpía mediante la ecuación de estado de Trebble-Bishnoi"""
        Ho=self._Ho(T)
        Delta=self.TB_H_exceso(T, P)
        return unidades.Enthalpy(Delta+Ho)

    def TB_S_exceso(self, T, P):
        """Método de cálculo de la entropía de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        H=self.TB_H_exceso(T, P)
        f=self.TB_Fugacidad(T, P)
        return unidades.SpecificHeat(H.Jg/T-R/self.peso_molecular*log(f.atm/P), "JgK")

    def TB_Entropia(self, T, P):
        """Método de cálculo de la entropía mediante la ecuación de estado de Trebble-Bishnoi"""
        So=self._so(T)
        Delta=self.TB_S_exceso(T, P)
        return unidades.SpecificHeat(Delta+So)

    def TB_Cv_exceso(self, T, P):
        """Método de cálculo de la capacidad calorífica a volumen constante de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        v=self.TB_V(T, P)
        z=P*v/R_atml/T
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        tita=abs(t)**0.5
        A=a*P/R_atml**2/T**2
        B=b*P/R_atml/T
        u=1+c/b
        delta=v**2+(b+c)*v-b*c-d**2
        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        da=-q1*a/self.Tc
        dda=q1**2*a/self.Tc**2
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
            ddb=-q2*b/beta/T**2
        else:
            db=0
            ddb=0

        dt=-db/b**2*(6*c+2*c**2/b+8*d**2/b)
        dtita=abs(dt)/20
        if t>=0:
            lamda=log((2*z+B*(u-tita))/(2*z+B*(u+tita)))
            dlamda=(db-db*tita-b*dtita)/(2*v+b+c-b*tita)-(db+db*tita+b*dtita)/((2*v+b+c+b*tita))
        else:
            lamda=2*arctan((2*z+u*B)/B/tita)-pi
            dlamda=2/(1+((2*v+b+c)/b/tita)**2)*(db/b/tita-(2*v+b+c)*(db/b**2/tita+dtita/b/tita**2))

        Cv=1/b/tita*(dlamda*(a-da*T)-lamda*dda*T-lamda*(a-da*T)*(db/b+dtita/tita))+(ddb*T+db)*(-R_atml*T/(v-b)+a/b**2/t*((v*(3*c+b)-b*c+c**2-2*d**2)/delta+(3*c+b)*lamda/b/tita))+db*T*(-R_atml/(v-b)-R_atml*T*db/(v-b)**2+1/b**2/t*(da-2*a*db/b-a*dt/t)*((v*(3*c+b)-b*c+c**2-2*d**2)/delta+(3*c+b)*lamda/b/tita)+a/b**2/t*(db*(v-c)*(v**2-2*c*v-c**2+d**2)/delta**2+db*lamda/b/tita+(3*c+b)/b/tita*(dlamda-lamda*(db/b+dtita/tita))))
        return unidades.SpecificHeat(Cv*101325/1000/self.peso_molecular, "JkgK")

    def TB_Cv(self, T, P):
        """Método de cálculo de la capacidad calorifica a volumen constante mediante la ecuación de estado de Trebble-Bishnoi"""
        Cvo=self.Cv_ideal(T)
        Delta=self.TB_Cv_exceso(T, P)
        return unidades.SpecificHeat(Delta+Cvo)

    def TB_Cp_exceso(self, T, P):
        """Método de cálculo de la capacidad calorífica a presión constante de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        v=self.TB_V(T, P)
        Cv=self.TB_Cv_exceso(T, P)

        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        delta=v**2+(b+c)*v-b*c-d**2
        da=-q1*a/self.Tc
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
        else:
            db=0

        dpdt=R_atml/(v-b)+R*T*db/(v-b)**2-da/delta+a*(v-c)*db/delta**2
        dpdv=-R_atml*T/(v-b)**2+a*(2*v+b+c)/delta**2
        Cp=-R-T*dpdt**2/dpdv*101325/1000+Cv.JgK*self.peso_molecular

        return unidades.SpecificHeat(Cp/self.peso_molecular, "JgK")

    def TB_Cp(self, T, P):
        """Método de cálculo de la capacidad calorifica a presión constante mediante la ecuación de estado de Trebble-Bishnoi"""
        Cpo=self.Cp_ideal(T)
        Delta=self.TB_Cp_exceso(T, P)
        return unidades.SpecificHeat(Delta+Cpo)

    def TB_Joule_Thomson(self, T, P):
        """Método de cálculo del coeficiente de Joule-Thomson mediante la ecuación de estado de Trebble-Bishnoi"""
        v=self.TB_V(T, P)
        Cp=self.TB_Cp(T, P)
        a, b, c, d, q1, q2=self.TB_lib(T, P)

        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        delta=v**2+(b+c)*v-b*c-d**2
        da=-q1*a/self.Tc
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
        else:
            db=0
        dpdt=R_atml/(v-b)+R*T*db/(v-b)**2-da/delta+a*(v-c)*db/delta**2
        dpdv=-R_atml*T/(v-b)**2+a*(2*v+b+c)/delta**2

        return -(T*dpdt+v*dpdv)/Cp/dpdv


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = TB(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = TB(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))

