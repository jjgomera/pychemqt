#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from math import exp, log

from scipy.optimize import fsolve
from scipy import constants as k

from lib import unidades
from lib.eos import EoS
from lib.EoS.Cubic import SRK
from lib.physics import R_atml
from lib.bip import Kij
from lib.mezcla import RhoL_CostaldMix


# Component paramters given in Starling, [1]_
# The order of parameters are:
# Bo, Ao, Co, gamma, b, a, alpha, c, Do, d, Eo
# Units of this costants are in imperial system
dat = {
    2: (0.723251, 7520.29, 2.71092e8, 1.48640, 0.925404, 2574.89, 0.468828,
        4.37222e8, 1.07737e10, 4.74891e4, 3.01122e10),
    3: (0.826059, 13439.3, 29.5195e8, 2.99656, 3.11206, 22404.5, 0.909681,
        68.1826e8, 25.747e10, 70.2189e4, 1468.19e10),
    4: (0.964762, 18634.7, 79.6178e8, 4.56182, 5.46248, 40066.4, 2.01402,
        274.461e8, 45.3708e10, 1505.20e4, 2560.53e10),
    5: (1.87890, 37264.0, 101.413e8, 7.11486, 8.58663, 47990.7, 4.23987,
        406.763e8, 85.3176e10, 2168.63e4, 8408.60e10),
    6: (1.56588, 32544.7, 137.436e8, 7.54122, 9.14066, 71181.8, 4.00985,
        700.044e8, 33.3159e10, 3642.38e4, 230.902e10),
    7: (1.27752, 35742.0, 228.430e8, 11.7384, 19.8384, 204344., 6.16154,
        1290.83e8, 142.115e10, 3492.20e4, 2413.26e10),
    8: (2.44417, 51108.2, 223.931e8, 11.8593, 16.6070, 162185., 7.06702,
        1352.86e8, 101.769e10, 3885.21e4, 3908.60e10),
    10: (2.66233, 45333.1, 526.067e8, 14.8720, 29.4983, 434517., 9.70230,
         3184.12e8, 552.158e10, 3274.60e4, 62643.3e10),
    11: (3.60493, 77826.9, 615.662e8, 24.7604, 27.4415, 359087, 21.8782,
         3748.76e8, 777.123e10, 835.115e4, 636.251e10),
    12: (4.86965, 81690.6, 996.546e8, 21.9888, 10.5907, 131646., 34.5124,
         6420.53e8, 790.575e10, 18590.6e4, 3464.19e10),
    22: (0.747945, 12133.9, 16.3203e8, 2.27971, 2.62914, 15988.1, 0.589158,
         40.9725e8, 5.17563e10, 90.3550e4, 1.61706e10),
    23: (0.114457, 6051.36, 97.4762e8, 4.07919, 7.64114, 81880.4, 1.36532,
         294.141e8, 70.5921e10, 541.935e4, 3412.50e10),
    46: (0.677022, 4185.05, 1.37936e8, 1.10011, 0.833470, 1404.59, 0.302696,
         0.844317e8, 1.95183e10, 3.11894e4, 121.648e10),
    49: (0.394117, 6592.03, 29.5902e8, 1.64916, 0.97443, 5632.85, 0.395525,
         27.4668e8, 40.9151e10, 5.99297e4, 1.02898e10),
    50: (0.297508, 10586.3, 21.1496e8, 1.20447, 2.53315, 2.05110, 0.165961,
         43.6132e8, 4.86518e10, 1.99731e4, 3.93226e10)
         }


class BWRS(EoS):
    r"""Benedict-Webb-Rubin equation of state modified by Starling [1]_

    .. math::
        \begin{align*}
        p = \rho RT + \left(B_0RT-A_0-\frac{C_0}{T^2}+\frac{D_0}{T^3}-
        \frac{C_0}{T^2}\right)\rho^2 +  \left(bRT-a-\frac{d}{T}\right)\rho^3 \\
        {} + \alpha a\left(a+\frac{d}{T}\right)\rho^6 + \frac{c\rho^3}{T^2}
        \left(1+\gamma \rho^2\right) \exp\left(-\gamma\rho^2\right)
        \end{align*}

    This equation use 11 compound specified parameters, available for 15 common
    subsgtances, but there are generalized correlation of these parameters as
    function of critic tempereture, critic density and acentric factor, given
    in [2]_.

    .. math::
        \begin{array}[t]{l}
        \rho_c B_0 = A_1 + B_1\omega\\
        \rho_c \frac{A_0}{RT_c} = A_2 + B_2\omega\\
        \rho_c \frac{C_0}{RT_c^3} = A_3 + B_3\omega\\
        \rho_c^2 \gamma = A_4 + B_4\omega\\
        \rho_c^2 b = A_5 + B_5\omega\\
        \rho_c^2 \frac{a}{RT_c} = A_6 + B_6\omega\\
        \rho_c^3 \alpha = A_7 + B_7\omega\\
        \rho_c^2 \frac{c}{RT_c^3} = A_8 + B_8\omega\\
        \rho_c \frac{D_0}{RT_c^4} = A_9 + B_9\omega\\
        \rho_c^2 \frac{d}{RT_c^2} = A_{10} + B_{10}\omega\\
        \rho_c \frac{E_0}{RT_c^3} = A_{11} + B_{11}\omega
        \exp\left(-3.8\omega\right)\\
        \end{array}

    The mixing rules to get the mixture parameters are:

    .. math::
        \begin{array}[t]{l}
        A_0 = \sum_i \sum_j x_i x_j A_{0i}^{1/2} A_{0j}^{1/2}
        \left(1-k_ij\right)\\
        B_0 = \sum_i x_i B_{0i}\\
        C_0 = \sum_i \sum_j x_i x_j C_{0i}^{1/2} C_{0j}^{1/2}
        \left(1-k_ij\right)^3\\
        D_0 = \sum_i \sum_j x_i x_j D_{0i}^{1/2} D_{0j}^{1/2}
        \left(1-k_ij\right)^4\\
        E_0 = \sum_i \sum_j x_i x_j E_{0i}^{1/2} E_{0j}^{1/2}
        \left(1-k_ij\right)^5\\
        \alpha = \left(\sum_i x_i \alpha_i^{1/3}\right)^3\\
        \gamma = \left(\sum_i x_i \gamma_i^{1/2}\right)^2\\
        a = \left(\sum_i x_i a_i^{1/3}\right)^3\\
        b = \left(\sum_i x_i b_i^{1/3}\right)^3\\
        c = \left(\sum_i x_i c_i^{1/3}\right)^3\\
        d = \left(\sum_i x_i d_i^{1/3}\right)^3\\
        \end{array}

    This model is very accurate for VLE for light normal hydrocarbons.

    The model include too the low reduced temperatures extensión from [2]_ with
    four aditional compound specific parameters

    .. math::
        \begin{align*}
        p = \rho RT + \left(B_0RT-A_0-\frac{C_0}{T^2}+\frac{D_0}{T^3}-
        \frac{C_0}{T^2}\right)\rho^2 +
        \left(bRT-a-\frac{d}{T}-\frac{e}{T^4}-\frac{f}{T^{23}}\right)\rho^3 \\
        {} + \alpha \left(a+\frac{d}{T}+\frac{e}{T^4}+\frac{f}{T^{23}}\right)
        \rho^6 + \left(\frac{c}{T^2}+\frac{g}{T^8}+\frac{h}{T^{17}}\right)
        \rho^3 \left(1+\gamma \rho^2\right) \exp\left(-\gamma\rho^2\right)
        \end{align*}

    with generalized correlation of this new parameters and with this mixing
    rules

    .. math::
        \begin{array}[t]{l}
        e = \left(\sum_i x_i e_i^{1/3}\right)^3\\
        f = \left(\sum_i x_i f_i^{1/3}\right)^3\\
        g = \sum_i x_i g_i\\
        h = \sum_i x_i h_i\\
        \end{array}

    """
    __title__ = "Benedict-Webb-Rubin-Starling (1973)"
    __status__ = "BWRS"

    __doi__ = (
      {
        "autor": "Starling, K.E.",
        "title": "Fluid Thermodynamics Properties of Light Petroleum Systems",
        "ref": "Gulf Publishing Company, 1973",
        "doi": ""},
      {
        "autor": "Han, M.S., Starling, K.E.",
        "title": "Thermo Data Refined for LPG. Part 14. Mixtures",
        "ref": "Hydrocarbon Processing 51(5) (1972) 129",
        "doi": ""},
      {
        "autor": "Nishiumi, H., Saito, S.",
        "title": "An Improved Generalized BWR Equation of State Applicable to "
                 "Low Reduced Temperatures",
        "ref": "J. Chem. Eng. Japan 8(5) (1975) 356-360",
        "doi": "10.1252/jcej.8.356"})

    def __init__(self, T, P, mezcla, **kwargs):
        """
        For use the extended version use the "extended" parameters with
        value True
        """
        if "extended" not in kwargs:
            kwargs["extended"] = False
        EoS.__init__(self, T, P, mezcla, **kwargs)

        self.kij = Kij(mezcla.ids, "BWRS")

        # Apply generalized correlation for calculate compound parameters
        Aoi = []
        Boi = []
        Coi = []
        Doi = []
        Eoi = []
        ai = []
        bi = []
        ci = []
        di = []
        alfai = []
        gammai = []
        ei = []
        fi = []
        gi = []
        hi = []
        for cmp in self.componente:
            kw = self._lib(cmp)

            Aoi.append(kw["Ao"])
            Boi.append(kw["Bo"])
            Coi.append(kw["Co"])
            Doi.append(kw["Do"])
            Eoi.append(kw["Eo"])
            ai.append(kw["a"])
            bi.append(kw["b"])
            ci.append(kw["c"])
            di.append(kw["d"])
            alfai.append(kw["alpha"])
            gammai.append(kw["gamma"])
            ei.append(kw["e"])
            fi.append(kw["f"])
            gi.append(kw["g"])
            hi.append(kw["h"])

        self.Aoi = Aoi
        self.Boi = Boi
        self.Coi = Coi
        self.Doi = Doi
        self.Eoi = Eoi
        self.ai = ai
        self.bi = bi
        self.ci = ci
        self.di = di
        self.alfai = alfai
        self.gammai = gammai
        self.ei = ei
        self.fi = fi
        self.gi = gi
        self.hi = hi

        self.x, self.Zl, self.Zg, self.xi, self.yi, self.Ki = self._Flash()

        if self.Zl:
            self.Vl = unidades.MolarVolume(self.Zl*k.R*T/P, "m3mol")   # l/mol
            self.rhoL = unidades.MolarDensity(1/self.Vl)

            kwl = self._mix(self.xi)
            self.HexcL = self._Hexc(self.Zg, T, **kwl)
        else:
            self.Vl = None
            self.rhoL = None
            self.HexcL = None

        if self.Zg:
            self.Vg = unidades.MolarVolume(self.Zg*k.R*T/P, "m3mol")  # l/mol
            self.rhoG = unidades.MolarDensity(1/self.Vg)

            kwg = self._mix(self.yi)
            self.HexcG = self._Hexc(self.Zg, T, **kwg)
        else:
            self.Vg = None
            self.rhoG = None
            self.HexcG = None

        # print("H_exc = ", self.H_exc)
        # print("Z = ", self.Z)
        # print('ρm = ', rhom)
        # print('ρ = ', rhom*mezcla.M)

    def _lib(self, cmp):
        """Library for compound specified calculation from generalized equation
        of Starling"""

        # Using molar base
        rhoc = 1/cmp.Vc/cmp.M
        w = cmp.f_acent

        if cmp.id in dat:
            # Using the compound specifie parameters given in Starling, [1]_
            Bo, Ao, Co, gamma, b, a, alpha, c, Do, d, Eo = dat[cmp.id]

            # The original data is in imperial units, converting to base unit
            Bo *= k.foot**3/k.pound
            Ao *= k.psi/1e3*(k.foot**3/k.pound)**2
            Co *= k.psi/1e3*k.Rankine**2*(k.foot**3/k.pound)**2
            gamma *= (k.foot**3/k.pound)**2
            b *= (k.foot**3/k.pound)**2
            a *= k.psi/1e3*(k.foot**3/k.pound)**3
            alpha *= (k.foot**3/k.pound)**3
            c *= k.psi/1e3*k.Rankine**2*(k.foot**3/k.pound)**3
            Do *= k.psi/1e3*k.Rankine**3*(k.foot**3/k.pound)**2
            d *= k.psi/1e3*k.Rankine*(k.foot**3/k.pound)**3
            Eo *= k.psi/1e3*k.Rankine**4*(k.foot**3/k.pound)**2

            # Conversion checked with
            # McFee, D.G., Mueller, K.H., Lielmezs, J.
            # Comparison of Benedict-Webb-Rubin, Starling and Lee-Kesler
            # Equations of State for Use in P-V-T Calculations
            # Thermochimica Acta 54(1-2) (1982) 9-25
            # 10.1016/0040-6031(82)85060-0

        else:
            # Using the generalized correlation given in Han and Starling, [2]_

            A = [None, 0.443690, 1.28438, 0.356306, 0.544979, 0.528629,
                 0.484011, 0.0705233, 0.504087, 0.0307452, 0.0732828, 0.006450]
            B = [None, 0.115449, -0.920731, 1.70871, -0.270896, 0.349261,
                 0.754130, -0.044448, 1.32245, 0.179433, 0.463492, -0.022143]

            Bo = (A[1]+B[1]*w)/rhoc
            Ao = (A[2]+B[2]*w)*k.R*cmp.Tc/rhoc
            Co = (A[3]+B[3]*w)*k.R*cmp.Tc**3/rhoc
            gamma = (A[4]+B[4]*w)/rhoc**2
            b = (A[5]+B[5]*w)/rhoc**2
            a = (A[6]+B[6]*w)*k.R*cmp.Tc/rhoc**2
            alpha = (A[7]+B[7]*w)/rhoc**3
            c = (A[8]+B[8]*w)*k.R*cmp.Tc**3/rhoc**2
            Do = (A[9]+B[9]*w)*k.R*cmp.Tc**4/rhoc
            d = (A[10]+B[10]*w)*k.R*cmp.Tc**2/rhoc**2
            Eo = (A[11]+B[11]*w*exp(-3.8*cmp.f_acent)) * \
                k.R*cmp.Tc**5/rhoc

        if self.kwargs["extended"]:
            e = (4.65593e-3 - 3.07393e-2*w +
                 5.58125e-2*w**2 - 3.40721e-3*exp(-7.72753*w-45.3152*w**2)) * \
                k.R*cmp.Tc**5/rhoc**2
            f = (0.697e-13 + 8.08e-13*w - 16e-13*w**2 -
                 0.363078e-13*exp(30.9009*w-283.68*w**2)) * \
                k.R*cmp.Tc**24/rhoc**2
            g = (2.2e-5-1.065e-4*w + 1.09e-5*exp(-26.024*w)) * \
                k.R*cmp.Tc**9/rhoc**2
            h = (-2.4e-11 + 11.8e-11*w - 2.05e-11*exp(-21.52*w)) * \
                k.R*cmp.Tc**18/rhoc**2
        else:
            e = 0
            f = 0
            g = 0
            h = 0

        kw = {}
        kw["Bo"] = Bo
        kw["Ao"] = Ao
        kw["Co"] = Co
        kw["gamma"] = gamma
        kw["b"] = b
        kw["a"] = a
        kw["alpha"] = alpha
        kw["c"] = c
        kw["Do"] = Do
        kw["d"] = d
        kw["Eo"] = Eo
        kw["e"] = e
        kw["f"] = f
        kw["g"] = g
        kw["h"] = h

        return kw

    def _mix(self, zi):
        """Mixing rules"""
        Ao = Co = Do = Eo = Bo = a = b = c = d = alfa = gamma = 0
        e = f = g = h = 0

        # Parameters without interaction parameter
        for i, x in enumerate(zi):
            Bo += x*self.Boi[i]
            a += x*self.ai[i]**(1/3)
            b += x*self.bi[i]**(1/3)
            c += x*self.ci[i]**(1/3)
            d += x*self.di[i]**(1/3)
            alfa += x*self.alfai[i]**(1/3)
            gamma += x*self.gammai[i]**0.5
            e += x*self.ei[i]**(1/3)
            f += x*self.fi[i]**(1/3)
            g += x*self.gi[i]
            h += x*self.hi[i]
        a = a**3
        b = b**3
        c = c**3
        d = d**3
        alfa = alfa**3
        gamma = gamma**2
        e = e**3
        f = f**3

        # Parameters with interaction parameter
        for i, (xi, kiji) in enumerate(zip(zi, self.kij)):
            for j, (xj, kij) in enumerate(zip(zi, kiji)):
                Ao += xi*xj * self.Aoi[i]**0.5 * self.Aoi[j]**0.5*(1-kij)
                Co += xi*xj * self.Coi[i]**0.5 * self.Coi[j]**0.5*(1-kij)**3
                Do += xi*xj * self.Doi[i]**0.5 * self.Doi[j]**0.5*(1-kij)**4
                Eo += xi*xj * self.Eoi[i]**0.5 * self.Eoi[j]**0.5*(1-kij)**5
        kw = {}
        kw["Ao"] = Ao
        kw["Bo"] = Bo
        kw["Co"] = Co
        kw["Do"] = Do
        kw["Eo"] = Eo
        kw["a"] = a
        kw["b"] = b
        kw["c"] = c
        kw["d"] = d
        kw["alfa"] = alfa
        kw["gamma"] = gamma
        kw["e"] = e
        kw["f"] = f
        kw["g"] = g
        kw["h"] = h
        return kw

    def _Z(self, zi=None, T=None, P=None, rho0=None, **kw):

        if not zi:
            zi = self.zi
        if not T:
            T = self.T
        if not P:
            P = self.P
        if not rho0:
            rho0 = 0

        if "Ao" not in kw:
            kw = self._mix(zi)

        Ao = kw["Ao"]
        Bo = kw["Bo"]
        Co = kw["Co"]
        Do = kw["Do"]
        Eo = kw["Eo"]
        a = kw["a"]
        b = kw["b"]
        c = kw["c"]
        d = kw["d"]
        alfa = kw["alfa"]
        gamma = kw["gamma"]
        e = kw["e"]
        f = kw["f"]
        g = kw["g"]
        h = kw["h"]

        def func(rho):
            return self.P.kPa - k.R*T*rho - \
                (Bo*k.R*T-Ao-Co/T**2+Do/T**3-Eo/T**4)*rho**2 - \
                (b*k.R*T - a - d/T - e/T**4 - f/T**23)*rho**3 - \
                alfa*(a + d/T + e/T**4 + f/T**23)*rho**6 - \
                (c/T**2 + g/T**8 + h/T**17)*rho**3 * \
                (1+gamma*rho**2)*exp(-gamma*rho**2)

        rho = fsolve(func, rho0)
        Z = self.P.kPa/rho/k.R/T
        return Z

    def _fug(self, xi, yi, T, P):
        """Fugacities of component in mixture calculation

        Parameters
        ----------
        xi : list
            Molar fraction of component in liquid phase, [-]
        yi : list
            Molar fraction of component in vapor phase, [-]

        Returns
        -------
        tital : list
            List with liquid phase component fugacities
        titav : list
            List with vapour phase component fugacities
        """
        kw = self._mix(xi)

        Tci = mezcla._arraylize("Tc")
        Vci = mezcla._arraylize("Vc")
        wi = mezcla._arraylize("f_acent")
        Mi = mezcla._arraylize("M")
        rho0 = RhoL_CostaldMix(T, mezcla.fraccion, Tci, wi, Vci, Mi)
        Zl = self._Z(xi, rho0=rho0, **kw)
        tital = self._fugacity(Zl, xi, **kw)
        Zv = self._Z(yi, rho0=0, **kw)
        titav = self._fugacity(Zv, yi, **kw)
        return tital, titav

    def _fugacity(self, Z, x, **kw):
        rho = self.P.kPa/Z/k.R/self.T

        Bo = kw["Bo"]
        a = kw["a"]
        b = kw["b"]
        c = kw["c"]
        d = kw["d"]
        alfa = kw["alfa"]
        gamma = kw["gamma"]
        e = kw["e"]
        f = kw["f"]
        g = kw["g"]
        h = kw["h"]

        tita = []
        for i, (xi, kiji) in enumerate(zip(x, self.kij)):
            suma = 0
            for j, (xj, kij) in enumerate(zip(x, kiji)):
                suma += xj*(
                    -(self.Aoi[j]*self.Aoi[i])**0.5*(1-kij) -
                    (self.Coi[j]*self.Coi[i])**0.5*(1-kij)**3/self.T**2 +
                    (self.Doi[j]*self.Doi[i])**0.5*(1-kij)**4/self.T**3 -
                    (self.Eoi[j]*self.Eoi[i])**0.5*(1-kij)**5/self.T**4)

            rhs = k.R*self.T*log(rho*k.R*self.T) + \
                rho*(Bo+self.Boi[i])*k.R*self.T + 2*rho*suma + \
                3*rho**2/2*((b**2*self.bi[i])**(1/3)*k.R*self.T -
                            (a**2*self.ai[i])**(1/3) -
                            (d**2*self.di[i])**(1/3)/self.T -
                            (e**2*self.ei[i])**(1/3)/self.T**4 -
                            (f**2*self.fi[i])**(1/3)/self.T**23) + \
                3*alfa*rho**5/5*((a**2*self.ai[i])**(1/3) +
                                 (d**2*self.di[i])**(1/3)/self.T +
                                 (e**2*self.ei[i])**(1/3)/self.T**4 +
                                 (f**2*self.fi[i])**(1/3)/self.T**23) + \
                3*rho**5/5*(a+d/self.T+e/self.T**4+f/self.T**23) * \
                    (alfa**2*self.alfai[i])**(1/3) + \
                rho**2*(3*(c**2*self.ci[i])**(1/3)/self.T**2 +
                        (self.gi[i]+2*g)/self.T**8 +
                        (self.hi[i]+2*h)/self.T**17)*(
                            (1/gamma/rho**2-(1/gamma/rho**2+0.5) *
                             exp(-gamma*rho**2))) -\
                (2*(c/self.T**2+g/self.T**8+h/self.T**17)*self.gammai[i]**0.5 /
                 gamma**1.5)*(
                    (1-exp(-gamma*rho**2)) *
                    (1+gamma*rho**2+gamma**2*rho**4/2))
            tita.append(exp(rhs/k.R/self.T)*xi)
        return tita

    def _Hexc(self, Z, T, **kw):
        rho = self.P.kPa/Z/k.R/self.T

        Ao = kw["Ao"]
        Bo = kw["Bo"]
        Co = kw["Co"]
        Do = kw["Do"]
        Eo = kw["Eo"]
        a = kw["a"]
        b = kw["b"]
        c = kw["c"]
        d = kw["d"]
        alfa = kw["alfa"]
        gamma = kw["gamma"]
        e = kw["e"]
        f = kw["f"]
        g = kw["g"]
        h = kw["h"]

        Hexc = (Bo*k.R*T - 2*Ao - 4*Co/T**2 + 5*Do/T**3 - 6*Eo/T**4)*rho + \
            (b*k.R*T - 1.5*a - 2*d/T - 7*e/2/T**4 - 13*f/T**23)*rho**2 + \
            alfa*(6/5*a + 7*d/5/T + 2*e/T**4 + 29*f/5/T**23)*rho**5 + \
            c/gamma/T**2*(
                3-(3+gamma/2*rho**2-gamma**2*rho**4)*exp(-gamma*rho**2)) + \
            g/gamma/T**8*(
                9-(9+7/2*gamma*rho**2-gamma**2*rho**4)*exp(-gamma*rho**2)) + \
            h/gamma/T**17*(
                18-(18*gamma*rho**2+gamma**2*rho**4)*exp(-gamma*rho**2))

        return Hexc


_all = [BWRS]


if __name__ == "__main__":
    from lib.mezcla import Mezcla

    # mix = Mezcla(2, ids=[10, 38, 22, 61], caudalUnitarioMolar=[0.3, 0.5, 0.05, 0.15])
    # eq = BWRS(340, 101325, mix)

    # mix = Mezcla(2, ids=[62], caudalUnitarioMolar=[1])
    # eq = BWRS(300, 101325, mix)
    # print(eq.V)
    # mix = Mezcla(2, ids=[50], caudalUnitarioMolar=[1])
    # eq = BWRS(200, 101325, mix)
    # print('%0.1f' % (eq.Vl.ccmol))
    # print('V = %0.1f' % (1/eq.Vl.m3kmol*eq.mezcla.M))
    # from lib.mEoS import H2S
    # st = H2S(T=200, P=101325)
    # print("MEoS ρ: ", st.rho)

    # mix = Mezcla(2, ids=[12, 40], caudalUnitarioMolar=[0.324, 0.676])
    # eq = BWRS(470, 1.4e6, mix)

    # mix = Mezcla(2, ids=[1001], caudalUnitarioMolar=[1])
    # eq1 = BWRS(250, 2e6, mix)
    # eq2 = BWRS(450, 1e7, mix)

    # mix = Mezcla(2, ids=[2, 3], caudalUnitarioMolar=[0.85, 0.15])
    # eq = BWRS(288.15, 6e6, mix)

    # mix = Mezcla(2, ids=[153, 41], caudalUnitarioMolar=[0.5, 0.5])
    # eq = BWRS(323.15, 25e3, mix)

    # mix = Mezcla(2, ids=[4, 2], caudalUnitarioMolar=[0.766, 0.234])
    # eq = BWRS(260, 1.7e6, mix)
    # print(eq.rhoL*eq.mezcla.M, eq.rhoL*40.93)
    # print(eq.mezcla.M)
    # print(eq.rhoG*eq.mezcla.M)
    # eq = BWRS(420, 1e7, mix)
    # print(eq.rhoL*eq.mezcla.M)
    # from lib.EoS.Cubic import PR, SRK
    # eq = PR(260, 1.7e6, mix)
    # print(eq.rhoL*eq.mezcla.M, eq.rhoL*40.93)
    # eq = SRK(260, 1.7e6, mix)
    # print(eq.rhoL*eq.mezcla.M, eq.rhoL*40.93)
    # eq = PR(420, 1e7, mix)
    # print(eq.rhoG*eq.mezcla.M)
    # eq = SRK(420, 1e7, mix)
    # print(eq.rhoG*eq.mezcla.M)

    # mix = Mezcla(2, ids=[2, 3, 4, 6, 5, 8, 46, 49, 50, 22],
                 # caudalUnitarioMolar=[1]*10)
    # eq = BWRS(293.15, 8.769e6, mix)

    # mezcla = Mezcla(1, ids=[4, 40], caudalUnitarioMasico=[26.92, 73.08])
    # P = unidades.Pressure(410.3, "psi")
    # T = unidades.Temperature(400, "F")
    # eq = BWRS(T, P, mezcla)

    # mezcla = Mezcla(1, ids=[10, 38, 22, 61], caudalUnitarioMasico=[0.3, 0.5, 0.05, 0.15])
    # eq = BWRS(340, 101325, mezcla)

    mezcla = Mezcla(3, ids=[4, 2], caudalMasico=1, fraccionMolar=[0.766, 0.234])
    eq = BWRS(200, 101325, mezcla, extended=False)
    print("q = ", eq.x)
    print("x = ", eq.xi)
    print("y = ", eq.yi)
    print("K = ", eq.Ki)
    print("Hexc: ", eq.HexcG, eq.HexcL)

    # mezcla = Mezcla(3, ids=[10], caudalMasico=1, fraccionMolar=[1])
    # eq = BWRS(340, 1.7e5, mezcla)
    # mezcla = Mezcla(3, ids=[10], caudalMasico=1, fraccionMolar=[1])
    # eq = BWRS(400, 1.7e5, mezcla)
