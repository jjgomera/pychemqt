#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


from math import exp, log, log10

from scipy.constants import R

from lib.eos import EoS
from lib.EoS.Cubic import RK


class Grayson_Streed(EoS):
    r"""
    Hybrid model for vapor-liquid equilibria in hydrocarbon mixtures using
    three diferent parameters, given in [1]_.

    .. math::
        K_i = \frac{\nu_i^o\gamma_i}{\Phi_i^V}

    The fugacity coefficient of pure liquid is calculate4d with a Curl-Pitzer
    corresponding state correlation:

    .. math::
        \begin{array}[t]{l}
        \log{\nu_i^o} = \log{\nu_i^{(0)}} + \omega\log{\nu_i^{(1)}}\\
        \log{\nu_i^{(0)}} = A_0 + \frac{A_1}{T_r} + A_2T_r + A_3T_r^2 +
        A_4T_r^3 + \left(A_5+A_6T_r+A_7T_r^2\right)P_r +
        \left(A_8+A_9T_r\right)P_r^2 - \log{P_r}\\
        \log{\nu_i^{(1)}} = B_1 + B_1T_r-\frac{B_3}{T_r} + B_4T_r^3 -
        B_5\left(P_r-0.6\right)\\
        \end{array}

    The modified parameters of correlation are given in [2]_

    Liquid solution are considered regualar solutions, and the liquid activity
    coefficient are calculated using the Scatchard-Hildebrand correlation:

    .. math::
        \begin{array}[t]{l}
        \ln{\gamma_i} = \frac{V_i\left(\delta_i-\bar{\delta}\right)^2}{RT}\\
        \bar{\delta} = \frac{\sum_i{x_iV_i\delta_i}}{\sum_j{x_jV_j}}\\
        \end{array}

    Optionally can use the Flory-Huggins extension as is given in [3]_:

    .. math::
        \begin{array}[t]{l}
        \ln{\gamma_i} = \frac{V_i\left(\delta_i-\bar{\delta}\right)^2}{RT}
        + \ln{\Theta_i} + 1 - \Theta_i\\
        \Theta_i = \frac{V_i}{\sum_i{x_iV_i}}\\
        \end{array}

    The fugacity coefficient of vapor is calculated using the
    :doc:`Redlich-Kwong <lib.EoS.Cubic.RK>` cubic equation of state.

    This model is applicable to systems of non-polar hydrocarbons, for vapor-
    liquid equilibrium. It's use the value of gas density from Redlich-Kwong
    equation, but don't support the liquid density or the enthalpy calculation.

    Parameters
    ----------
    flory : boolean
        Use the Flory-Huggins extension to regular solutions

    Examples
    --------
    Example 7.2 from [4]_, pag 279, liquid-vapor flash equilibrium of a mixture
    of hydrogen, methane, benzene and toluene.

    >>> from lib.corriente import Mezcla
    >>> from lib import unidades
    >>> zi = [0.31767, 0.58942, 0.07147, 0.02144]
    >>> mix = Mezcla(2, ids=[1, 2, 40, 41], caudalUnitarioMolar=zi)
    >>> P = unidades.Pressure(485, "psi")
    >>> T = unidades.Temperature(100, "F")
    >>> eq = Grayson_Streed(T, P, mix, flory=0)
    >>> "β = %0.4f" % eq.x
    'β = 0.9106'
    >>> "xi = %0.4f %0.4f %0.4f %0.4f" % tuple(eq.xi)
    'xi = 0.0043 0.0576 0.7096 0.2286'
    >>> "yi = %0.4f %0.4f %0.4f %0.4f" % tuple(eq.yi)
    'yi = 0.3485 0.6417 0.0088 0.0011'
    """

    __title__ = "Grayson Streed (1961)"
    __status__ = "GS"

    __doi__ = (
      {
        "autor": "Chao, K.C., Seader, J.D.",
        "title": "A General Correlation of Vapor-Liquid Equilibria in "
                 "Hydrocarbon Mixtures",
        "ref": "AIChE J. 7(4) (1961) 598-605",
        "doi": "10.1002/aic.690070414"},
      {
        "autor": "Grayson, H.G., Streed, C.W.",
        "title": "Vapor-Liquid Equilibria for High Temperature, High Pressure "
                 "Hydrogen-Hydrocarbon Systems",
        "ref": "6th World Petroleum Congress, Frankfurt am Main, Germany, "
               "19-26 June (1963) 169-181",
        "doi": ""},
      {
        "autor": "Walas, S.M.",
        "title": "Phase Equiibria in Chemical Engineering",
        "ref": "Butterworth-Heinemann, 1985",
        "doi": "10.1016/C2013-0-04304-6"},
      {
        "autor": "Henley, E.J., Seader, J.D.",
        "title": "Equilibrium-Stage Separation Operations in Chemical "
                 "Engineering",
        "ref": "John Wiley & Sons, 1981",
        "doi": ""}
      )

    def __init__(self, T, P, mezcla, **kwargs):
        EoS.__init__(self, T, P, mezcla, **kwargs)
        self.rk = RK(T, P, mezcla)

        self.x, Zl, Zg, self.xi, self.yi, self.Ki = self._Flash()

        self.Zg = self.rk.Zg
        self.rhoG = self.rk.rhoG

        self.Zl = None
        self.rhoL = None
        # print("q = ", self.x)
        # print("x = ", self.xi)
        # print("y = ", self.yi)
        # print("K = ", self.Ki)

    def _nio(self, T, P):
        """Liquid fugacity coefficient"""

        nio = []
        for cmp in self.componente:
            Tr = T/cmp.Tc
            Pr = P/cmp.Pc

            # Modified Parameters from [2]_
            if cmp.id == 1:  # Hydrogen
                A = [1.50709, 2.74283, -.0211, .00011, 0, .008585, 0, 0, 0, 0]
            elif cmp.id == 2:  # Methane
                A = [1.36822, -1.54831, 0, 0.02889, -0.01076, 0.10486,
                     -0.02529, 0, 0, 0]
            else:
                A = [2.05135, -2.10899, 0, -0.19396, 0.02282, 0.08852, 0,
                     -0.00872, -0.00353, 0.00203]

            # Eq 3
            logn0 = A[0] + A[1]/Tr + A[2]*Tr + A[3]*Tr**2 + A[4]*Tr**3 + \
                (A[5]+A[6]*Tr+A[7]*Tr**2)*Pr + (A[8]+A[9]*Tr)*Pr**2 - log10(Pr)

            # Eq 4
            logn1 = -4.23893 + 8.65808*Tr - 1.2206/Tr - 3.15224*Tr**3 - \
                0.025*(Pr-0.6)

            # Eq 2
            nio.append(10**(logn0 + cmp.f_acent_mod*logn1))
            # The correlation use the modified acentric factor table

        return nio

    def _gi(self, xi, T):
        """Liquid activity coefficient"""
        Vm = 0
        for x, cmp in zip(xi, self.componente):
            Vm += x*cmp.wilson.m3mol

        phi = []
        for cmp in self.componente:
            phi.append(cmp.wilson.m3mol/Vm)

        sum1 = 0
        sum2 = 0
        for x, cmp in zip(xi, self.componente):
            sum1 += x*cmp.wilson.m3mol*cmp.SolubilityParameter
            sum2 += x*cmp.wilson.m3mol
        d_ = sum1/sum2

        # Eq 5
        gi = []
        for cmp, phii in zip(self.componente, phi):
            # Scatchard-Hildebrand regular solution activity-coefficient
            g = cmp.wilson.m3mol*(cmp.SolubilityParameter-d_)**2/R/T

            # Flory-Huggins extension
            if self.kwargs.get("flory", 0):
                g += log(phii) + 1 - phii

            gi.append(exp(g))

        return gi

    def _fug(self, xi, yi, T, P):
        gi = self._gi(xi, T)
        nio = self._nio(T, P)
        tital = [n*g for n, g in zip(nio, gi)]

        self.rk._cubicDefinition(T)
        Bi = [bi*P/R/T for bi in self.rk.bi]
        Ai = [ai*P/(R*T)**2 for ai in self.rk.ai]
        Z = self.rk._Z(yi, T, P)[-1]
        a, b, delta, epsilon = self.rk._GEOS(yi)
        B = b*P/R/T
        A = a*P/(R*T)**2
        titav = self.rk._fugacity(Z, yi, A, B, Ai, Bi)

        return tital, titav

    def _Z(self, xi, T, P):
        return None,


_all = [Grayson_Streed]


if __name__ == "__main__":
    from lib.corriente import Mezcla
    from lib import unidades

    # Example 7.2, pag 153
    # Method pag 107
    # mezcla = Mezcla(2, ids=[1, 2, 40, 41],
    #                 caudalUnitarioMolar=[0.31767, 0.58942, 0.07147, 0.02144])
    # P = unidades.Pressure(485, "psi")
    # T = unidades.Temperature(100, "F")
    # eq = Grayson_Streed(T, P, mezcla, flory=0)

    # mix = Mezcla(2, ids=[2, 3, 4, 6, 5, 8, 46, 49, 50, 22],
    #              caudalUnitarioMolar=[1]*10)
    # eq = Grayson_Streed(293.15, 5e6, mix, flory=0)

    # Example 4.2, pag 89
    mezcla = Mezcla(1, ids=[4, 40], caudalUnitarioMasico=[26.92, 73.08])
    P = unidades.Pressure(410.3, "psi")
    T = unidades.Temperature(400, "F")
    eq = Grayson_Streed(T, P, mezcla)
    print(eq.__dict__)

    # Example 6.7, Wallas, pag 342, dew point calculation
    # mezcla = Mezcla(2, ids=[23, 5], caudalUnitarioMolar=[0.607, 0.393])
    # P = unidades.Pressure(20, "atm")
    # eq = Grayson_Streed(400, P, mezcla, flory=0)
    # print(eq._Dew_T(P))
    # eq = Grayson_Streed(300, P, mezcla, flory=0)
    # print(eq._Dew_T(P))
    # print(eq._Bubble_T(P))

    # mix = Mezcla(2, ids=[2, 3, 4, 6, 5, 8, 46, 49, 50, 22],
    #              caudalUnitarioMolar=[1]*10)
    # eq = Grayson_Streed(293.15, 8.769e6, mix, flory=0)
