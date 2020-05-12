#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from math import exp

from scipy.constants import R

from lib.EoS.Cubic.RK import RK


class RKTwu(RK):
    r"""Equation of state of Redlich-Kwong with a modified alpha temperature
    dependence by Twu, (1995) [1]_.

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)}\\
        a = 0.427480263354\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.086640349965\frac{RT_c}{P_c}\\
        \alpha = alpha^{(0)} + \omega\left(\alpha^{(1)}-\alpha^{(0)}\right)\\
        \alpha^{(0)} = T_r^{-0.201158} \exp{0.141599\left(1-T_r^{2.29528}
        \right)}\\
        \alpha^{(1)} = T_r^{-0.660145} \exp{0.500315\left(1-T_r^{2.63165}
        \right)}\\
        \end{array}
    """
    __title__ = "Twu-Redlich-Kwong (1995)"
    __status__ = "RKTwu"
    __doi__ = {
        "autor": "Twu, C.H., Coon, J.E., Cunningham, J.R.",
        "title": "A New Generalized Alpha Function for a Cubic Equation of "
                 "State Part 2. Redlich-Kwong equation",
        "ref": "Fluid Phase Equilibria 105 (1995) 61-69",
        "doi": "10.1016/0378-3812(94)02602-w"},

    def _lib(self, cmp, T):
        """Modified parameteres correlations"""
        a = 0.42748023354*R**2*cmp.Tc**2/cmp.Pc
        alfa = self._alpha(cmp, T)
        b = 0.086640349965*R*cmp.Tc/cmp.Pc
        return a*alfa, b

    def _alpha(self, cmp, T):
        """Modified α expression"""
        Tr = T/cmp.Tc

        if Tr <= 1:
            alpha0 = Tr**(-0.201158)*exp(0.141599*(1-Tr**2.29528))      # Eq 17
            alpha1 = Tr**(-0.660145)*exp(0.500315*(1-Tr**2.63165))      # Eq 18
        else:
            alpha0 = Tr**(-1.10)*exp(0.441411*(1-Tr**(-1.30)))          # Eq 19
            alpha1 = Tr**(-2.31278)*exp(0.03258*(1-Tr**(-10.3128)))     # Eq 20

        # Eq 15
        alpha = alpha0 + cmp.f_acent*(alpha1-alpha0)

        return alpha


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = RKTwu(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = RKTwu(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
