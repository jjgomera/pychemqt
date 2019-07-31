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

from lib.EoS.Cubic.PR import PR


class PRTwu(PR):
    r"""Peng-Robinson cubic equation of state with a modified dependence of 
    temperature by Twu et al.

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.457235528921\frac{R^2T_c^{2.5}}{P_c}\alpha\\
        b = 0.0777960739039\frac{RT_c}{P_c}\\
        \alpha = alpha^{(0)} + \omega\left(\alpha^{(1)-\alpha^{(0)}\right)\\
        \alpha^{(0)} = T_r^{-0.171813} \exp{0.125283\left(1-T_r^{1.77634}
        \right)\\
        \alpha^{(1)} = T_r^{-0.607352} \exp{0.511614\left(1-T_r^{2.20517}
        \right)\\
        \end{array}

    Examples
    --------
    Example 4.3 from 1_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = PRTwu(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2039 86.8'
    >>> eq = PRTwu(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '84.1'

    It give better result than in example references
    """

    __title__="PR Twu (1995)"
    __status__="PRTwu"
    __doi__ = {
        "autor": "Twu, C.H., Coon, J.E., Cunningham, J.R.",
        "title": "A New Generalized Alpha Function for a Cubic Equation of "
                 "State Part 1. Peng-Robinson equation",
        "ref": "Fluid Phase Equilibria 105 (1995) 49-59",
        "doi": "10.1016/0378-3812(94)02601-v"},

    OmegaA = 0.457235528921 
    OmegaB = 0.0777960739039

    def _alfa(self, cmp, T):
        Tr = T/cmp.Tc

        alpha0 = Tr**(-0.171813)*exp(0.125283*(1-Tr**1.77634))          # Eq 11
        alpha1 = Tr**(-0.607352)*exp(0.511614*(1-Tr**2.20517))          # Eq 12

        # Eq 9
        alpha = alpha0 + cmp.f_acent*(alpha1-alpha0)

        return 0, alpha


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRTwu(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRTwu(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
