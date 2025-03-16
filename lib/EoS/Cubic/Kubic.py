#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


from scipy.constants import R

from lib.EoS.cubic import Cubic


class Kubic(Cubic):
    r"""Kubic cubic equation of state, [1]_

    This is a two parameter cubic equation of state with both parameters
    temperature dependent

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a(T)}{\left(V+c\right)^2}\\
        a(T) = \frac{27}{64} \frac{R^2T_c^2}{P_c}
        \left(\alpha^0\left(T_r\right)+\omega'\alpha^1\left(T_r\right)\right)\\
        \alpha^0 = -0.1514T_r + 0.7895 + \frac{0.3314}{T_r} +
        \frac{0.029}{T_r^2} + \frac{0.0015}{T_r^7}\\
        \alpha^1 = -0.237T_r - \frac{0.7846}{T_r} + \frac{1.0026}{T_r^2} +
        \frac{0.019}{T_r^7}\\
        b = \left(0.082-0.0715\omega'\right)\frac{RT_c}{P_c}\\
        c = \left(0.043\gamma^0\left(T_r\right)+0.0713\omega'\gamma^1
        \left(T_r\right)\right)\frac{RT_c}{P_c}\\
        \gamma^0 = 4.275051 - \frac{8.878889}{T_r} + \frac{8.508932}{T_r^2} -
        \frac{3.481408}{T_r^3} + \frac{0.576312}{T_r^4}\\
        \gamma^1 = 12.856404 - \frac{34.744125}{T_r} + \frac{37.433095}{T_r^2}
        - \frac{18.059421}{T_r^3} + \frac{3.514050}{T_r^4}\\
        \omega' = 0.000756+0.90984\omega+0.16226\omega^2+0.14549\omega^3\\
        \end{array}
    """

    __title__ = "Kubic (1982)"
    __status__ = "Kubic"
    __doi__ = {
        "autor": "Kubic, W.L.",
        "title": "A Modification of the Martin Equation of State for "
                 "Calculating Vapour-Liquid Equilibria",
        "ref": "Fluid Phase Equilibria 9 (1982) 79-97",
        "doi": "10.1016/0378-3812(82)85006-1"},

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""

        ai = []
        bi = []
        ci = []
        for cmp in self.componente:
            a, b, c = self._lib(cmp, T)

            ai.append(a)
            bi.append(b)
            ci.append(c)

        self.ai = ai
        self.bi = bi
        self.ci = ci

    def _GEOS(self, xi):
        am, bm, cm = self._mixture(None, xi, [self.ai, self.bi, self.ci])

        delta = 2*cm
        epsilon = cm**2
        return am, bm, delta, epsilon

    def _lib(self, cmp, T):

        Tr = T/cmp.Tc

        # Eq 26
        omega = 0.000756 + 0.90984*cmp.f_acent + 0.16226*cmp.f_acent**2 + \
            0.14549*cmp.f_acent**3

        # Eq 24
        g0 = 4.275051 - 8.878889/Tr + 8.508932/Tr**2 - 3.481408/Tr**3 + \
            0.576312/Tr**4

        # Eq 25
        g1 = 12.856404 - 34.744125/Tr + 37.433095/Tr**2 - 18.059421/Tr**3 + \
            3.514050/Tr**4

        # Eq 20
        a0 = -0.1514*Tr + 0.7895 + 0.3314/Tr + 0.029/Tr**2 + 0.0015/Tr**7

        # Eq 21
        a1 = -0.237*Tr - 0.7846/Tr + 1.0026/Tr**2 + 0.019/Tr**7

        b = (0.082 - 0.0713*omega)*R*cmp.Tc/cmp.Pc                      # Eq 22
        c = (0.043*g0 + 0.0713*omega*g1)*R*cmp.Tc/cmp.Pc                # Eq 23
        a = 27/64*(a0+omega*a1)*R**2*cmp.Tc**2/cmp.Pc                   # Eq 19

        return a, b, c


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = Kubic(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = Kubic(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
