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


from scipy.constants import R

from lib.EoS.cubic import Cubic


class RK(Cubic):
    r"""Equation of state of Redlich-Kwong (1949), [1]_.

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)\sqrt{T}}\\
        a = 0.42747\frac{R^2T_c^{2.5}}{P_c}\\
        b = 0.08664\frac{RT_c}{P_c}\\
        \end{array}

    Examples
    --------
    Example 4.3 from [2]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = RK(300, 9.9742e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '101.4'
    >>> eq = RK(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vg.ccmol)
    '97.3'
    """

    __title__ = "Redlich-Kwong (1949)"
    __status__ = "RK"
    __doi__ = (
      {
        "autor": "Redlich, O., Kwong, J.N.S.",
        "title": "On the Thermodynamnics of Solutions V. An Equation of State."
                 " Fugacities of Gaseous Solutions",
        "ref": "Chem. Rev. 44 (1949) 233-244",
        "doi": "10.1021/cr60137a013"},
      {
         "autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""})

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""

        # Schmidt-Wenzel factorization of terms
        self.u = 1
        self.w = 0

        ai = []
        bi = []
        for cmp in self.componente:
            a, b = self._lib(cmp, T)

            ai.append(a)
            bi.append(b)

        self.ai = ai
        self.bi = bi

    def _lib(self, cmp, T):
        a0 = 0.42747*R**2*cmp.Tc**2/cmp.Pc
        alfa = (T/cmp.Tc)**-0.5
        b = 0.08664*R*cmp.Tc/cmp.Pc
        return a0*alfa, b

    def _GEOS(self, xi):
        am, bm = self._mixture("RK", xi, [self.ai, self.bi])

        delta = bm
        epsilon = 0
        return am, bm, delta, epsilon


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    from lib import unidades

    # mix = Mezcla(2, ids=[10, 38, 22, 61],
    #              caudalUnitarioMolar=[0.3, 0.5, 0.05, 0.15])
    # eq = RK(340, 101325, mix)

    mezcla = Mezcla(2, ids=[1, 2, 40, 41],
                    caudalUnitarioMolar=[0.31767, 0.58942, 0.07147, 0.02144])
    P = unidades.Pressure(485, "psi")
    T = unidades.Temperature(100, "F")
    eq = RK(T, P, mezcla, flory=1)
