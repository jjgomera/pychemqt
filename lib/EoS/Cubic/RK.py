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


from scipy.constants import R

from lib.bip import Kij, Mixing_Rule
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
    Example 4.3 from 1_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = RK(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2085 101.4'
    >>> eq = RK(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '97.3'
    """

    __title__ = "Redlich-Kwong (1949)"
    __status__ = "RK"
    __doi__ = {
        "autor": "Redlich, O., Kwong, J.N.S.",
        "title": "On the Thermodynamnics of Solutions V. An Equation of State."
                 " Fugacities of Gaseous Solutions",
        "ref": "Chem. Rev. 44 (1949) 233-244",
        "doi": "10.1021/cr60137a013"},

    def __init__(self, T, P, mezcla):
        ai = []
        bi = []
        for componente in mezcla.componente:
            a, b = self._lib(componente, T)
            ai.append(a)
            bi.append(b)

        self.kij = Kij(mezcla.ids)
        a, b = Mixing_Rule(mezcla.fraccion, [ai, bi], self.kij)

        self.ai = ai
        self.bi = bi
        self.b = b
        self.tita = a
        self.delta = b
        self.epsilon = 0

        super(RK, self).__init__(T, P, mezcla)

    def _lib(self, cmp, T):
        a = 0.42747*R**2*cmp.Tc**2/cmp.Pc
        alfa = (T/cmp.Tc)**-0.5
        b = 0.08664*R*cmp.Tc/cmp.Pc
        return a*alfa, b


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = RK(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = RK(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
