#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from scipy.constants import R

from lib.EoS.cubic import Cubic


class vdW(Cubic):
    r"""Equation of state of van der Waals (1873), [1]_.

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V^2}\\
        a = 0.421875\frac{R^2T_c^2}{P_c}\\
        b = 0.125\frac{RT_c}{P_c}\\
        \end{array}

    This equation is not accuracy and it's implemented only by its historic
    relevance because it's the first cubic equation of state.

    Examples
    --------
    Example 4.3 from [2]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = vdW(300, 9.9742e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '145.4'
    >>> eq = vdW(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vg.ccmol)
    '135.5'
    """

    __title__ = "van der Waals (1890)"
    __status__ = "vdW"
    __doi__ = (
        {"autor": "van der Waals, J.D.",
         "title": "Over de Continuiteit van den Gas- En Vloestoftoestand",
         "ref": "Dissertation, Leiden University, Leiden, Niederlande, 1873",
         "doi": ""},
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""})

    def _cubicDefinition(self, T=None):
        """Definition of individual components coefficients"""

        # Schmidt-Wenzel factorization of terms
        self.u = 0
        self.w = 0

        ai = []
        bi = []
        for cmp in self.componente:
            a = 0.421875*R**2*cmp.Tc**2/cmp.Pc
            b = 0.125*R*cmp.Tc/cmp.Pc
            ai.append(a)
            bi.append(b)

        self.ai = ai
        self.bi = bi

    def _GEOS(self, xi):
        am, bm = self._mixture(None, xi, [self.ai, self.bi])

        delta = 0
        epsilon = 0

        return am, bm, delta, epsilon


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = vdW(300, 9.9742e5, mix)
    print(eq.Vl.ccmol)
    eq = vdW(300, 42.477e5, mix)
    print(eq._Z())
    print('%0.1f' % (eq.Vg.ccmol))
