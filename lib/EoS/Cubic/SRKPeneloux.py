#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from scipy.constants import R

from lib import unidades
from lib.EoS.Cubic import SRK


class SRKPeneloux(SRK):
    r"""Equation of state of Soave-Redlich-Kwong (1972) with volume translation
    as define by Peneloux et al. [1]_

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{(V+b)\left(V+b+2c\right)}\\
        a = 0.42747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.08664\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-Tr^{0.5}\right)\\
        m = 0.48 + 1.574\omega-0.176\omega^2\\
        c = 0.40768 \frac{RT_c}{P_c} \left(0.29441-z_{RA}\right)\\
        \end{array}

    Examples
    --------
    Example 4.3 from [2]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = SRKPeneloux(300, 9.9742e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '93.1'
    >>> eq = SRKPeneloux(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vg.ccmol)
    '89.7'

    # Tiny desviation
    '2061 94.2'
    '90.6'
    """
    __title__ = "SRK-Peneloux (1982)"
    __status__ = "SRKPeneloux"
    __doi__ = (
      {
        "autor": "Péneloux, A., Rauzy, E., Fréze, R.",
        "title": "A Consistent Correction for Redlich-Kwong-Soave Volumes",
        "ref": "Fluid Phase Equilibria 8 (1982) 7-23",
        "doi": "10.1016/0378-3812(82)80002-2"},
      {
         "autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""})

    def _volumeCorrection(self):
        """Apply volume correction to the rhoL property"""
        c = 0
        for cmp, x in zip(self.componente, self.xi):
            c += 0.40768*R*cmp.Tc/cmp.Pc*(0.29441-cmp.rackett)           # Eq 8

        if self.rhoL:
            v = self.Vl.m3mol-c
            self.rhoL = unidades.MolarDensity(1/v, "molm3")
            self.Vl = unidades.MolarVolume(v, "m3mol")
        if self.rhoG:
            v = self.Vg.m3mol-c
            self.rhoG = unidades.MolarDensity(1/v, "molm3")
            self.Vg = unidades.MolarVolume(v, "m3mol")


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = SRKPeneloux(300, 9.9742e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
    eq = SRKPeneloux(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vg.ccmol))
