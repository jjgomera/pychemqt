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


from math import exp

from scipy.constants import R

from lib import unidades
from lib.EoS.Cubic import PR78


dat = {
    2: (-0.59970, -12.09400),
    3: (-0.29550, -14.93840),
    4: (-0.21240, -14.82030),
    6: (-0.13320, -15.82780),
    8: (-0.07520, -17.50000),
    10: (-0.00510, -19.36890),
    11: (0.00740, -19.63760),
    12: (0.05260, -21.08080),
    13: (0.07550, -22.27520),
    14: (0.08380, -23.07570),
    15: (0.10990, -26.80350),
    16: (0.10990, -26.80350),
    24: (-0.20500, -16.00870),
    5: (-0.18380, -15.48850),
    25: (-0.11930, -17.21510),
    26: (-0.15940, -16.61190),
    7: (-0.08400, -17.81840),
    36: (-0.13870, -16.91350),
    38: (-0.13870, -16.91350),
    100: (-0.11930, -17.21510),
    645: (0.18086, -25.50000),
    220: (-0.12570, -16.95870),
    216: (-0.18240, -16.29710),
    217: (-0.28070, -15.09170),
    215: (-0.23950, -15.59130),
    636: (-0.26540, -15.73780),
    218: (-0.24880, -15.90000),
    642: (-0.10340, -17.44440),
    643: (0.06780, -21.37850),
    232: (-0.18050, -16.27380),
    231: (-0.18040, -16.38390),
    229: (-0.13730, -19.00000),
    236: (-0.29510, -15.38060),
    # 1631: (-0.06320, -18.61910),
    # 1630: (-0.10810, -18.05000),
    # 1629: (-0.04860, -19.24500),
    # 1231: (-0.07480, -17.98790),
    # 1235: (-0.12380, -17.46360),
    1001: (0.03360, -21.33390),
    # 1633: (-0.11220, -17.54060),
    241: (-0.03780, -18.74630),
    243: (0.08890, -22.50000),
    245: (0.11150, -23.14840),
    671: (-0.32250, -15.18110),
    # 1,1,1,2,3,3,3-Hepta-fluoroethane: (-0.01990, -19.41640),
    # 1,1,1,3,3-Penta-fluoroethane: (-0.04360, -19.11100),
    # 1,1,2,2,3-Penta-fluoroethane: (-0.08560, -17.80890),
    # 1,1,1,2,3,3-Penta-fluoroethane: (-0.19980, -16.07300),
    # 1873: (-0.04130, -18.82970),
    225: (0.13670, -26.75710),
    692: (-0.24370, -15.55100),
    112: (-0.95420, -10.88110),
    119: (-0.01330, -19.32650),
    440: (-0.01330, -19.32650),
    22: (-0.23010, -15.70700),
    57: (-0.10100, -17.51680),
    23: (-0.18150, -16.31030),
    58: (-0.03890, -18.72320),
    28: (-0.08400, -17.81840),
    27: (-0.18150, -16.31030),
    29: (-0.08400, -17.81840),
    30: (-0.10150, -15.31780),
    31: (-0.21740, -15.10380),
    35: (-0.01330, -19.32650),
    56: (0.01950, -20.23130),
    83: (0.06870, -22.04110),
    367: (0.10670, -24.15240),
    452: (-0.18150, -16.31030),
    41: (-0.00160, -19.62810),
    45: (0.00920, -19.92970),
    43: (0.04650, -21.13620),
    42: (0.00920, -19.92970),
    44: (0.03800, -20.83460),
    71: (0.01950, -20.23130),
    470: (0.12180, -25.35890),
    98: (-0.79760, -11.57240),
    994: (-0.56390, -12.99720),
    46: (-0.74680, -11.96970),
    47: (-0.63680, -12.43320),
    208: (-0.65000, -12.41490),
    105: (-0.41900, -13.19250),
    99: (-0.59360, -12.69080),
    49: (-0.16000, -16.43410),
    102: (-0.36530, -17.24100),
    48: (-1.10870, -10.27790),
    51: (-0.06800, -18.12000),
    111: (0.06870, -22.04110),
    63: (0.14910, -26.07030),
    62: (0.19940, -30.05360),
    101: (-0.50110, -13.29410),
    178: (0.02900, -20.53300)}


class PRLinDuan(PR78):
    r"""Volume translation modification for Peng-Robinson equation of state by
    Lin and Duan, [1]_, in this equation the volumen translation is
    temperature-dependent

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V+c-b}-\frac{a}{(V+c)^2 + 2b\left(V+c\right) - b}\\
        a = 0.45724\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-Tr^{0.5}\right)\\
        m = 0.37464 + 1.54226\omega - 0.26992\omega^2\\
        c = c_cf\left(T_r\right)\\
        c_c = \left(0.3074-Z_c\right)\frac{RT_c}{P_c}\\
        f\left(T_r\right) = \beta + \left(1-\beta\right)
        \exp\left(\gamma\left|1-T_r\right|\right)\\
        \end{array}

    β and γ are compounds specific parameters fitted to experimental density
    data, if no available the paper give generalized correlation for this

    .. math::
        \begin{array}[t]{l}
        \beta = -2.8431\exp\left(-64.2184(0.3074-Z_c\right) + 0.1735\\
        \gamma = -99.2558 + 301.6201Z_c\\
        \end{array}
    """
    __title__ = "SRK-Peneloux (1982)"
    __status__ = "SRKPeneloux"
    __doi__ = {
        "autor": "Lin, H., Duan, Y.-Y.",
        "title": "Empirical Correction to the Peng-Robinson Equation of "
                 "State for the Saturated Region",
        "ref": "Fluid Phase Equilibria 233 (2005) 194-203",
        "doi": "10.1016/j.fluid.2005.05.008"},

    def _volumeCorrection(self):
        """Apply volume correction to the rhoL property"""

        c = 0
        for cmp, x in zip(self.componente, self.xi):
            if cmp.id in dat:
                b, g = dat[cmp.id]
            else:
                # Generalized correlation
                b = -2.8431*exp(-64.2184*(0.3074-self.Zc))+0.1735       # Eq 12
                g = -99.2558+301.6201*self.Zc                           # Eq 13

            Tr = self.T/cmp.Tc
            f = b+(1-b)*exp(g*abs(1-Tr))                                # Eq 10

            cc = (0.3074-cmp.Zc)*R*cmp.Tc/cmp.Pc                        # Eq 9
            c += cc * f                                                 # Eq 8

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
    eq = PRLinDuan(300, 9.9742e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
    eq = PRLinDuan(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vg.ccmol))
