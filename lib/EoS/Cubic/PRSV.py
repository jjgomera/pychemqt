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

from lib.bip import Kij
from lib.EoS.cubic import Cubic, CubicHelmholtz


dat = {
    46 : 0.01996,
    47 : 0.01512,
    49 : 0.04285,
    63 : 0.00100,
    62 : -0.06635,
    104 : 0.01989,

    2 : -0.00159,
    3 : 0.02669,
    22 : 0.04400,
    4 : 0.03136,
    6 : 0.03443,
    8 : 0.03946,
    9 : 0.04303,
    10 : 0.05104,
    11 : 0.04648,
    12 : 0.04464,
    13 : 0.04104,
    14 : 0.04510,
    15 : 0.02919,
    16 : 0.05426,
    17 : 0.04157,
    18 : 0.02686,
    19 : 0.01892,
    20 : 0.02665,
    21 : 0.04048,
    90 : 0.08291,
    38 : 0.07023,
    # : 0.01805,  # Bicyclohexyl

    40 : 0.07019,
    41 : 0.03849,
    45 : 0.03994,
    44 : 0.01277,
    182 : 0.01173,
    70 : 0.02715,
    75 : -0.01384,
    185 : 0.03297,
    191 : -0.01842,
    192 : -0.01639,
    194 : 0.11487,
    406 : 0.05955 ,
    # : -0.01393,  # 9,10-Dihydrophenanthrene

    140 : -0.00888,
    448 : 0.00554,
    304 : 0.01681,
    # : 0.03558,  # 2-pentanone
    # : 0.04113,  # 3-pentanone
    771 : 0.00984,
    # : 0.02321,    # 3-hexanone
    # : 0.04005,    # Dimethylbutanone
    # : 0.02731,    # 2-heptanone
    # : 0.02002,    # 5-nonanone

    117 : -0.16816,
    134 : -0.03374,
    146 : 0.21419,
    145 : 0.23264,
    160 : 0.33431,
    159 : 0.39045,
    314 : 0.37200,
    316 : 0.43099,
    313 : 0.36781,
    335 : -0.00237,
    360 : 0.82940,
    396 : 0.80898,

    133 : 0.05717,
    486 : 0.16948,
    # : 0.02300,     # Methyl n-propyl ether
    # : 0.04123,     # Methyl i-propyl ether
    # : 0.01622,     # Methyl n-butyl ether
    456 : 0.05129,
    318 : -0.01668,
    778 : -0.03162,
    337 : 0.03751,
    344 : 0.01610,

    226 : -0.10299,
    125 : -0.13991,
    130 : -0.19724,
    144 : 0.18999,
    688 : -0.42503,
    270 : 0.14326,
    271 : 0.06001,
    # : -0.09508,    # 2-Methoxypropionitrile
    # : 0.13440,      # 2-Methyl-2-Propylamine
    281 : 0.03961,
    295 : 0.06946,
    164 : -0.03471,
    # : 0.11367,     # N-Methylpyrrolidone
    319 : 0.02752,
    796 : -0.00901,
    346 : 0.24705,
    # : 0.06043,   # Thianaphthene
}

class PRSV(Cubic):
    r"""Peng-Robinson cubic equation of state with a modified dependence of 
    temperature by Stryjek-Vera
    
    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^{2.5}}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + k\left(1-Tr^{0.5}\right)\\
        k = k_0+k_1\left(1+\sqrt{T_r}\right)\left(0.7-T_r\right)\\
        k_0 = 0.378893+1.4897153\omega-0.17131848\omega^2+0.0196554\omega^3\\
        \end{array}

    :math:`k_1` is a parameter characteristic to eatch compound

    Examples
    --------
    Example 4.3 from 1_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = PRSV(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2040 86.8'
    >>> eq = PRSV(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '84.2'
    """

    __title__="PRSV (1986)"
    __status__="PRSV"
    __doi__ = {
        "autor": "Stryjek, R., Vera, J.H.",
        "title": "PRSV: An Improved Peng—Robinson Equation of State for Pure "
                 "Compounds and Mixtures",
        "ref": "Can. J. Chem. Eng. 1986, 64: 323–333",
        "doi":  "10.1002/cjce.5450640224"},

    def __init__(self, T, P, mezcla):

        self.T = T
        self.P = P
        self.mezcla = mezcla

        ai=[]
        bi=[]
        for cmp in mezcla.componente:
            Tr = T/cmp.Tc

            # Eq 11
            ko = 0.378893 + 1.4897153*cmp.f_acent - \
                0.17131848*cmp.f_acent**2 + 0.0196554*cmp.f_acent**3

            if cmp.id in dat and Tr < 1:
                k1 = dat[cmp.id]
            else:
                k1 = 0
            k = ko + k1*(1+Tr**0.5)*(0.7-Tr)

            alfa = (1 + k*(1-Tr**0.5))**2
            ao = 0.457235*R**2*cmp.Tc**2/cmp.Pc
            b = 0.077796*R*cmp.Tc/cmp.Pc

            ai.append(ao*alfa)
            bi.append(b)

        am, bm = self._mixture(None, mezcla.ids, [ai, bi])

        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = 2*bm
        self.epsilon = -bm**2

        # tdadt=0
        # for i in range(len(mezcla.componente)):
            # for j in range(len(mezcla.componente)):
                # tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])
        # self.dTitadT=tdadt
        # self.u=2
        # self.w=-1

        super(PRSV, self).__init__(T, P, mezcla)


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRSV(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRSV(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))

