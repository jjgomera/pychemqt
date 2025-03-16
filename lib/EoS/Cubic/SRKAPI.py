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


from lib.EoS.Cubic.SRK import SRK


# Table 8D1.2
S1 = {
    62: 1.243997,
    63: 0.975515,
    219: 0.592450,
    227: 0.529899,
    137: 0.763226,
    117: 1.828343,
    134: 1.678665,
    146: 0.169684,
    145: 0.140334,
    160: 0.293950,
    159: 0.703883,
    450: 0.601957,
    161: 0.745244,
    456: 0.956082,
    # 1363: 0.886894,
    626: 0.905443,
    337: 1.025312}

# Table 8D1.3
S2 = {
    # Nonhydrocarbons
    1: -0.025891,
    62: -0.201789,
    63: -0.087598,
    50: 0.010699,
    46: -0.011016,
    48: -0.025280,
    49: -0.004474,
    98: -0.011721,

    # Paraffins
    2: -0.012223,
    3: -0.012416,
    4: -0.003791,
    6: 0.003010,
    5: 0.006209,
    8: -0.000636,
    7: -0.005345,
    9: 0.000695,
    10: -0.007459,
    52: -0.003898,
    53: -0.009786,
    54: -0.006195,
    55: -0.004579,
    11: -0.003031,
    79: -0.002132,
    80: -0.002806,
    432: -0.003786,
    433: -0.001338,
    434: 0.000811,
    435: -0.004779,
    436: 0.000168,
    12: -0.000821,
    592: -0.000599,
    593: -0.002424,
    594: -0.001943,
    595: -0.001621,
    596: -0.001728,
    597: -0.001744,
    81: -0.002094,
    590: -0.004039,
    591: -0.002870,
    541: -0.001434,
    82: -0.004045,
    599: -0.003849,
    600: -0.004313,
    603: 0.035734,
    598: -0.002745,
    601: -0.002499,
    602: -0.004598,
    13: 0.005435,
    371: -0.009265,
    368: 0.012650,
    370: -0.000017,
    372: 0.002328,
    373: -0.002088,
    374: -0.001878,
    375: -0.005733,
    14: 0.003324,
    15: -0.012698,
    16: -0.001931,
    17: 0.010085,
    19: -0.033625,
    20: -0.002741,
    90: -0.003612,
    92: -0.014476,

    # Naphthenes
    258: -0.017047,
    705: -0.014555,
    36: -0.003383,
    37: -0.001461,
    59: -0.002891,
    93: -0.013551,
    94: -0.006564,
    95: -0.009015,
    96: -0.001337,
    97: -0.017759,
    179: 0.012036,
    583: -0.008055,
    582: -0.014972,
    578: -0.022620,
    579: -0.024026,
    580: -0.026470,
    412: 0.003368,
    38: -0.004637,
    39: -0.000172,
    87: 0.007389,
    88: 0.014432,
    577: 0.007290,
    89: 0.010622,
    85: 0.014903,
    86: 0.007723,
    184: 0.021283,
    366: -0.033961,
    190: 0.028388,
    390: 0.025957,
    391: -0.000764,
    414: -0.003096,
    356: -0.004421,
    584: 0.467088,

    # Olefins
    22: -0.002805,
    23: -0.006163,
    24: 0.001178,
    25: -0.005287,
    26: -0.016832,
    27: 0.000612,
    29: -0.003670,
    30: -0.002330,
    31: -0.001574,
    32: -0.003510,
    34: -0.011281,
    33: -0.004255,
    35: -0.003219,
    552: 0.003966,

    553: 0.004891,
    554: -0.005113,
    555: 0.000441,
    767: 0.005190,
    556: 0.011483,
    # 1352: -0.012116,
    768: -0.000858,
    557: 0.034720,
    558: -0.024994,
    559: 0.001675,
    560: 0.000704,
    562: 0.016761,
    561: 0.004265,
    563: -0.008468,
    766: 0.015833,
    56: -0.006221,
    # 1395: 0.018598,
    # : 0.005120,  cis-3-Methyl-3-hexene
    # : -0.001891,  2,4-Dimethyl-1-pentene
    # : 0.008488,  2,4-Dimethyl-2-pentene
    # : -0.005706,  cis-4,4-Dimethyl-2-pentene
    # 1398: 0.000205,
    # : 0.006036,  3-Methyl-2-ethyl-1-butene
    83: -0.002307,
    # : 0.004362,  2,4,4-Trimethyl-I-pentene
    367: -0.000564,
    392: -0.005588,
    399: 0.010786,
    402: 0.008653,
    408: 0.005573,
    411: 0.008227,
    413: 0.000957,
    415: 0.010177,
    420: 0.021612,

    # Diolefins and acetylenes
    28: -0.009137,
    296: 0.011814,
    724: 0.043579,
    297: 0.052253,
    298: 0.057869,
    # 1308: 0.023742,
    300: 0.026266,
    61: 0.023298,
    # 1338: 0.010481,
    324: -0.016851,
    # 1346: 0.000656,
    65: 0.001831,
    66: -0.002337,
    151: -0.109294,
    # 1349: 0.022376,

    # Aromatics
    40: -0.000318,
    41: -0.005125,
    45: -0.004227,
    43: -0.005645,
    42: -0.006569,
    44: -0.010556,
    70: -0.007551,
    71: -0.008698,
    75: -0.014458,
    76: -0.010675,
    77: -0.014912,
    73: 0.036093,
    72: 0.047885,
    74: 0.036093,
    # : -0.036510,  1-methyl-3-ethenyl benzene
    78: -0.005215,
    377: -0.013419,
    378: 0.018480,
    379: 0.018084,
    861: 0.010750,
    862: 0.011193,
    383: -0.006384,
    381: 0.004533,
    382: 0.000849,
    384: -0.005169,
    614: 0.019944,
    615: -0.014791,
    885: 0.022869,
    911: 0.011198,
    178: 0.019670,

    # Diaromatics and other hydrocarbon rings
    386: -0.006020,
    387: 0.011339,
    376: 0.018025,
    181: -0.028120,
    185: -0.005006,
    191: 0.056699,
    192: -0.000472,
    # 1495: -0.014329,
    # 1490: 0.045513,
    888: -0.020675,
    69: -0.006625,
    194: 0.000398,
    406: 0.000798,

    # Sulfur compounds
    227: 0.141244,
    137: 0.003516,
    219: 0.026420,

    # Oxygenated compounds
    117: -0.430885,
    134: -0.216396,
    146: 1.188769,
    145: 1.269059,
    160: 1.005612,
    159: 0.503560,
    450: 0.600508,
    161: 0.484002,
    456: -0.053074,
    # : 0.046280,  Ethyl-tertbutyl ether
    626: 0.010549,
    337: -0.033494}


# TODO: Add interaction coefficient from [1]_


class SRKAPI(SRK):
    r"""Cubic equation of state derived from Soave-Redlich-Kwong and described
    in [1]_, Procedure 8D1.1.

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)}\\
        a = 0.42747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.08664\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + S_1\left(1-\sqrt{Tr}\right) +
        S_2\frac{\left(1-\sqrt{Tr}\right)}{\sqrt{T_r}}\\
        S_1 = 0.48508 + 1.55171\omega - 0.15613\omega^2\\
        \end{array}

    :math:`S_1` generalized correlation is given in [2]_

    This method define too its custom mixture interaction parameters from [3]_
    and [4]_
    """

    __title__ = "SRK-API (1993)"
    __status__ = "SRK-API"
    __doi__ = (
        {"autor": "API",
         "title": "Technical Data book: Petroleum Refining 6th Edition",
         "ref": "",
         "doi": ""},
        {"autor": "Graboski, M.S., Daubert, T.E.",
         "title": "A Modified Soave Equation of State for Phase Equilibrium "
                  "Calculations. 1. Hydrocarbon Systems",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 17(4) (1978) 443-448",
         "doi": "10.1021/i260068a009"},
        {"autor": "Graboski, M.S., Daubert, T.E.",
         "title": "A Modified Soave Equation of State for Phase Equilibrium "
                  "Calculations. 2. Systems Containing CO₂, H₂S, N₂ and CO",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 17(4) (1978) 448-454",
         "doi": "10.1021/i260068a010"},
        {"autor": "Graboski, M.S., Daubert, T.E.",
         "title": "A Modified Soave Equation of State for Phase Equilibrium "
                  "Calculations. 3. Systems Containing Hydrogen",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 18(2) (1979) 300-306",
         "doi": "10.1021/i260070a022"})

    def _alfa(self, cmp, T):
        """Special alpha temperature dependence as give in [1]_"""

        Tr = T/cmp.Tc

        if cmp.id in S2:
            s2 = S2[cmp.id]
        else:
            s2 = 0

        if cmp.id in S1:
            s1 = S1[cmp.id]
        else:
            # Eq 8D1.1-7
            s1 = 0.48508 + 1.555171*cmp.f_acent - 0.15613*cmp.f_acent**2

        alfa = (1 + s1*(1-Tr**0.5) + s2*(1-Tr**0.5)/Tr**0.5)**2   # Eq 8D1.1-6
        return alfa


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = SRKAPI(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = SRKAPI(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
