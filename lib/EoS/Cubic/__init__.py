#!/usr/bin/python3
# -*- coding: utf-8 -*-

r'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Other equations don't implemented


Clausius (1881)
---------------

Historical equation of state

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{T \left(V+c\right)^2}\\
    \Omega_a = 27/64\\
    \Omega_b = Z_c-0.25\\
    \Omega_c = 3/8-Z_c\\
    \end{array}

Clausius, R., “Ueber das Verhaiten der Kohlensaure in Bezug auf Druck,
Volumen und Temperatur”, Ann. Phys. Chem. 9 (1880) 337-359

Wilson (1964)
-------------

Modified RK temperature dependence

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)}\\
    \Omega_a = \Omega_{ac}\left(1+m\left(\frac{1}{T_r}-1\right)\right)T_r\\
    m = 1.57 + 1.62\omega\\
    \end{array}

Wilson, G.M. Vapor-Liquid Equilibria, Correlations by Means of a Modified
Redlich-Kwong Equation of State. Adv. Cryog. Eng. 9(D2) (1964) 168-176.


HPW (1976)
----------

Hederer, Peter, Wenzel equation of state, published like PR in 1976.

.. math::
    P = \frac{RT}{V-b}-\frac{aT^{\alpha}}{V\left(V+b\right)}\\

α, a and b are compound specific parameters, the paper give correlations for
parameter for several homologous series like alkanes, alkenes or alkynes, but
there isn't any general correlation.

Hederer, H., Peter, S., Wenzel, H. Calculation of Thermodynamic Properties from
a Modified Redlich-Kwong Equation of State. Chem. Eng. J., 11 (1976) 183-190,
http://dx.doi.org/10.1016/0300-9467(76)80039-1


vdW Adachi (1984)
-----------------

Adachi modification to van der Waals original correlation using a logaritmic
temperature dependence for a

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{V^2}\\
    a(T) = \frac{27}{64}10^{m\left(1-T_r\right)}\\
    m = 0.228165 + 0.791981\omega - 0.648552*\omega^2 + 0.654505*\omega^3\\
    b = 0.125\frac{RT_c}{P_c}\\
    \end{array}


Valid only by VLE calculation, not applicable in a general use

Adachi, Y., Lu, B.C.-Y. Simplest Equation of State for Vapor-Liquid Equilibrium
calculation: a Modification of the van der Walls Equation. AIChE J. 30(6)
(1984) 991-993, http://dx.doi.org/10.1002/aic.690300619


Usdin-McAuliffe (1976)
----------------------
Incomplete study to improve SRK72 EoS liquid density prediction

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{V\left(V+d\right)}\\
    a = \frac{R^2T_c^2}{A_cP_c}\\
    b = \frac{RT_c}{B_cP_c}\\
    d = \frac{RT_c}{D_cP_c}\\
    \alpha^{0.5} = 1 + m\left(1-\sqrt{T_r}\right)\\
    \end{array}

If T_r ≤ 0.7:

.. math::
    m = 0.48049 + 4.516\omega Z_c^* + \left(0.67713\left(\omega-0.35\right)
    -0.02\right)x\left(T_r-0.7\right)

If 0.7 < T_r ≤ 1.0:

.. math::
    m = 0.48049 + 4.516\omega Z_c^* + \left(37.7846\omega{Z_c^*}^3+0.78662
    \right)\left(T_r-0.7\right)^2

:math:`D_c` is the most positive root of

.. math::
    D_c^3 + D_c^2\left(6Z_c-1\right) + D_c\left(4Z_c-1\right)3Z_c +
    \left(8Z_c-3\right)Z_c^2 = 0

:math:`B_c` and :math:`A_c` are calculated known :math:`D_c` from

.. math::
    \begin{array}[t]{l}
    B_c = D_c + 3Z_c - 1\\
    A_c = Z_c^3/B_c
    \end{array}

Usdin, E., McAuliffe, J.C. A One Parameter Family of Equations of State. Chem.
Eng. Sci., 31(ll) (1976) 1077-1084,
http://dx.doi.org/10.1016/0009-2509(76)87030-3


Each equation is specially suitable for different compounds, for example, the
Schmidt-Wenzel (SW) equation (1980) and the Adachi-Lu-Sugie (ALS) equation
(1983) are good for methane to n-decane. The Yu-Lu (YL) equation (1987) was
designed for asymmetric nonpolar mixtures, but not for polar substances. The
Iwai-Margerum-Lu (IML) equation ( 1987) was developed for polar substances, but
not suitable for nonpolar substances with large molecular weight.


'''


from lib.EoS.Cubic.vdW import vdW
from lib.EoS.Cubic.RK import RK
from lib.EoS.Cubic.SRK import SRK
from lib.EoS.Cubic.SRKPeneloux import SRKPeneloux
from lib.EoS.Cubic.MSRK import MSRK
from lib.EoS.Cubic.SRKAPI import SRKAPI
from lib.EoS.Cubic.PR import PR
from lib.EoS.Cubic.PR78 import PR78
from lib.EoS.Cubic.PRMathiasCopeman import PRMathiasCopeman
from lib.EoS.Cubic.PRGasem import PRGasem
from lib.EoS.Cubic.PRMelhem import PRMelhem
from lib.EoS.Cubic.PRAlmeida import PRAlmeida
from lib.EoS.Cubic.PRSV import PRSV
from lib.EoS.Cubic.PRSV2 import PRSV2
from lib.EoS.Cubic.PRTwu import PRTwu
from lib.EoS.Cubic.PRYuLu import PRYuLu

from lib.EoS.Cubic.ALS1983 import ALS1983

from lib.EoS.Cubic.PT import PT
from lib.EoS.Cubic.PTV import PTV
from lib.EoS.Cubic.PTVC import PTVC

from lib.EoS.Cubic.TB import TB
from lib.EoS.Cubic.TBS import TBS


_all = [vdW,
        RK,
        SRK, SRKPeneloux, MSRK, SRKAPI,
        PR, PR78, PRMathiasCopeman,
        PRGasem, PRMelhem, PRAlmeida, PRSV, PRSV2, PRTwu, PRYuLu,
        ALS1983,
        PT, PTV, PTVC,
        TB, TBS
        ]
# _all= [SRK_API]


# Add references from equations hardcoded in __doi__ property
__doi__ = {}
for obj in _all:
    if obj.__doi__:
        __doi__[obj.__name__] = obj.__doi__
