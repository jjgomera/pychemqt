#!/usr/bin/python3
# -*- coding: utf-8 -*-

r'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


from lib.EoS.Cubic.vdW import vdW
from lib.EoS.Cubic.RK import RK
from lib.EoS.Cubic.RKTwu import RKTwu
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
from lib.EoS.Cubic.PRZV import PRZV
from lib.EoS.Cubic.PRTwu import PRTwu
from lib.EoS.Cubic.PRYuLu import PRYuLu
from lib.EoS.Cubic.PRLinDuan import PRLinDuan

from lib.EoS.Cubic.ALS1983 import ALS1983

from lib.EoS.Cubic.PT import PT
from lib.EoS.Cubic.PTV import PTV
from lib.EoS.Cubic.PTVC import PTVC

from lib.EoS.Cubic.TB import TB
from lib.EoS.Cubic.TBS import TBS

from lib.EoS.Cubic.Nasrifar import Nasrifar
from lib.EoS.Cubic.Kubic import Kubic


_all = [vdW,
        RK, RKTwu,
        SRK, SRKPeneloux, MSRK, SRKAPI,
        PR, PR78, PRMathiasCopeman,
        PRGasem, PRMelhem, PRAlmeida, PRSV, PRSV2, PRZV, PRTwu, PRYuLu,
        PRLinDuan,
        ALS1983,
        PT, PTV, PTVC,
        TB, TBS,
        Nasrifar, Kubic
        ]
# _all= [SRK_API]


# Add references from equations hardcoded in __doi__ property
__doi__ = {}
for obj in _all:
    if obj.__doi__:
        __doi__[obj.__name__] = obj.__doi__
