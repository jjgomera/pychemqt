#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


from lib.meos import MEoS, MEoSBlend

# Noble Gases
from lib.mEoS import He
from lib.mEoS import Ne
from lib.mEoS import Ar
from lib.mEoS import Kr
from lib.mEoS import Xe

# Gases
from lib.mEoS import H2
from lib.mEoS import D2
from lib.mEoS import pH2
from lib.mEoS import oH2
from lib.mEoS import N2
from lib.mEoS import O2
from lib.mEoS import F2
from lib.mEoS import H2O
from lib.mEoS import D2O
from lib.mEoS import CO2
from lib.mEoS import CO
from lib.mEoS import N2O
from lib.mEoS import SO2
from lib.mEoS import COS
from lib.mEoS import NH3
from lib.mEoS import H2S

# Alkanes

from lib.mEoS import CH4
from lib.mEoS import C2
from lib.mEoS import C3
from lib.mEoS import nC4
from lib.mEoS import iC4
from lib.mEoS import nC5
from lib.mEoS import neoC5
from lib.mEoS import iC5
from lib.mEoS import nC6
from lib.mEoS import iC6
from lib.mEoS import nC7
from lib.mEoS import nC8
from lib.mEoS import iC8
from lib.mEoS import nC9
from lib.mEoS import nC10
from lib.mEoS import nC11
from lib.mEoS import nC12

# Naphthenes
from lib.mEoS import Cyclopropane
from lib.mEoS import Cyclopentane
from lib.mEoS import Cyclohexane
from lib.mEoS import C1Cyclohexane
from lib.mEoS import C3Cyclohexane

# Alkenes
from lib.mEoS import Benzene
from lib.mEoS import Toluene
from lib.mEoS import oXylene
from lib.mEoS import pXylene
from lib.mEoS import pXylene
from lib.mEoS import EthylBenzene
from lib.mEoS import Ethylene
from lib.mEoS import Propylene
from lib.mEoS import Butene_1
from lib.mEoS import iButene
from lib.mEoS import Cis_2_butene
from lib.mEoS import Trans_2_butene
from lib.mEoS import Propyne
from lib.mEoS import C1Oleate
from lib.mEoS import C1Linolenate
from lib.mEoS import C1Linoleate
from lib.mEoS import C1Palmitate
from lib.mEoS import C1Stearate

# Heteroatom
from lib.mEoS import Methanol
from lib.mEoS import Ethanol
from lib.mEoS import Acetone
from lib.mEoS import DME
from lib.mEoS import DEE
from lib.mEoS import DMC
from lib.mEoS import NF3
from lib.mEoS import SF6
from lib.mEoS import HCl

# CFCs
from lib.mEoS import CF3I
from lib.mEoS import C4F10
from lib.mEoS import C5F12
from lib.mEoS import R11
from lib.mEoS import R12
from lib.mEoS import R13
from lib.mEoS import R14
from lib.mEoS import R21
from lib.mEoS import R22
from lib.mEoS import R23
from lib.mEoS import R32
from lib.mEoS import R40
from lib.mEoS import R41
from lib.mEoS import R113
from lib.mEoS import R114
from lib.mEoS import R115
from lib.mEoS import R116
from lib.mEoS import R123
from lib.mEoS import R124
from lib.mEoS import R125
from lib.mEoS import R134a
from lib.mEoS import R141b
from lib.mEoS import R142b
from lib.mEoS import R143a
from lib.mEoS import R152a
from lib.mEoS import R161
from lib.mEoS import R218
from lib.mEoS import R227ea
from lib.mEoS import R236ea
from lib.mEoS import R236fa
from lib.mEoS import R245ca
from lib.mEoS import R245fa
from lib.mEoS import R365mfc
from lib.mEoS import RC318
from lib.mEoS import R1234yf
from lib.mEoS import R1234ze
from lib.mEoS import R1216
from lib.mEoS import R1233zd
from lib.mEoS import RE143a
from lib.mEoS import RE245cb2
from lib.mEoS import RE245fa2
from lib.mEoS import RE347mcc

# Siloxanes
from lib.mEoS import D4
from lib.mEoS import D5
from lib.mEoS import D6
from lib.mEoS import MDM
from lib.mEoS import MD2M
from lib.mEoS import MD3M
from lib.mEoS import MD4M
from lib.mEoS import MM

# PseudoCompound
from lib.mEoS import Air
from lib.mEoS import R404a
from lib.mEoS import R407c
from lib.mEoS import R410a
from lib.mEoS import R507a

# TODO: Add Novec 649 from REFPROP

__all__ = MEoS.__subclasses__()[1:] + MEoSBlend.__subclasses__()
id_mEoS = [i.id for i in __all__]


if __name__ == "__main__":
        import doctest
#    import timeit
#    def test():
        for module in __all__:
            if module.__module__ != "R245fa":
                continue
            print(module.__module__)
            inst = module()
            for eq in inst.eq[0:1]:
                if "__test__" in eq:
                    inst.__doc__ += eq["__test__"]
            if inst._viscosity is not None:
                for eq in inst._viscosity:
                    if "__test__" in eq:
                        inst.__doc__ += eq["__test__"]
            if inst._thermal is not None:
                for eq in inst._thermal:
                    if "__test__" in eq:
                        inst.__doc__ += eq["__test__"]
            doctest.run_docstring_examples(inst, globs={module.__module__: module})
#    timeit.timeit("test()", setup="from __main__ import test", number=3)

# TODO: Add 1-propanol from 10.1016_j.fluid.2004.06.028
# TODO: Add 2-propanol from 10.1063/1.3112608
