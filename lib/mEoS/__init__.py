#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS, MEoSBlend

# Noble Gases
from . import He
from . import Ne
from . import Ar
from . import Kr
from . import Xe

# Gases
from . import H2
from . import D2
from . import pH2
from . import oH2
from . import N2
from . import O2
from . import F2
from . import H2O
from . import D2O
from . import CO2
from . import CO
from . import N2O
from . import SO2
from . import COS
from . import NH3
from . import H2S

# Alkanes

from . import CH4
from . import C2
from . import C3
from . import nC4
from . import iC4
from . import nC5
from . import neoC5
from . import iC5
from . import nC6
from . import iC6
from . import nC7
from . import nC8
from . import iC8
from . import nC9
from . import nC10
from . import nC11
from . import nC12

# Naphthenes
from . import Cyclopropane
from . import Cyclopentane
from . import Cyclohexane
from . import C1Cyclohexane
from . import C3Cyclohexane

# Alkenes
from . import Benzene
from . import Toluene
from . import oXylene
from . import pXylene
from . import pXylene
from . import EthylBenzene
from . import Ethylene
from . import Propylene
from . import Butene_1
from . import iButene
from . import Cis_2_butene
from . import Trans_2_butene
from . import Propyne
from . import C1Oleate
from . import C1Linolenate
from . import C1Linoleate
from . import C1Palmitate
from . import C1Stearate

# Heteroatom
from . import Methanol
from . import Ethanol
from . import Acetone
from . import DME
from . import DEE
from . import DMC
from . import NF3
from . import SF6
from . import HCl

# CFCs
from . import CF3I
from . import C4F10
from . import C5F12
from . import R11
from . import R12
from . import R13
from . import R14
from . import R21
from . import R22
from . import R23
from . import R32
from . import R40
from . import R41
from . import R113
from . import R114
from . import R115
from . import R116
from . import R123
from . import R124
from . import R125
from . import R134a
from . import R141b
from . import R142b
from . import R143a
from . import R152a
from . import R161
from . import R218
from . import R227ea
from . import R236ea
from . import R236fa
from . import R245ca
from . import R245fa
from . import R365mfc
from . import RC318
from . import R1234yf
from . import R1234ze
from . import R1216
from . import R1233zd
from . import RE143a
from . import RE245cb2
from . import RE245fa2
from . import RE347mcc

# Siloxanes
from . import D4
from . import D5
from . import D6
from . import MDM
from . import MD2M
from . import MD3M
from . import MD4M
from . import MM

# PseudoCompound
from . import Air
from . import R404a
from . import R407c
from . import R410a
from . import R507a

# TODO: Add Novec 649 from REFPROP

__all__ = MEoS.__subclasses__()[1:] + MEoSBlend.__subclasses__()
id_mEoS = [i.id for i in __all__]


if __name__ == "__main__":
        import doctest
#    import timeit
#    def test():
        for module in __all__:
            if module.__module__ != "D2O":
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
