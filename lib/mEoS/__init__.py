#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS, MEoSBlend

# Noble Gases
import He
import Ne
import Ar
import Kr
import Xe

# Gases
import H2
import D2
import pH2
import oH2
import N2
import O2
import F2
import H2O
import D2O
import CO2
import CO
import N2O
import SO2
import COS
import NH3
import H2S

# Alkanes
import CH4
import C2
import C3
import nC4
import iC4
import nC5
import neoC5
import iC5
import nC6
import iC6
import nC7
import nC8
import nC9
import nC10
import nC12

# Naphthenes
import Cyclopropane
import Cyclopentane
import Cyclohexane
import C1Cyclohexane
import C3Cyclohexane

# Alkenes
import Benzene
import Toluene
import Ethylene
import Propylene
import Butene_1
import iButene
import Cis_2_butene
import Trans_2_butene
import Propyne
import C1Oleate
import C1Linolenate
import C1Linoleate
import C1Palmitate
import C1Stearate

# Heteroatom
import Methanol
import Ethanol
import Acetone
import DME
import DMC
import NF3
import SF6

# CFCs
import CF3I
import C4F10
import C5F12
import R11
import R12
import R13
import R14
import R21
import R22
import R23
import R32
import R41
import R113
import R114
import R115
import R116
import R123
import R124
import R125
import R134a
import R141b
import R142b
import R143a
import R152a
import R161
import R218
import R227ea
import R236ea
import R236fa
import R245ca
import R245fa
import R365mfc
import RC318
import R1234yf
import R1234ze

# Siloxanes
import D4
import D5
import D6
import MDM
import MD2M
import MD3M
import MD4M
import MM

# PseudoCompound
import Air
import R404a
import R407c
import R410a
import R507a

__all__ = MEoS.__subclasses__()[1:] + MEoSBlend.__subclasses__()
id_mEoS = [i.id for i in __all__]


if __name__ == "__main__":
        import doctest
#    import timeit
#    def test():
        for module in __all__:
            if module.__module__ != "nC12":
                continue
            print module.__module__
            inst = module()
            for eq in inst.eq:
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
