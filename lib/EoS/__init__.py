#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from PyQt5.QtWidgets import QApplication

from lib.EoS import cubic
from lib.EoS.cubic import alfa
from . import BWRS
from . import BWRSoave
from . import Cubic
from . import Grayson_Streed
from . import Lee_Kesler
from . import virial


K = Cubic._all + BWRS._all + Lee_Kesler._all + Grayson_Streed._all + virial._all
K_name = [k.__title__.split(" (")[0] for k in K]
K_status = [k.__status__ for k in K]
H = Cubic._all + BWRS._all + Lee_Kesler._all
H_name = [h.__title__.split(" (")[0] for h in H]
H_status = [h.__status__ for h in H]

mix = ("van der Waals", "Stryjek-Vera", "Panagiotopoulos", "Melhem")
cp_ideal = (QApplication.translate("pychemqt", "Ideal"), "DIPPR")

__all__ = [BWRS, BWRSoave, Cubic, Grayson_Streed, Lee_Kesler, virial]

# Add references
# each submodule must define its custom __doi__
__doi__ = {}
for obj in __all__:
    if "__doi__" in obj.__dict__:
        __doi__[obj.__name__] = obj.__doi__

    # Add references in modules with only a equation
    if len(obj._all) == 1:
        for obj2 in obj._all:
            if "__doi__" in obj2.__dict__:
                __doi__[obj.__name__] = dict(enumerate(obj2.__doi__))

# Add cubic library
__doi__["lib.EoS.cubic"] = cubic.__doi__


# Equations classification

# Virial-type
# Anderko 1990

# Virial
# BWR
# BWR-Starling-Han
# BWR-Nishiumi

# Cubic-empirical EoS vdW derivatives
# vdW
# RK
# SRK
# PR
# PTV

# Non-cubic EoS vdW derivatives
# Carnahan-Starling
# BACK
# Heiling-Franck
# Dieters
# Soave-quartic

# Wei-Sadus 2000
# Molecular-based EoS, for chain-molecules
# PHCT
# SPHCT
# PACT
# TPT
# PHSC

# MOlecular-based EoS, for associating fluids
# APACT
# SAFT
# SSAFT
# CPA
# AEOS
