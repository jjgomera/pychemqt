#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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



# Generate the MEoS-*.rst files with list of references

import os

from lib.mEoS import __all__


# Generate index file
txt2 = ""
txt = "lib.mEoS"
txt += os.linesep + "========" + os.linesep + os.linesep
txt += "The list of available fluid with high quality multiparameter equations"
txt += " is automatically updated here:" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in __all__:
    txt += "    lib.mEoS.%s" % mod.__name__ + os.linesep
    txt2 += "* :doc:`%s <lib.mEoS.%s>` (%s)" % (mod.__name__, mod.__name__,
                                                mod.name) + os.linesep

with open("docs/lib.mEoS.rst", "w") as file:
    file.write(txt)
with open("docs/lib.mEoSlst.rst", "w") as file:
    file.write(txt2)

# Generate the individual compounds files
for mod in __all__:
    txt = "lib.mEoS.%s" % mod.__name__
    txt += os.linesep + "="*len(txt) + os.linesep + os.linesep

    # Section for fluid info
    txt += "Fluid info" + os.linesep
    txt += "----------" + os.linesep + os.linesep
    txt += "* CAS Number: " + mod.CASNumber + os.linesep
    txt += "* Formula: " + mod.formula + os.linesep
    if mod.synonym:
        txt += "* Synonym: " + mod.synonym + os.linesep

    txt += "* Molecular weigth: %s g/mol" % mod.M + os.linesep
    txt += "* Tc: " + mod.Tc.str + os.linesep
    txt += "* Pc: " + mod.Pc.get_str("MPa") + os.linesep
    txt += "* ρc: " + mod.rhoc.str + os.linesep
    txt += "* Tt: " + mod.Tt.str + os.linesep
    txt += "* Tb: " + mod.Tb.str + os.linesep
    txt += "* Acentric factor: %s" % mod.f_acent + os.linesep
    txt += "* Dipole moment: %s" % mod.momentoDipolar.get_str("Debye")
    txt += os.linesep + os.linesep

    txt += "Equation of state" + os.linesep
    txt += "-----------------" + os.linesep + os.linesep
    for eq in mod.eq:
        rf = eq["__doi__"]
        txt += "* %s; %s. %s" % (rf["autor"], rf["title"], rf["ref"])
        if rf["doi"]:
            txt += ", http://dx.doi.org/%s" % rf["doi"]
        txt += os.linesep
    txt += os.linesep + os.linesep

    if mod._viscosity:
        txt += "Viscosity" + os.linesep
        txt += "---------" + os.linesep + os.linesep
        for eq in mod._viscosity:
            rf = eq["__doi__"]
            txt += "* %s; %s. %s" % (rf["autor"], rf["title"], rf["ref"])
            if rf["doi"]:
                txt += ", http://dx.doi.org/%s" % rf["doi"]
            txt += os.linesep + os.linesep

    if mod._thermal:
        txt += "Thermal Conductivity" + os.linesep
        txt += "--------------------" + os.linesep + os.linesep
        for eq in mod._thermal:
            rf = eq["__doi__"]
            txt += "* %s; %s. %s" % (rf["autor"], rf["title"], rf["ref"])
            if rf["doi"]:
                txt += ", http://dx.doi.org/%s" % rf["doi"]
            txt += os.linesep + os.linesep

    imageFname = "docs/images/%s.png" % mod.__name__
    if os.path.isfile(imageFname):
        txt += os.linesep + "Calculation example" + os.linesep
        txt += "-------------------"
        txt += os.linesep + os.linesep
        txt += "Using the first option for equation of state, we can get this "
        txt += "diagram plots with the liquid-gas saturation region:"
        txt += os.linesep + os.linesep
        txt += ".. image:: images/%s.png" % mod.__name__
        txt += os.linesep + os.linesep
        txt += "The diagram is generated with "
        txt += ":download:`this module <../.script/plotMEoS.py>` "
        txt += "running with the compound name as parameter or edited in file"
        txt += os.linesep + os.linesep
        txt += ".. code-block:: bash" + os.linesep + os.linesep
        txt += "    python3 plotMEoS.py %s" % mod.__name__
        txt += os.linesep + os.linesep

    with open("docs/lib.mEoS.%s.rst" % mod.__name__, "w") as file:
        file.write(txt)
