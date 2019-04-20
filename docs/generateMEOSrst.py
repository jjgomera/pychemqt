#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Generate the MEoS-*.rst files with list of references

import os

import lib


# Generate index file
txt2 = ""
txt = "lib.mEoS"
txt += os.linesep + "========" + os.linesep + os.linesep
txt += "The list of available fluid with high quality multiparameter equations"
txt += " is automatically updated here:" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in lib.mEoS.__all__:
    txt += "    lib.mEoS.%s" % mod.__name__ + os.linesep
    txt2 += "* :doc:`%s <lib.mEoS.%s>` (%s)" % (mod.__name__, mod.__name__,
                                                mod.name) + os.linesep

with open("docs/lib.mEoS.rst", "w") as file:
    file.write(txt)
with open("docs/lib.mEoSlst.rst", "w") as file:
    file.write(txt2)

# Generate the individual compounds files
for mod in lib.mEoS.__all__:
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
    txt += os.linesep

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
            txt += os.linesep

    if mod._thermal:
        txt += "Thermal Conductivity" + os.linesep
        txt += "--------------------" + os.linesep + os.linesep
        for eq in mod._thermal:
            rf = eq["__doi__"]
            txt += "* %s; %s. %s" % (rf["autor"], rf["title"], rf["ref"])
            if rf["doi"]:
                txt += ", http://dx.doi.org/%s" % rf["doi"]
            txt += os.linesep

    imageFname = "docs/images/%s.png" % mod.__name__
    if os.path.isfile(imageFname):
        txt += os.linesep + "Calculation example" + os.linesep
        txt += "-------------" + os.linesep
        txt += os.linesep + os.linesep
        txt += "Using the first option for equation of state, we can get this "
        txt += "diagram plots with the liquid-gas saturation region:"
        txt += os.linesep + os.linesep
        txt += ".. image:: images/%s.png" % mod.__name__
        txt += os.linesep + os.linesep

    with open("docs/lib.mEoS.%s.rst" % mod.__name__, "w") as file:
        file.write(txt)
