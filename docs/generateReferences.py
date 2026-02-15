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


# Generate the *-ref.rst files with list of references

import os

import equipment
import lib
import plots
import tools
import UI

# List with all references
total = []

# Lib module
# Generate index file
txt = "lib package" + os.linesep
txt += "===========" + os.linesep + os.linesep
txt += "Submodules" + os.linesep
txt += "----------" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in lib.__all__:
    txt += "    lib.%s" % mod + os.linesep

with open("docs/lib.rst", "w") as file:
    file.write(txt)


for library in lib.__all__:
    # Import module to intronspection
    __import__("lib.%s" % library)
    module = lib.__getattribute__(library)

    # Special case for mEoS library
    if library in ("mEoS", "newComponent"):
        for cmp in module.__doi__:
            for eq in module.__doi__[cmp]:
                rf = module.__doi__[cmp][eq]
                if rf not in total:
                    total.append(rf)
        continue

    # Special case for EoS library
    if library == "EoS":
        for cmp in module.__doi__:
            if cmp == "lib.EoS.Cubic":
                for eq in module.__doi__[cmp]:
                    lst = module.__doi__[cmp][eq]
                    for rf in lst:
                        if rf not in total:
                            total.append(rf)
            else:
                for eq in module.__doi__[cmp]:
                    rf = module.__doi__[cmp][eq]
                    if rf not in total:
                        total.append(rf)
        continue

    # General case for simple libraries
    # Make lib.rst schemas
    with open("docs/lib.%s.rst" % library, "w") as file:
        print("lib.%s module" % library, file=file)
        print("="*(len(library)+4+7), file=file)
        print("", file=file)
        print(".. automodule:: lib.%s" % library, file=file)

        if hasattr(module, "__doi__") and module.__doi__:
            print("", file=file)
            print(".. include:: lib.%s_ref.rst" % library, file=file)

    if hasattr(module, "__doi__") and module.__doi__:
        with open("docs/lib.%s_ref.rst" % library, "w") as file:
            print("References", file=file)
            print("----------", file=file)
            for id, rf in module.__doi__.items():
                id = str(id)
                if library == "iapws97" and id != "1":
                    print(" * [%s] %s; %s. %s" % (
                        id, rf["autor"], rf["title"], rf["ref"]), file=file)
                else:
                    print(".. [%s] %s; %s. %s" % (
                        id, rf["autor"], rf["title"], rf["ref"]), file=file)
                if rf not in total:
                    total.append(rf)

# Plots module
# Generate index file
txt = "plots package" + os.linesep
txt += "=============" + os.linesep + os.linesep
txt += "Submodules" + os.linesep
txt += "----------" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in plots.__all__:
    txt += "    plots.%s" % mod + os.linesep

with open("docs/plots.rst", "w") as file:
    file.write(txt)

# Generate each module documentation file
for lib in plots.__all__:
    # Make plots.rst schemas
    with open("docs/plots.%s.rst" % lib, "w") as file:
        print("plots.%s module" % lib, file=file)
        print("="*(len(lib)+6+7), file=file)
        print("", file=file)
        print(".. automodule:: plots.%s" % lib, file=file)


# Tools module
# Generate index file
txt = "tools package" + os.linesep
txt += "=============" + os.linesep + os.linesep
txt += "Submodules" + os.linesep
txt += "----------" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in tools.__all__:
    txt += "    tools.%s" % mod + os.linesep

with open("docs/tools.rst", "w") as file:
    file.write(txt)

# Generate each module documentation file
for tool in tools.__all__:
    # Make tools.rst schemas
    with open("docs/tools.%s.rst" % tool, "w") as file:
        print("tools.%s module" % tool, file=file)
        print("="*(len(tool)+6+7), file=file)
        print("", file=file)
        print(".. automodule:: tools.%s" % tool, file=file)

        if tool == "UI_Tables":
            for sub in ("chooseFluid", "library", "plot", "reference", "table"):
                print(".. automodule:: tools.%s.%s" % (tool, sub), file=file)


# UI module
# Generate index file
txt = "UI package" + os.linesep
txt += "==========" + os.linesep + os.linesep
txt += "Submodules" + os.linesep
txt += "----------" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in UI.__all__:
    txt += "    UI.%s" % mod + os.linesep

with open("docs/UI.rst", "w") as file:
    file.write(txt)

# Generate each module documentation file
for lib in UI.__all__:
    # Make UI.rst schemas
    with open("docs/UI.%s.rst" % lib, "w") as file:
        print("UI.%s module" % lib, file=file)
        print("="*(len(lib)+3+7), file=file)
        print("", file=file)
        print(".. automodule:: UI.%s" % lib, file=file)


# Equipment module
# Generate index file
txt = "equipment package" + os.linesep
txt += "=================" + os.linesep + os.linesep
txt += "Submodules" + os.linesep
txt += "----------" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in equipment.__all__:
    txt += "    equipment.%s" % mod + os.linesep
for mod in equipment.widget.__all__:
    txt += "    equipment.widget.%s" % mod + os.linesep

with open("docs/equipment.rst", "w") as file:
    file.write(txt)

# Generate each module documentation file
for equip in equipment.equipments:
    # Make equipment.rst schemas
    with open("docs/equipment.%s.rst" % equip.__name__, "w") as file:
        print("equipment.%s module" % equip.__name__, file=file)
        print("="*(len(equip.__name__)+10+7), file=file)
        print("", file=file)
        print(".. autoclass:: equipment.%s" % equip.__name__, file=file)

        if equip.__doi__:
            print("", file=file)
            print(".. include:: equipment.%s_ref.rst" % equip.__name__, file=file)

    if equip.__doi__:
        with open("docs/equipment.%s_ref.rst" % equip.__name__, "w") as file:
            print("References", file=file)
            print("----------", file=file)
            for id, rf in enumerate(equip.__doi__):
                id = str(id+1)
                print(".. [%s] %s; %s. %s" % (
                    id, rf["autor"], rf["title"], rf["ref"]), file=file)

                if rf not in total:
                    total.append(rf)

# Generate each module documentation file
for wdg in equipment.widget.__all__:

    __import__("equipment.widget.%s" % wdg)
    module = equipment.widget.__getattribute__(wdg)

    # Make equipment.widget.rst schemas
    with open("docs/equipment.widget.%s.rst" % wdg, "w") as file:
        print("equipment.widget.%s module" % wdg, file=file)
        print("="*(len(wdg)+10+13), file=file)
        print("", file=file)
        print(".. automodule:: equipment.widget.%s" % wdg, file=file)

        if module.__doi__:
            print("", file=file)
            print(".. include:: equipment.widget.%s_ref.rst" % wdg, file=file)

    if module.__doi__:
        with open("docs/equipment.widget.%s_ref.rst" % wdg, "w") as file:
            print("References", file=file)
            print("----------", file=file)
            for id, rf in module.__doi__.items():
                print(".. [%s] %s; %s. %s" % (
                    id, rf["autor"], rf["title"], rf["ref"]), file=file)

                if rf not in total:
                    total.append(rf)


# Generate global references file
with open("docs/references.rst", "w") as file:
    print("References", file=file)
    print("----------", file=file)

    id = 0
    for lnk in sorted(total, key=lambda lnk: str.upper(lnk["autor"])
                      if lnk["autor"] else str.upper(lnk["title"])):
        if lnk["autor"] or lnk["title"] or lnk["ref"]:
            id += 1
            ref = f"{id}. "
            if lnk["autor"]:
                ref += f"{lnk['autor']}; "
            if lnk["title"]:
                ref += f"{lnk['title']}. "
            if lnk["ref"]:
                ref += f"{lnk['ref']}"

            if lnk["doi"]:
                ref += f", http://dx.doi.org/{lnk['doi']}"

            print(ref, file=file)
