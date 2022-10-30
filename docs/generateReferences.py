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


# Generate the *-ref.rst files with list of references

import os

import lib
import UI

total = []
for library in lib.__all__:
    # Import module to intronspection
    __import__("lib.%s" % library)
    module = lib.__getattribute__(library)

    # Special case for mEoS library
    if library == "mEoS":
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
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :private-members:", file=file)
        print("    :show-inheritance:", file=file)
        print("    :member-order: bysource", file=file)

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

# Global references file
with open("docs/references.rst", "w") as file:
    print("References", file=file)
    print("----------", file=file)

    id = 0
    for lnk in sorted(total, key=lambda lnk: str.upper(lnk["autor"])):
        if lnk["autor"] or lnk["title"] or lnk["ref"]:
            id += 1
            ref = "%i. " % id
            if lnk["autor"]:
                ref += "%s; " % lnk["autor"]
            if lnk["title"]:
                ref += "%s. " % lnk["title"]
            if lnk["ref"]:
                ref += "%s" % lnk["ref"]

            if lnk["doi"]:
                ref += ", http://dx.doi.org/%s" % lnk["doi"]

            print(ref, file=file)

# UI module
# Generate index file
txt = "UI package" + os.linesep
txt += "==========" + os.linesep + os.linesep
txt += "Submodules" + os.linesep
txt += "----------" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 2" + os.linesep + os.linesep

for mod in UI.__all__:
    txt += "    UI.%s" % mod + os.linesep

with open("docs/UI.rst", "w") as file:
    file.write(txt)

# Generate each module documentation file
for ui in UI.__all__:
    # Make UI.rst schemas
    with open("docs/UI.%s.rst" % ui, "w") as file:
        print("UI.%s module" % ui, file=file)
        print("="*(len(ui)+4+7), file=file)
        print("", file=file)
        print(".. automodule:: UI.%s" % ui, file=file)
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :private-members:", file=file)
        print("    :show-inheritance:", file=file)
        print("    :member-order: bysource", file=file)
