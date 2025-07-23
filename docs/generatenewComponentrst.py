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


# Generate the EoS-*.rst files with list of references

import os

from lib.newComponent import __all__, _methods


# Generate index file
txt = "lib.newComponent module" + os.linesep
txt += "=======================" + os.linesep + os.linesep
txt += ".. include:: newComponent.rst" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in _methods:
    txt += f"    lib.newComponent.{mod.__name__}" + os.linesep

with open("docs/lib.newComponent.rst", "w") as file:
    file.write(txt)

# Generate documentation file for _base file
with open("docs/lib.newComponent._base.rst", "w") as file:
    print("lib.newComponent._base module", file=file)
    print("=============================", file=file)
    print("", file=file)
    print(".. automodule:: lib.newComponent._base", file=file)
    print("    :members:", file=file)
    print("    :undoc-members:", file=file)
    print("    :show-inheritance:", file=file)
    print("    :member-order: bysource", file=file)

# Generate files for each methods
# for mod in _methods:
for fname, mod in zip(__all__, _methods):

    library = mod.__name__
    with open(f"docs/lib.newComponent.{library}.rst", "w") as file:
        print(f"{library}", file=file)
        print("="*(len(library)+1), file=file)
        print("", file=file)
        print(f".. autoclass:: lib.newComponent.{fname}.{library}", file=file)
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :show-inheritance:", file=file)
        print("    :member-order: bysource", file=file)

        try:
            doi = mod.__doi__
        except AttributeError:
            continue

        print("", file=file)
        print("References", file=file)
        print("----------", file=file)

        count = 1

        for lnk in doi:
            ref = f".. [{count}] {lnk['autor']}; {lnk['title']}, {lnk['ref']}"
            if lnk["doi"]:
                ref += f", http://dx.doi.org/{lnk['doi']}"
            count += 1
            print(ref, file=file)
