#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Batch editing of intro information in modules
# Run each year for update date and other possible uses


import sys
from pathlib import Path

from tools import firstrun


old = firstrun.__doc__

new = """Pychemqt, Chemical Engineering Process simulator
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


# Find all python files in project folder
files = Path(sys.argv[0]).parent.parent.absolute().glob("**/*.py")

for archivo in files:

    with open(archivo, "r", encoding="utf-8") as readfile:
        code = "".join(readfile.readlines())

    try:
        code.index(old)
    except ValueError:
        print(f"Docs don't found in {archivo}")

    with open(archivo, "w", encoding="utf-8") as writefile:
        writefile.write(code.replace(old, new))
