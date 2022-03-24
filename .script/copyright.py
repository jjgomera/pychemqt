#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Batch editing of intro information in modules
# Run each year for update date and other possible uses


from glob import glob
from tools import firstrun


old = firstrun.__doc__

files = glob("/home/jjgomera/Programacion/pychemqt/**/*.py", recursive=True)
files.append("/home/jjgomera/Programacion/pychemqt/pychemqt.py")

new = """Pychemqt, Chemical Engineering Process simulator
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


for archivo in files:

    with open(archivo, "r") as readfile:
        code = "".join(readfile.readlines())

    try:
        code.index(old)
    except ValueError:
        print("Error in %s" % archivo)

    with open(archivo, "w") as writefile:
        writefile.write(code.replace(old, new))
