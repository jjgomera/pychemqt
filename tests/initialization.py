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


from configparser import ConfigParser
import os
import shutil
import sys


# Add pychemqt folder to python path
sys.path.insert(0, os.path.abspath('.'))

# Define pychemqt environment
os.environ["pychemqt"] = os.path.abspath('.')

conf_dir = os.path.expanduser("~") + os.sep + ".pychemqt" + os.sep

# Checking config folder
if not os.path.isdir(conf_dir):
    os.mkdir(conf_dir)

try:
    open(conf_dir + "pychemqt.log", 'x')
except FileExistsError:  # noqa
    pass


# Checking config files
from tools import firstrun  # noqa

# Checking config file
default_Preferences = firstrun.Preferences()
change = False
if not os.path.isfile(conf_dir + "pychemqtrc"):
    default_Preferences.write(open(conf_dir + "pychemqtrc", "w"))
    Preferences = default_Preferences
    change = True
else:
    # Check Preferences options to find set new options
    Preferences = ConfigParser()
    Preferences.read(conf_dir + "pychemqtrc")
    for section in default_Preferences.sections():
        if not Preferences.has_section(section):
            Preferences.add_section(section)
            change = True
        for option in default_Preferences.options(section):
            if not Preferences.has_option(section, option):
                value = default_Preferences.get(section, option)
                Preferences.set(section, option, value)
                change = True
    if change:
        Preferences.write(open(conf_dir + "pychemqtrc", "w"))

# FIXME: This file might not to be useful but for now I use it to save project
# configuration data
if not os.path.isfile(conf_dir + "pychemqtrc_temporal"):
    Config = firstrun.config()
    Config.write(open(conf_dir + "pychemqtrc_temporal", "w"))

# Checking costindex
if not os.path.isfile(conf_dir + "CostIndex.dat"):
        orig = os.path.join(os.environ["pychemqt"], "dat", "costindex.dat")
        with open(orig) as cost_index:
            lista = cost_index.readlines()[-1].split(" ")
            with open(conf_dir + "CostIndex.dat", "w") as archivo:
                for data in lista:
                    archivo.write(data.replace(os.linesep, "") + os.linesep)

# Checking currency rates
origen = os.path.join(os.environ["pychemqt"], "dat", "moneda.dat")
shutil.copy(origen, conf_dir + "moneda.dat")

# Checking database with custom components
if not os.path.isfile(conf_dir + "databank.db"):
    firstrun.createDatabase(conf_dir + "databank.db")
