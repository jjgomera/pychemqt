#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Library to implement a library for chemical elements as mendeleev table
###############################################################################


import os
import sqlite3

from numpy import linspace, logspace, log
from PyQt5.QtCore import QLocale

from lib.utilities import colors


# Connection to database with element data
connection = sqlite3.connect(os.path.join(
    os.environ["pychemqt"], "dat", "elemental.db"))
databank = connection.cursor()

# Load system locale to implement a custon translation system (non qt)
locale = QLocale.system().name().upper()
if "_" in locale:
    locale = locale.split("_")[0]

databank.execute("PRAGMA table_info(TRANSLATION)")
translation = []
for i, name, type_, other, other2, primary_key in databank:
    if "name_" in name:
        translation.append(name.split("_")[-1])
if locale in translation:
    tr_available = True
else:
    tr_available = False


def cleanFloat(flo):
    if flo:
        try:
            value = float(flo)
        except ValueError:
            value = float(flo.split("(")[1].split(",")[0])
    else:
        value = 0
    return value

color_serie = ["#DDDDDD", "#795681", "#B92D2D", "#B8873A", "#D7C848",
               "#94738F", "#6186AC", "#88AE62", "#949692", "#BF924E",
               "#C44343"]
color_phase = ["#DDDDDD", "#BB8F4A", "#7BB245", "#5D82A8"]
NUMERIC_VALUES = ["density_Solid", "density_Liq", "density_Gas", "date",
                  "atomic_mass", "atomic_volume", "atomic_radius",
                  "covalent_radius", "vanderWaals_radius", "electronegativity",
                  "electron_affinity", "first_ionization", "Tf", "Tb",
                  "Heat_f", "Heat_b", "Cp", "k", "T_debye"]


def _configValues(Preferences):
    PROP = Preferences.get("Applications", "elementalColorby")
    NUM = Preferences.getint("Applications", "elementalDefinition")
    LOG = Preferences.getboolean("Applications", "elementalLog")

    PMIN = None
    PMAX = None
    if PROP == "phase":
        CATEGORIES = ["", "Solid", "Liquid", "Gas"]
        COLORS = color_phase
    elif PROP in NUMERIC_VALUES:
        databank.execute("SELECT %s FROM ELEMENTS" % PROP)
        PMAX = 0
        for st, in databank:
            value = cleanFloat(st)
            if value > PMAX:
                PMAX = value

        if LOG:
            PMIN = 1
            CATEGORIES = logspace(log(PMIN), log(PMAX), NUM)
        else:
            PMIN = 0
            CATEGORIES = linspace(PMIN, PMAX, NUM)
        COLORS = colors(NUM, scale=True)
    elif PROP == "Element":
        CATEGORIES = []
        COLORS = []
    else:
        q = "SELECT %s, COUNT(*) c FROM ELEMENTS GROUP BY %s HAVING c > 0" % (
            PROP, PROP)
        databank.execute(q)
        CATEGORIES = []
        for category, count in databank:
            CATEGORIES.append(category)
        if PROP == "serie":
            COLORS = color_serie
        else:
            COLORS = colors(len(CATEGORIES))

    return CATEGORIES, PROP, COLORS, PMAX


class Elemental(object):
    """Element class with data"""

    def __init__(self, id):
        """@param id atomic number of element

        The class implement this properties
           id: atomic number
           name
           altname
           symbol
           serie
           group
           period
           block
           density_Solid
           density_Liq
           density_Gas
           appearance
           date
           country
           discover
           etymology
           atomic_mass
           atomic_volume
           atomic_radius
           covalent_radius
           vanderWaals_radius
           ionic_radii
           lattice_type
           space_group
           lattice_edges
           lattice_angles
           electron_configuration
           oxidation
           electronegativity
           electron_affinity
           first_ionization
           Tf
           Tb
           Heat_f
           Heat_b
           Cp
           k
           T_debye
           color
           notes
        """

        if id > 118:
            id = 118
        databank.execute("SELECT * FROM ELEMENTS WHERE id=='%i'" % id)
        data = databank.fetchone()

        self.id = int(data[0])
        self.altname = data[2]
        self.symbol = data[3]
        self.serie = data[4]
        self.group = int(data[5])
        self.period = int(data[6])
        self.block = data[7]
        self.density_Solid = self._unit(data[8])
        self.density_Liq = self._unit(data[9])
        self.density_Gas = self._unit(data[10])
        self.appearance = data[11]
        self.date = data[12]
        self.country = data[13]
        self.discover = data[14]
        self.etymology = data[15]
        self.atomic_mass = self._unit(data[16])
        self.atomic_volume = self._unit(data[17])
        self.atomic_radius = self._unit(data[18])
        self.covalent_radius = self._unit(data[19])
        self.vanderWaals_radius = self._unit(data[20])
        self.ionic_radii = data[21]
        self.lattice_type = data[22]
        self.space_group = data[23]
        self.lattice_edges = eval(data[24])
        self.lattice_volume = self.lattice_edges[0]*self.lattice_edges[1] * \
            self.lattice_edges[2] / 1e9
        self.lattice_angles = eval(data[25])
        self.electron_configuration = data[26]
        self.oxidation = data[27]
        self.electronegativity = self._unit(data[28])
        self.electron_affinity = self._unit(data[29])
        self.first_ionization = self._unit(data[30])
        self.Tf = self._unit(data[31])
        self.Tb = self._unit(data[32])
        if not self.Tf or not self.Tb:
            self.phase = ""
        elif self.Tf > 273.15:
            self.phase = "Solid"
        elif self.Tb < 273.15:
            self.phase = "Gas"
        else:
            self.phase = "Liquid"
        self.Heat_f = self._unit(data[33])
        self.Heat_b = self._unit(data[34])
        self.Cp = self._unit(data[35])
        self.k = self._unit(data[36])
        self.T_debye = self._unit(data[37])
        self.color = data[38]
        self.notes = data[39]

        # Translation
        self.name = data[1]
        if tr_available:
            qu = "SELECT name_%s FROM TRANSLATION WHERE id==%i" % (locale, id)
            databank.execute(qu)
            tr_name = databank.fetchone()[0]
            if tr_name:
                self.name = tr_name

        # Isotopes
        query = "SELECT * FROM ISOTOPES WHERE atomic_number==?" + \
            "ORDER BY mass_number"
        databank.execute(query, (self.id, ))
        self.isotopes = []
        for data in databank:
            self.isotopes.append((int(data[4]), data[2], data[3]))

    def _unit(self, str):
        aproximate = False
        try:
            value = float(str)
        except:
            if not str:
                value = None
            elif str[-1] == ")":
                value = float(str.split("(")[1].split(",")[0])
                aproximate = True
        if aproximate:
            value.code = "stimated"
        return value
