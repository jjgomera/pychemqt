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
# Library to define functions necessary to run at first pychemqt run, so don't
# let import other pychemqt library to avoid import error when no config
# files availables.
#   -createDatabase: Create empty database
###############################################################################


from configparser import ConfigParser
import json
import sqlite3
import sys
import urllib.request

# It must be defined previously to avoid to early import of libraries
# See end of lib/unidades.py to know how to get this list, check when new
# magnitude are added
magnitudes = [
    'Acceleration', 'Angle', 'Area', 'CakeResistance', 'Currency', 'Density',
    'DenLiq', 'DenGas', 'DensityPressure', 'DensityTemperature', 'Diffusivity',
    'KViscosity', 'DipoleMoment', 'PotencialElectric', 'Energy', 'Work',
    'Enthalpy', 'EnthalpyDensity', 'EnthalpyPressure', 'Entropy', 'Force',
    'Fouling', 'Frequency', 'V2V', 'HeatFlux', 'HeatTransfCoef', 'Length',
    'ParticleDiameter', 'Thickness', 'PipeDiameter', 'Head', 'Mass',
    'MassFlow', 'Mol', 'MolarDensity', 'MolarEnthalpy', 'MolarFlow',
    'MolarSpecificHeat', 'MolarVolume', 'PackingDP', 'EnergyFlow', 'Power',
    'Pressure', 'DeltaP', 'InvPressure', 'PressureTemperature',
    'PressureDensity', 'SolubilityParameter', 'SpecificHeat',
    'SpecificEntropy', 'SpecificVolume', 'Speed', 'Tension', 'Temperature',
    'DeltaT', 'InvTemperature', 'TemperaturePressure', 'ThermalConductivity',
    'SpecificVolume_square', 'Time', 'UA', 'Viscosity', 'Volume', 'VolLiq',
    'VolGas', 'VolFlow', 'QLiq', 'QGas', 'Dimensionless']

# See end of equipment.__init__.py to know how to get this list, check when new
# fully functional are added
equipos = ['Divider', 'Valve', 'Mixer', 'Pump', 'Compressor', 'Turbine',
           'Pipe', 'Flash', 'ColumnFUG', 'Heat_Exchanger', 'Shell_Tube',
           'Hairpin', 'Fired_Heater', 'Ciclon', 'GravityChamber', 'Baghouse',
           'ElectricPrecipitator', 'Dryer', 'Scrubber', 'Spreadsheet',
           'Reactor']


def which(program):
    """Function to detect program availability in system and return path"""
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return ""


calculator = ""
if sys.platform == "win32":
    calculator = "calc.exe"
else:
    for programa in ["qalculate", "gcalctool", "kcalc"]:
        ejecutable = which(programa)
        if ejecutable:
            calculator = ejecutable
            break


editor = ""
if sys.platform == "win32":
    editor = "notepad.exe"
else:
    for programa in ["gedit", "leafpad", "geany", "kate", "kwrite", "vim",
                     "vi", "emacs", "nano", "pico"]:
        ejecutable = which(programa)
        if ejecutable:
            editor = ejecutable
            break

shell = ""
if sys.platform == "win32":
    shell = ""
else:
    shell = which("xterm")
# TODO: De momento solo soporta xterm
#    for programa in ["xterm", "gnome-terminal", "kterminal", "lxterminal", "xfce4-terminal", "terminator"]:
#        ejecutable=which(programa)
#        if ejecutable:
#            shell=ejecutable
#            break


def Preferences():
    """Function to define a first preferences file"""
    config = ConfigParser()

    # General
    config.add_section("General")
    config.set("General", "Color_Resaltado", "#ffff00")
    config.set("General", "Color_ReadOnly", "#eaeaea")
    config.set("General", "Recent_Files", "10")
    config.set("General", "Load_Last_Project", "True")
    config.set("General", "Tray", "False")

    # PFD
    config.add_section("PFD")
    config.set("PFD", "x", "800")
    config.set("PFD", "y", "600")
    config.set("PFD", "Color_Entrada", "#c80000")
    config.set("PFD", "Color_Salida", "#0000c8")
    config.set("PFD", "Color_Stream", "#000000")
    config.set("PFD", "Width", "1.0")
    config.set("PFD", "Union", "0")
    config.set("PFD", "Miter_limit", "2.0")
    config.set("PFD", "Punta", "0")
    config.set("PFD", "Guion", "0")
    config.set("PFD", "Dash_offset", "0.0")

    # Tooltip
    config.add_section("Tooltip")
    config.set("Tooltip", "Show", "True")
    config.set("Tooltip", "SI", "False")
    config.set("Tooltip", "CGS", "False")
    config.set("Tooltip", "AltSI", "False")
    config.set("Tooltip", "English", "False")
    config.set("Tooltip", "Metric", "False")
    for i, magnitud in enumerate(magnitudes[:-1]):
        config.set("Tooltip", magnitud, "[0,1]")

    # TooltipEntity
    config.add_section("TooltipEntity")
    config.set("TooltipEntity", "Corriente", "[0,1]")
    for equipo in equipos:
        config.set("TooltipEntity", equipo, "[0,1]")

    # NumericFactor
    config.add_section("NumericFormat")
    for magnitud in magnitudes:
        kwarg = {'total': 0, 'signo': False, 'decimales': 4, 'format': 0}
        config.set("NumericFormat", magnitud, str(kwarg))

    # Petro
    config.add_section("petro")
    config.set("petro", "M", "0")
    config.set("petro", "critical", "0")
    config.set("petro", "vc", "0")
    config.set("petro", "f_acent", "0")
    config.set("petro", "Tb", "0")
    config.set("petro", "SG", "0")
    config.set("petro", "n", "0")
    config.set("petro", "Zc", "0")
    config.set("petro", "PNA", "0")
    config.set("petro", "H", "0")
    config.set("petro", "curve", "0")

    # Applications
    config.add_section("Applications")
    config.set("Applications", "Calculator", calculator)
    config.set("Applications", "TextViewer", editor)
    config.set("Applications", "Shell", shell)
    config.set("Applications", "ipython", "False")
    config.set("Applications", "maximized", "False")
    config.set("Applications", "foregroundColor", "#ffffff")
    config.set("Applications", "backgroundColor", "#000000")
    config.set("Applications", "elementalColorby", "serie")
    config.set("Applications", "elementalDefinition", "10")
    config.set("Applications", "elementalLog", "False")

    # mEoS
    config.add_section("MEOS")
    config.set("MEOS", "coolprop", "False")
    config.set("MEOS", "refprop", "False")
    config.set("MEOS", "saturation"+"Color", "#000000")
    config.set("MEOS", "saturation"+"alpha", "255")
    config.set("MEOS", "saturation"+"lineWidth", "1.0")
    config.set("MEOS", "saturation"+"lineStyle", "-")
    config.set("MEOS", "saturation"+"marker", "None")
    config.set("MEOS", "saturation"+"markersize", "3")
    config.set("MEOS", "saturation"+"markerfacecolor", "#ff0000")
    config.set("MEOS", "saturation"+"markeredgewidth", "1")
    config.set("MEOS", "saturation"+"markeredgecolor", "#000000")
    config.set("MEOS", "grid", "False")
    config.set("MEOS", "definition", "1")
    lineas = ["Isotherm", "Isobar", "Isoenthalpic", "Isoentropic", "Isochor",
              "Isoquality"]
    for linea in lineas:
        config.set("MEOS", linea+"Start", "0")
        config.set("MEOS", linea+"End", "0")
        config.set("MEOS", linea+"Step", "0")
        config.set("MEOS", linea+"Custom", "True")
        config.set("MEOS", linea+"List", "")
        if linea != "Isoquality":
            config.set("MEOS", linea+"Critic", "True")
        config.set("MEOS", linea+"Color", "#000000")
        config.set("MEOS", linea+"alpha", "255")
        config.set("MEOS", linea+"lineWidth", "0.5")
        config.set("MEOS", linea+"lineStyle", "-")
        config.set("MEOS", linea+"marker", "None")
        config.set("MEOS", linea+"markersize", "3")
        config.set("MEOS", linea+"markerfacecolor", "#ff0000")
        config.set("MEOS", linea+"markeredgewidth", "1")
        config.set("MEOS", linea+"markeredgecolor", "#000000")

        config.set("MEOS", linea+"Label", "False")
        config.set("MEOS", linea+"Variable", "False")
        config.set("MEOS", linea+"Units", "False")
        config.set("MEOS", linea+"Position", "50")

    # Psychr
    config.add_section("Psychr")
    config.set("Psychr", "chart", "True")
    config.set("Psychr", "virial", "False")
    config.set("Psychr", "coolprop", "False")
    config.set("Psychr", "refprop", "False")

    config.set("Psychr", "saturation"+"Color", "#000000")
    config.set("Psychr", "saturation"+"alpha", "255")
    config.set("Psychr", "saturation"+"lineWidth", "0.5")
    config.set("Psychr", "saturation"+"lineStyle", "-")
    config.set("Psychr", "saturation"+"marker", "None")
    config.set("Psychr", "saturation"+"markersize", "3")
    config.set("Psychr", "saturation"+"markerfacecolor", "#ff0000")
    config.set("Psychr", "saturation"+"markeredgewidth", "1")
    config.set("Psychr", "saturation"+"markeredgecolor", "#000000")

    config.set("Psychr", "crux"+"Color", "#0000ff")
    config.set("Psychr", "crux"+"alpha", "255")
    config.set("Psychr", "crux"+"lineWidth", "0.5")
    config.set("Psychr", "crux"+"lineStyle", "-")
    config.set("Psychr", "crux"+"marker", "None")
    config.set("Psychr", "crux"+"markersize", "3")
    config.set("Psychr", "crux"+"markerfacecolor", "#ff0000")
    config.set("Psychr", "crux"+"markeredgewidth", "1")
    config.set("Psychr", "crux"+"markeredgecolor", "#000000")

    lineas = ["IsoTdb", "IsoW", "IsoHR", "IsoTwb", "Isochor"]
    values = [
        {"start": 274.0, "end": 330.0, "step": 1.0, "color": "#000000",
         "linewidth": 0.5, "linestyle": ":", "label": "False",
         "units": "False", "position": 50},
        {"start": 0.0, "end": 0.04, "step": 0.001, "color": "#000000",
         "linewidth": 0.5, "linestyle": ":", "label": "False",
         "units": "False", "position": 50},
        {"start": 10.0, "end": 100.0, "step": 10.0, "color": "#000000",
         "linewidth": 0.5, "linestyle": "--", "label": "True",
         "units": "True", "position": 85},
        {"start": 250.0, "end": 320.0, "step": 1.0, "color": "#aa0000",
         "linewidth": 0.8, "linestyle": ":", "label": "False",
         "units": "False", "position": 90},
        {"start": 0.8, "end": 1.0, "step": 0.01, "color": "#00aa00",
         "linewidth": 0.8, "linestyle": ":", "label": "False",
         "units": "False", "position": 90}]
    for linea, value in zip(lineas, values):
        config.set("Psychr", linea+"Start", str(value["start"]))
        config.set("Psychr", linea+"End", str(value["end"]))
        config.set("Psychr", linea+"Step", str(value["step"]))
        config.set("Psychr", linea+"Custom", "False")
        config.set("Psychr", linea+"List", "")
        config.set("Psychr", linea+"Color", str(value["color"]))
        config.set("Psychr", linea+"alpha", "255")
        config.set("Psychr", linea+"lineWidth", str(value["linewidth"]))
        config.set("Psychr", linea+"lineStyle", str(value["linestyle"]))
        config.set("Psychr", linea+"marker", "None")
        config.set("Psychr", linea+"markersize", "3")
        config.set("Psychr", linea+"markerfacecolor", "#ff0000")
        config.set("Psychr", linea+"markeredgewidth", "1")
        config.set("Psychr", linea+"markeredgecolor", "#000000")
        config.set("Psychr", linea+"Label", str(value["label"]))
        config.set("Psychr", linea+"Units", str(value["units"]))
        config.set("Psychr", linea+"Position", str(value["position"]))
        config.set("Psychr", linea+"variable", str(False))

    # Moody
    config.add_section("Moody")
    config.set("Moody", "fanning", "False")
    config.set("Moody", "method", "0")
    config.set("Moody", "ed", "[0, 1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, "
               "4e-4, 6e-4, 8e-4, .001, .0015, .002, .003, .004, .006, .008, "
               ".01, .0125, .015, .0175, .02, .025, .03, .035, .04, .045, "
               ".05, .06, .07]")

    config.set("Moody", "line"+"Color", "#000000")
    config.set("Moody", "line"+"alpha", "255")
    config.set("Moody", "line"+"lineWidth", "0.5")
    config.set("Moody", "line"+"lineStyle", "-")
    config.set("Moody", "line"+"marker", "None")
    config.set("Moody", "line"+"markersize", "3")
    config.set("Moody", "line"+"markerfacecolor", "#ff0000")
    config.set("Moody", "line"+"markeredgewidth", "1")
    config.set("Moody", "line"+"markeredgecolor", "#000000")

    config.set("Moody", "crux"+"Color", "#0000ff")
    config.set("Moody", "crux"+"alpha", "255")
    config.set("Moody", "crux"+"lineWidth", "0.5")
    config.set("Moody", "crux"+"lineStyle", "-")
    config.set("Moody", "crux"+"marker", "None")
    config.set("Moody", "crux"+"markersize", "3")
    config.set("Moody", "crux"+"markerfacecolor", "#ff0000")
    config.set("Moody", "crux"+"markeredgewidth", "1")
    config.set("Moody", "crux"+"markeredgecolor", "#000000")

    # Standing-Katz
    config.add_section("Standing_Katz")
    config.set("Standing_Katz", "method", "0")
    config.set("Standing_Katz", "Prmin", "0.8")
    config.set("Standing_Katz", "Prmax", "10")
    config.set("Standing_Katz", "Tr", "[1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35,"
               "1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.]")

    config.set("Standing_Katz", "line"+"Color", "#000000")
    config.set("Standing_Katz", "line"+"alpha", "255")
    config.set("Standing_Katz", "line"+"lineWidth", "0.5")
    config.set("Standing_Katz", "line"+"lineStyle", "-")
    config.set("Standing_Katz", "line"+"marker", "None")
    config.set("Standing_Katz", "line"+"markersize", "3")
    config.set("Standing_Katz", "line"+"markerfacecolor", "#ff0000")
    config.set("Standing_Katz", "line"+"markeredgewidth", "1")
    config.set("Standing_Katz", "line"+"markeredgecolor", "#000000")

    config.set("Standing_Katz", "crux"+"Color", "#0000ff")
    config.set("Standing_Katz", "crux"+"alpha", "255")
    config.set("Standing_Katz", "crux"+"lineWidth", "0.5")
    config.set("Standing_Katz", "crux"+"lineStyle", "-")
    config.set("Standing_Katz", "crux"+"marker", "None")
    config.set("Standing_Katz", "crux"+"markersize", "3")
    config.set("Standing_Katz", "crux"+"markerfacecolor", "#ff0000")
    config.set("Standing_Katz", "crux"+"markeredgewidth", "1")
    config.set("Standing_Katz", "crux"+"markeredgecolor", "#000000")

    # Openbabel
    config.add_section("Openbabel")
    config.set("Openbabel", "BondColor", "#000000")
    config.set("Openbabel", "BackColor", "#ffffff")
    config.set("Openbabel", "BackColorAlpha", "0")
    config.set("Openbabel", "AtomsColor", "True")
    config.set("Openbabel", "AtomsAll", "False")
    config.set("Openbabel", "AtomsEnd", "True")
    config.set("Openbabel", "AtomsNone", "False")
    config.set("Openbabel", "TighBond", "True")

    return config


def config():
    """Function to define a first project config file"""
    config = ConfigParser()

    # Components
    config.add_section("Components")
    config.set("Components", "Components", "[]")
    config.set("Components", "Solids", "[]")

    # Thermodynamics
    config.add_section("Thermo")
    config.set("Thermo", "K", "0")
    config.set("Thermo", "Alfa", "0")
    config.set("Thermo", "Mixing", "0")
    config.set("Thermo", "H", "0")
    config.set("Thermo", "Cp_ideal", "0")
    config.set("Thermo", "MEoS", "False")
    config.set("Thermo", "iapws", "False")
    config.set("Thermo", "GERG", "False")
    config.set("Thermo", "freesteam", "False")
    config.set("Thermo", "coolProp", "False")
    config.set("Thermo", "refprop", "False")

    # Transport
    config.add_section("Transport")
    config.set("Transport", "RhoL", "0")
    config.set("Transport", "Corr_RhoL", "0")
    config.set("Transport", "MuL", "0")
    config.set("Transport", "Corr_MuL", "0")
    config.set("Transport", "Corr_MuG", "0")
    config.set("Transport", "MuG", "0")
    config.set("Transport", "Tension", "0")
    config.set("Transport", "ThCondL", "0")
    config.set("Transport", "Corr_ThCondL", "0")
    config.set("Transport", "ThCondG", "0")
    config.set("Transport", "Pv", "0")
    config.set("Transport", "f_acent", "0")
    config.set("Transport", "RhoLMix", "0")
    config.set("Transport", "Corr_RhoLMix", "0")
    config.set("Transport", "MuLMix", "0")
    config.set("Transport", "MuGMix", "0")
    config.set("Transport", "Corr_MuGMix", "0")
    config.set("Transport", "ThCondLMix", "0")
    config.set("Transport", "ThCondGMix", "0")
    config.set("Transport", "Corr_ThCondGMix", "0")

    # Units
    config.add_section("Units")
    config.set("Units", "System", "0")
    for magnitud in magnitudes[:-1]:
        config.set("Units", magnitud, "0")

    # Resolution
    config.add_section("PFD")
    config.set("PFD", "x", "600")
    config.set("PFD", "y", "480")

    return config


def getrates(archivo):
    """Procedure to update change rates"""
    rates = {}
    date = 0

    url = "https://finance.yahoo.com/webservice/v1/symbols/allcurrencies/" \
        "quote?format=json"
    fh = urllib.request.urlopen(url)
    data = json.loads(fh.read().decode("utf-8"))

    for change in data["list"]["resources"]:
        tas = change["resource"]["fields"]
        if "USD" not in tas["name"]:
            continue
        name = tas["name"]
        if name != "USD":
            name = name.split("/")[1]
        rates[name.lower()] = float(tas["price"])

        # Set the last date number
        if int(tas["ts"]) > date:
            date = int(tas["ts"])

    rates["date"] = date
    json.dump(rates, open(archivo, "w"), indent=4)


def createDatabase(name):
    """Create empty database"""
    conn = sqlite3.connect(name)
    curs = conn.cursor()
    curs.execute("""
                 CREATE TABLE compuestos (
                 id  INTEGER PRIMARY KEY,
                 formula TEXT,
                 name TEXT,
                 M FLOAT,
                 tc          FLOAT,
                 pc          FLOAT,
                 vc          FLOAT,
                 API         FLOAT,
                 Cp_ideal_A    FLOAT,
                 Cp_ideal_B    FLOAT,
                 Cp_ideal_C    FLOAT,
                 Cp_ideal_D    FLOAT,
                 Cp_ideal_E    FLOAT,
                 Cp_ideal_F    FLOAT,
                 antoine_A   FLOAT,
                 antoine_B   FLOAT,
                 antoine_C   FLOAT,
                 henry_A     FLOAT,
                 henry_B     FLOAT,
                 henry_C     FLOAT,
                 henry_D     FLOAT,
                 visco_A     FLOAT,
                 visco_B     FLOAT,
                 tension_A       FLOAT,
                 tension_B       FLOAT,
                 rhoS_DIPPR_EQ   INTEGER,
                 rhoS_DIPPR_A   FLOAT,
                 rhoS_DIPPR_B   FLOAT,
                 rhoS_DIPPR_C   FLOAT,
                 rhoS_DIPPR_D   FLOAT,
                 rhoS_DIPPR_E  FLOAT,
                 rhoS_DIPPR_tmin   FLOAT,
                 rhoS_DIPPR_tmax   FLOAT,
                 rhoL_DIPPR_EQ   INTEGER,
                 rhoL_DIPPR_A   FLOAT,
                 rhoL_DIPPR_B   FLOAT,
                 rhoL_DIPPR_C   FLOAT,
                 rhoL_DIPPR_D   FLOAT,
                 rhoL_DIPPR_E  FLOAT,
                 rhoL_DIPPR_tmin   FLOAT,
                 rhoL_DIPPR_tmax   FLOAT,
                 Pv_DIPPR_EQ   INTEGER,
                 Pv_DIPPR_A   FLOAT,
                 Pv_DIPPR_B   FLOAT,
                 Pv_DIPPR_C   FLOAT,
                 Pv_DIPPR_D   FLOAT,
                 Pv_DIPPR_E  FLOAT,
                 Pv_DIPPR_tmin   FLOAT,
                 Pv_DIPPR_tmax   FLOAT,
                 Hv_DIPPR_EQ   INTEGER,
                 Hv_DIPPR_A   FLOAT,
                 Hv_DIPPR_B   FLOAT,
                 Hv_DIPPR_C   FLOAT,
                 Hv_DIPPR_D   FLOAT,
                 Hv_DIPPR_E  FLOAT,
                 Hv_DIPPR_tmin   FLOAT,
                 Hv_DIPPR_tmax   FLOAT,
                 CpS_DIPPR_EQ   INTEGER,
                 CpS_DIPPR_A   FLOAT,
                 CpS_DIPPR_B   FLOAT,
                 CpS_DIPPR_C   FLOAT,
                 CpS_DIPPR_D   FLOAT,
                 CpS_DIPPR_E  FLOAT,
                 CpS_DIPPR_tmin   FLOAT,
                 CpS_DIPPR_tmax   FLOAT,
                 CpL_DIPPR_EQ   INTEGER,
                 CpL_DIPPR_A   FLOAT,
                 CpL_DIPPR_B   FLOAT,
                 CpL_DIPPR_C   FLOAT,
                 CpL_DIPPR_D   FLOAT,
                 CpL_DIPPR_E  FLOAT,
                 CpL_DIPPR_tmin   FLOAT,
                 CpL_DIPPR_tmax   FLOAT,
                 CpG_DIPPR_EQ   INTEGER,
                 CpG_DIPPR_A   FLOAT,
                 CpG_DIPPR_B   FLOAT,
                 CpG_DIPPR_C   FLOAT,
                 CpG_DIPPR_D   FLOAT,
                 CpG_DIPPR_E  FLOAT,
                 CpG_DIPPR_tmin   FLOAT,
                 CpG_DIPPR_tmax   FLOAT,
                 muL_DIPPR_EQ   INTEGER,
                 muL_DIPPR_A   FLOAT,
                 muL_DIPPR_B   FLOAT,
                 muL_DIPPR_C   FLOAT,
                 muL_DIPPR_D   FLOAT,
                 muL_DIPPR_E  FLOAT,
                 muL_DIPPR_tmin   FLOAT,
                 muL_DIPPR_tmax   FLOAT,
                 muG_DIPPR_EQ   INTEGER,
                 muG_DIPPR_A   FLOAT,
                 muG_DIPPR_B   FLOAT,
                 muG_DIPPR_C   FLOAT,
                 muG_DIPPR_D   FLOAT,
                 muG_DIPPR_E  FLOAT,
                 muG_DIPPR_tmin   FLOAT,
                 muG_DIPPR_tmax   FLOAT,
                 ThcondL_DIPPR_EQ   INTEGER,
                 ThcondL_DIPPR_A   FLOAT,
                 ThcondL_DIPPR_B   FLOAT,
                 ThcondL_DIPPR_C   FLOAT,
                 ThcondL_DIPPR_D   FLOAT,
                 ThcondL_DIPPR_E  FLOAT,
                 ThcondL_DIPPR_tmin   FLOAT,
                 ThcondL_DIPPR_tmax   FLOAT,
                 ThcondG_DIPPR_EQ   INTEGER,
                 ThcondG_DIPPR_A   FLOAT,
                 ThcondG_DIPPR_B   FLOAT,
                 ThcondG_DIPPR_C   FLOAT,
                 ThcondG_DIPPR_D   FLOAT,
                 ThcondG_DIPPR_E  FLOAT,
                 ThcondG_DIPPR_tmin   FLOAT,
                 ThcondG_DIPPR_tmax   FLOAT,
                 tension_DIPPR_EQ   INTEGER,
                 tension_DIPPR_A   FLOAT,
                 tension_DIPPR_B   FLOAT,
                 tension_DIPPR_C   FLOAT,
                 tension_DIPPR_D   FLOAT,
                 tension_DIPPR_E  FLOAT,
                 tension_DIPPR_tmin   FLOAT,
                 tension_DIPPR_tmax   FLOAT,
                 dipole FLOAT,
                 V_liq FLOAT,
                 Rackett   FLOAT,
                 SG FLOAT,
                 f_acent    FLOAT,
                 SolubilityParameter   FLOAT,
                 watson      FLOAT,
                 MSRK_A    FLOAT,
                 MSRK_B    FLOAT,
                 Stiehl  FLOAT,
                 Tb FLOAT,
                 Tf  FLOAT,
                 CAS  TEXT,
                 alternateFormula TEXT,
                 UNIFAC  TEXT,
                 Dm FLOAT,
                 Eps_k   FLOAT,
                 UNIQUAC_area    FLOAT,
                 UNIQUAC_volumen FLOAT,
                 f_acent_MSRK FLOAT,
                 Hf FLOAT,
                 Gf FLOAT,
                 volumen_wilson  FLOAT,
                 NetHeating FLOAT,
                 GrossHeating FLOAT,
                 Synonyms TEXT,
                 volumen_caracteristico FLOAT,
                 calor_formacion_solido  FLOAT,
                 energia_libre_solido    FLOAT,
                 PolarParameter FLOAT,
                 smile   TEXT,
                 antoine_to FLOAT,
                 antoine_n FLOAT,
                 antoine_E FLOAT,
                 antoine_F FLOAT,
                 wagner_a FLOAT,
                 wagner_b FLOAT,
                 wagner_c FLOAT,
                 wagner_d FLOAT)
                 """)
    conn.commit()
    conn.close()
