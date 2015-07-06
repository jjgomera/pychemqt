#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to define functions necessary to run at first pychemqt run, so don't
# let import other pychemqt library to avoid import error when no config
# files availables.
###############################################################################

import sys
from configparser import ConfigParser
import urllib.request
import pickle

# It must be defined previously to avoid to early import of libraries
# See end of lib/unidades.py to know how to get this list, check when new
# magnitude are added
magnitudes = ['Temperature', 'DeltaT', 'Angle', 'Length', 'ParticleDiameter',
              'Thickness', 'PipeDiameter', 'Head', 'Area', 'Volume', 'VolLiq',
              'VolGas', 'Time', 'Frequency', 'Speed', 'Acceleration', 'Mass',
              'Mol', 'SpecificVolume', 'SpecificVolume_square', 'MolarVolume',
              'Density', 'DenLiq', 'DenGas', 'MolarDensity', 'Force',
              'Pressure', 'DeltaP', 'Energy', 'Work', 'Enthalpy',
              'MolarEnthalpy', 'Entropy', 'SpecificHeat', 'SpecificEntropy',
              'MolarSpecificHeat', 'EnergyFlow', 'Power', 'MassFlow',
              'MolarFlow', 'VolFlow', 'QLiq', 'QGas', 'Diffusivity',
              'KViscosity', 'HeatFlux', 'ThermalConductivity', 'UA',
              'HeatTransfCoef', 'Fouling', 'Tension', 'Viscosity',
              'SolubilityParameter', 'PotencialElectric', 'DipoleMoment',
              'CakeResistance', 'PackingDP', 'V2V', 'InvTemperature',
              'InvPressure', 'EnthalpyPressure', 'EnthalpyDensity',
              'TemperaturePressure', 'PressureTemperature', 'PressureDensity',
              'DensityPressure', 'DensityTemperature', 'Currency',
              'Dimensionless']
# See end of equipment.__init__.py to know how to get this list, check when new
# fully functional are added
equipos = ['Divider', 'Valve', 'Mixer', 'Pump', 'Compressor', 'Turbine',
           'Pipe', 'Flash', 'ColumnFUG', 'Heat_Exchanger', 'Shell_Tube',
           'Hairpin', 'Fired_Heater', 'Ciclon', 'GravityChamber', 'Baghouse',
           'ElectricPrecipitator', 'Dryer', 'Scrubber', 'Spreadsheet']


def which(program):
    """Function to detect program availability in systemi and return path"""
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
    config.set("petro", "molecular_weight", "0")
    config.set("petro", "critical", "0")
    config.set("petro", "vc", "0")
    config.set("petro", "f_acent", "0")
    config.set("petro", "t_ebull", "0")
    config.set("petro", "Zc", "0")
    config.set("petro", "PNA", "0")
    config.set("petro", "H", "0")
    config.set("petro", "curva", "0")

    # Applications
    config.add_section("Applications")
    config.set("Applications", "Calculator", calculator)
    config.set("Applications", "TextViewer", editor)
    config.set("Applications", "Shell", shell)
    config.set("Applications", "ipython", False)
    config.set("Applications", "maximized", False)
    config.set("Applications", "foregroundColor", "#ffffff")
    config.set("Applications", "backgroundColor", "#000000")

    # mEoS
    config.add_section("MEOS")
    config.set("MEOS", "coolprop", "False")
    config.set("MEOS", "refprop", "False")
    config.set("MEOS", "saturation"+"Color", "#000000")
    config.set("MEOS", "saturation"+"lineWidth", "1.0")
    config.set("MEOS", "saturation"+"lineStyle", "-")
    config.set("MEOS", "saturation"+"marker", "None")
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
        config.set("MEOS", linea+"lineWidth", "0.5")
        config.set("MEOS", linea+"lineStyle", "-")
        config.set("MEOS", linea+"marker", "None")

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
    config.set("Psychr", "saturation"+"lineWidth", "0.5")
    config.set("Psychr", "saturation"+"lineStyle", "-")
    config.set("Psychr", "saturation"+"marker", "None")
    lineas = ["IsoTdb", "IsoW", "IsoHR", "IsoTwb", "Isochor"]
    values = [
        {"start": 274.0, "end": 330.0, "step": 1.0, "color": "#000000", 
        "linewidth": 0.5, "linestyle": ":", "label": "False", "units": "False", 
        "position": 50}, 
        {"start": 0.0, "end": 0.04, "step": 0.001, "color": "#000000", 
        "linewidth": 0.5, "linestyle": ":", "label": "False", "units": "False", 
        "position": 50}, 
        {"start": 10.0, "end": 100.0, "step": 10.0, "color": "#000000", 
        "linewidth": 0.5, "linestyle": "--", "label": "True", "units": "True", 
        "position": 85}, 
        {"start": 250.0, "end": 320.0, "step": 1.0, "color": "#aa0000", 
        "linewidth": 0.8, "linestyle": ":", "label": "False", "units": "False", 
        "position": 90}, 
        {"start": 0.8, "end": 1.0, "step": 0.01, "color": "#00aa00", 
        "linewidth": 0.8, "linestyle": ":", "label": "False", "units": "False", 
        "position": 90}]
    for linea, value in zip(lineas, values):
        config.set("Psychr", linea+"Start", value["start"])
        config.set("Psychr", linea+"End", value["end"])
        config.set("Psychr", linea+"Step", value["step"])
        config.set("Psychr", linea+"Custom", "False")
        config.set("Psychr", linea+"List", "")
        config.set("Psychr", linea+"Color", value["color"])
        config.set("Psychr", linea+"lineWidth", value["linewidth"])
        config.set("Psychr", linea+"lineStyle", value["linestyle"])
        config.set("Psychr", linea+"marker", "None")
        config.set("Psychr", linea+"Label", value["label"])
        config.set("Psychr", linea+"Units", value["units"])
        config.set("Psychr", linea+"Position", value["position"])
    
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
    config.set("Transport", "MuG", "0")
    config.set("Transport", "Tension", "0")
    config.set("Transport", "ThCondL", "0")
    config.set("Transport", "Corr_ThCondL", "0")
    config.set("Transport", "ThCondG", "0")
    config.set("Transport", "Pv", "0")

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


def getrates(archivo):  # From Python Cookbook
    """Procedure to update change rates"""
    rates = {}
    url = "http://www.bankofcanada.ca/en/markets/csv/exchange_eng.csv"
    fh = urllib.request.urlopen(url)
    for line in fh:
        line = line.rstrip()
        if not line or line.startswith(("#", "Closing ")):
            continue
        fields = line.split(",")
        if line.startswith("Date "):
            date = fields[-1]
        elif line.startswith("U.S. dollar (close)"):
            pass
        else:
            value = float(fields[-1])
            rates[fields[1][1:].lower()] = value
    del rates["iexe0124"]
    del rates["iexe0125"]
    rates["cad"] = 1.
    for rate in rates:
        rates[rate] = rates[rate] / rates["usd"]
    rates["date"] = date
    pickle.dump(rates, open(archivo, "w"))
