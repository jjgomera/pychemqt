#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Library to define functions necessary to run at first pychemqt run, so don't
# let import other pychemqt library to avoid import error when no config
# files availables.
#   -createDatabase: Create empty database
###############################################################################


from configparser import ConfigParser
import datetime
import http.client
import json
import logging
import os
import sqlite3
import sys

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
           'ElectricPrecipitator', 'Dryer', 'Scrubber', 'Neumatic',
           'Spreadsheet', 'Reactor', 'Grinder']


def which(program):
    """Function to detect program availability in system and return path"""

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath = os.path.dirname(program)
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
    for programa in ["qalculate", "gnome-calculator", "gcalctool", "mate-calc",
                     "kcalc", "galculator", "deepin-calculator"]:
        ejecutable = which(programa)
        if ejecutable:
            calculator = ejecutable
            break


editor = ""
if sys.platform == "win32":
    editor = "notepad.exe"
else:
    for programa in ["gedit", "leafpad", "l3afpad", "featherpad", "mousepad",
                     "geany", "kate", "kwrite", "vim", "vi", "emacs", "nano",
                     "pico"]:
        ejecutable = which(programa)
        if ejecutable:
            editor = ejecutable
            break

shell = ""
if sys.platform == "win32":
    shell = ""
else:
    shell = which("xterm")
# TODO: For now only support xterm
# for program in ["xterm", "gnome-terminal", "kterminal", "lxterminal",
                # "xfce4-terminal", "terminator"]:
    # exe = which(program)
    # if exe:
        # shell = exe
        # break


def Preferences():
    """Function to define a first preferences file"""
    conf = ConfigParser()

    # General
    conf.add_section("General")
    conf.set("General", "Color_Resaltado", "#ffff00")
    conf.set("General", "Color_ReadOnly", "#eaeaea")
    conf.set("General", "Recent_Files", "10")
    conf.set("General", "Load_Last_Project", "True")
    conf.set("General", "Tray", "False")

    # PFD
    conf.add_section("PFD")
    conf.set("PFD", "x", "800")
    conf.set("PFD", "y", "600")
    conf.set("PFD", "brush", "1")
    conf.set("PFD", "brushColor", "#aaaaaa")
    conf.set("PFD", "Color_Entrada", "#c80000")
    conf.set("PFD", "Color_Salida", "#0000c8")
    conf.set("PFD", "Color_Stream", "#000000")
    conf.set("PFD", "Width", "1.0")
    conf.set("PFD", "Union", "0")
    conf.set("PFD", "Miter_limit", "2.0")
    conf.set("PFD", "Punta", "0")
    conf.set("PFD", "Guion", "0")
    conf.set("PFD", "Dash_offset", "0.0")
    conf.set("PFD", "Move_Factor", "5")

    # TooltipEntity
    conf.add_section("TooltipEntity")
    conf.set("TooltipEntity", "Corriente", "0,1")
    for equipo in equipos:
        conf.set("TooltipEntity", equipo, "0,1")

    # Tooltip
    conf.add_section("Tooltip")
    conf.set("Tooltip", "Show", "True")
    conf.set("Tooltip", "SI", "False")
    conf.set("Tooltip", "CGS", "False")
    conf.set("Tooltip", "AltSI", "False")
    conf.set("Tooltip", "English", "False")
    conf.set("Tooltip", "Metric", "False")
    for magnitud in magnitudes[:-1]:
        conf.set("Tooltip", magnitud, "0,1")

    # NumericFactor
    conf.add_section("NumericFormat")
    for magnitud in magnitudes:
        kwarg = {'total': 0, 'signo': False, 'decimales': 4, 'fmt': 0}
        conf.set("NumericFormat", magnitud, str(kwarg))

    # Petro
    conf.add_section("petro")
    conf.set("petro", "M", "0")
    conf.set("petro", "critical", "0")
    conf.set("petro", "vc", "0")
    conf.set("petro", "f_acent", "0")
    conf.set("petro", "Tb", "0")
    conf.set("petro", "SG", "0")
    conf.set("petro", "n", "0")
    conf.set("petro", "Zc", "0")
    conf.set("petro", "PNA", "0")
    conf.set("petro", "H", "0")
    conf.set("petro", "curve", "0")

    # Applications
    conf.add_section("Applications")
    conf.set("Applications", "Calculator", calculator)
    conf.set("Applications", "TextViewer", editor)
    conf.set("Applications", "PDF", "False")
    conf.set("Applications", "PDFExternal", "")
    conf.set("Applications", "Shell", shell)
    conf.set("Applications", "ipython", "False")
    conf.set("Applications", "maximized", "False")
    conf.set("Applications", "foregroundColor", "#ffffff")
    conf.set("Applications", "backgroundColor", "#000000")
    conf.set("Applications", "elementalColorby", "serie")
    conf.set("Applications", "elementalDefinition", "10")
    conf.set("Applications", "elementalLog", "False")

    # Plotting
    conf.add_section("Plot")
    conf.set("Plot", "style", "0")
    conf.set("Plot", "customize", "False")
    conf.set("Plot", "axes.axisbelow", "line")
    conf.set("Plot", "axes.formatter.limits", "-5, 5")
    conf.set("Plot", "axes.edgecolor", "#000000")
    conf.set("Plot", "axes.facecolor", "#ffffff")
    conf.set("Plot", "axes.formatter.min_exponent", "0")
    conf.set("Plot", "axes.formatter.offset_threshold", "4")
    conf.set("Plot", "axes.formatter.use_locale", "False")
    conf.set("Plot", "axes.formatter.use_mathtext", "False")
    conf.set("Plot", "axes.formatter.useoffset", "True")
    conf.set("Plot", "axes.grid", "False")
    conf.set("Plot", "axes.grid.axis", "both")
    conf.set("Plot", "axes.grid.which", "major")
    conf.set("Plot", "axes.labelcolor", "#000000")
    conf.set("Plot", "axes.labelpad", "4")
    conf.set("Plot", "axes.labelsize", "medium")
    conf.set("Plot", "axes.labelweight", "normal")
    conf.set("Plot", "axes.linewidth", "0.8")
    conf.set("Plot", "axes.spines.bottom", "True")
    conf.set("Plot", "axes.spines.left", "True")
    conf.set("Plot", "axes.spines.right", "True")
    conf.set("Plot", "axes.spines.top", "True")
    conf.set("Plot", "axes.titlecolor", "#000000")
    conf.set("Plot", "axes.titlelocation", "center")
    conf.set("Plot", "axes.titlepad", "6")
    conf.set("Plot", "axes.titlesize", "large")
    conf.set("Plot", "axes.titleweight", "normal")
    conf.set("Plot", "axes.titley", "1.0")
    conf.set("Plot", "axes.unicode_minus", "True")
    conf.set("Plot", "axes.xmargin", "0.05")
    conf.set("Plot", "axes.ymargin", "0.05")
    conf.set("Plot", "axes.zmargin", "0.05")
    conf.set("Plot", "axes3d.grid", "True")
    conf.set("Plot", "xaxis.labellocation", "center")
    conf.set("Plot", "yaxis.labellocation", "center")
    conf.set("Plot", "figure.autolayout", "False")
    conf.set("Plot", "figure.constrained_layout.h_pad", "0.04")
    conf.set("Plot", "figure.constrained_layout.hspace", "0.02")
    conf.set("Plot", "figure.constrained_layout.w_pad", "0.04")
    conf.set("Plot", "figure.constrained_layout.wspace", "0.02")
    conf.set("Plot", "figure.dpi", "100")
    conf.set("Plot", "figure.edgecolor", "#ffffff")
    conf.set("Plot", "figure.facecolor", "#ffffff")
    conf.set("Plot", "figure.frameon", "True")
    conf.set("Plot", "figure.labelsize", "large")
    conf.set("Plot", "figure.labelweight", "normal")
    conf.set("Plot", "figure.subplot.bottom", "0.11")
    conf.set("Plot", "figure.subplot.hspace", "0.2")
    conf.set("Plot", "figure.subplot.left", "0.13")
    conf.set("Plot", "figure.subplot.right", "0.9")
    conf.set("Plot", "figure.subplot.top", "0.88")
    conf.set("Plot", "figure.subplot.wspace", "0.2")
    conf.set("Plot", "figure.titlesize", "large")
    conf.set("Plot", "figure.titleweight", "normal")
    conf.set("Plot", "font.family", "sans-serif")
    conf.set("Plot", "font.size", "10.0")
    conf.set("Plot", "font.stretch", "normal")
    conf.set("Plot", "font.style", "normal")
    conf.set("Plot", "font.variant", "normal")
    conf.set("Plot", "font.weight", "normal")
    conf.set("Plot", "grid.alpha", "1.0")
    conf.set("Plot", "grid.color", "#b0b0b0")
    conf.set("Plot", "grid.linestyle", "-")
    conf.set("Plot", "grid.linewidth", "0.8")
    conf.set("Plot", "hatch.color", "#000000")
    conf.set("Plot", "hatch.linewidth", "1.0")
    conf.set("Plot", "legend.borderaxespad", "0.5")
    conf.set("Plot", "legend.borderpad", "0.4")
    conf.set("Plot", "legend.columnspacing", "2.0")
    conf.set("Plot", "legend.edgecolor", "#cccccc")
    conf.set("Plot", "legend.facecolor", "#ffffff")
    conf.set("Plot", "legend.fancybox", "True")
    conf.set("Plot", "legend.fontsize", "medium")
    conf.set("Plot", "legend.framealpha", "0.8")
    conf.set("Plot", "legend.frameon", "True")
    conf.set("Plot", "legend.handleheight", "0.7")
    conf.set("Plot", "legend.handlelength", "2.0")
    conf.set("Plot", "legend.handletextpad", "0.8")
    conf.set("Plot", "legend.labelcolor", "#000000")
    conf.set("Plot", "legend.labelspacing", "0.5")
    conf.set("Plot", "legend.loc", "best")
    conf.set("Plot", "legend.markerscale", "1.0")
    conf.set("Plot", "legend.numpoints", "1")
    conf.set("Plot", "legend.scatterpoints", "1")
    conf.set("Plot", "legend.shadow", "False")
    conf.set("Plot", "legend.title_fontsize", "xx-small")
    conf.set("Plot", "lines.antialiased", "True")
    conf.set("Plot", "lines.color", "#1f77b4")
    conf.set("Plot", "lines.dash_capstyle", "butt")
    conf.set("Plot", "lines.dash_joinstyle", "round")
    conf.set("Plot", "lines.linestyle", "-")
    conf.set("Plot", "lines.linewidth", "1.5")
    conf.set("Plot", "lines.marker", "None")
    conf.set("Plot", "lines.markeredgecolor", "#000000")
    conf.set("Plot", "lines.markeredgewidth", "1.0")
    conf.set("Plot", "lines.markerfacecolor", "#2ca02c")
    conf.set("Plot", "lines.markersize", "6.0")
    conf.set("Plot", "lines.scale_dashes", "True")
    conf.set("Plot", "lines.solid_capstyle", "projecting")
    conf.set("Plot", "lines.solid_joinstyle", "round")
    conf.set("Plot", "patch.antialiased", "True")
    conf.set("Plot", "patch.edgecolor", "#000000")
    conf.set("Plot", "patch.facecolor", "#1f77b4")
    conf.set("Plot", "patch.force_edgecolor", "False")
    conf.set("Plot", "patch.linewidth", "1.0")
    conf.set("Plot", "savefig.bbox", "standard")
    conf.set("Plot", "savefig.dpi", "100")
    conf.set("Plot", "savefig.edgecolor", "#ffffff")
    conf.set("Plot", "savefig.facecolor", "#ffffff")
    conf.set("Plot", "savefig.format", "png")
    conf.set("Plot", "savefig.pad_inches", "0.1")
    conf.set("Plot", "savefig.transparent", "False")
    conf.set("Plot", "xtick.alignment", "center")
    conf.set("Plot", "xtick.bottom", "True")
    conf.set("Plot", "xtick.color", "#000000")
    conf.set("Plot", "xtick.direction", "out")
    conf.set("Plot", "xtick.labelbottom", "True")
    conf.set("Plot", "xtick.labelcolor", "#000000")
    conf.set("Plot", "xtick.labelsize", "medium")
    conf.set("Plot", "xtick.labeltop", "False")
    conf.set("Plot", "xtick.major.bottom", "True")
    conf.set("Plot", "xtick.major.pad", "3.5")
    conf.set("Plot", "xtick.major.size", "3.5")
    conf.set("Plot", "xtick.major.top", "True")
    conf.set("Plot", "xtick.major.width", "0.8")
    conf.set("Plot", "xtick.minor.bottom", "True")
    conf.set("Plot", "xtick.minor.pad", "3.4")
    conf.set("Plot", "xtick.minor.size", "2.0")
    conf.set("Plot", "xtick.minor.top", "True")
    conf.set("Plot", "xtick.minor.visible", "False")
    conf.set("Plot", "xtick.minor.width", "0.6")
    conf.set("Plot", "xtick.top", "False")
    conf.set("Plot", "ytick.alignment", "center_baseline")
    conf.set("Plot", "ytick.color", "#000000")
    conf.set("Plot", "ytick.direction", "out")
    conf.set("Plot", "ytick.labelcolor", "#000000")
    conf.set("Plot", "ytick.labelleft", "True")
    conf.set("Plot", "ytick.labelright", "False")
    conf.set("Plot", "ytick.labelsize", "medium")
    conf.set("Plot", "ytick.left", "True")
    conf.set("Plot", "ytick.major.left", "True")
    conf.set("Plot", "ytick.major.pad", "3.5")
    conf.set("Plot", "ytick.major.right", "True")
    conf.set("Plot", "ytick.major.size", "3.5")
    conf.set("Plot", "ytick.major.width", "0.8")
    conf.set("Plot", "ytick.minor.left", "True")
    conf.set("Plot", "ytick.minor.pad", "3.4")
    conf.set("Plot", "ytick.minor.right", "True")
    conf.set("Plot", "ytick.minor.size", "2.0")
    conf.set("Plot", "ytick.minor.visible", "False")
    conf.set("Plot", "ytick.minor.width", "0.6")
    conf.set("Plot", "ytick.right", "False")

    # mEoS
    conf.add_section("MEOS")
    conf.set("MEOS", "coolprop", "False")
    conf.set("MEOS", "refprop", "False")
    conf.set("MEOS", "saturation"+"Color", "#000000")
    conf.set("MEOS", "saturation"+"alpha", "255")
    conf.set("MEOS", "saturation"+"lineWidth", "1.0")
    conf.set("MEOS", "saturation"+"lineStyle", "-")
    conf.set("MEOS", "saturation"+"marker", "None")
    conf.set("MEOS", "saturation"+"markersize", "3")
    conf.set("MEOS", "saturation"+"markerfacecolor", "#ff0000")
    conf.set("MEOS", "saturation"+"markeredgewidth", "1")
    conf.set("MEOS", "saturation"+"markeredgecolor", "#000000")
    conf.set("MEOS", "grid", "False")
    conf.set("MEOS", "definition", "1")
    lineas = ["Isotherm", "Isobar", "Isoenthalpic", "Isoentropic", "Isochor",
              "Isoquality"]
    for linea in lineas:
        conf.set("MEOS", linea+"Start", "0")
        conf.set("MEOS", linea+"End", "0")
        conf.set("MEOS", linea+"Step", "0")
        conf.set("MEOS", linea+"Custom", "True")
        conf.set("MEOS", linea+"List", "")
        if linea != "Isoquality":
            conf.set("MEOS", linea+"Critic", "True")
        conf.set("MEOS", linea+"Color", "#000000")
        conf.set("MEOS", linea+"alpha", "255")
        conf.set("MEOS", linea+"lineWidth", "0.5")
        conf.set("MEOS", linea+"lineStyle", "-")
        conf.set("MEOS", linea+"marker", "None")
        conf.set("MEOS", linea+"markersize", "3")
        conf.set("MEOS", linea+"markerfacecolor", "#ff0000")
        conf.set("MEOS", linea+"markeredgewidth", "1")
        conf.set("MEOS", linea+"markeredgecolor", "#000000")

        conf.set("MEOS", linea+"Label", "False")
        conf.set("MEOS", linea+"Variable", "False")
        conf.set("MEOS", linea+"Units", "False")
        conf.set("MEOS", linea+"Position", "50")

    conf.set("MEOS", "3Dmesh", "False")
    conf.set("MEOS", "3Dtype", "0")
    conf.set("MEOS", "3Dcolormap", "viridis")
    conf.set("MEOS", "3Dalphasurface", "150")
    conf.set("MEOS", "3Dcolor", "#000000")
    conf.set("MEOS", "3Dalpha", "150")
    conf.set("MEOS", "3Dlinewidth", "0.5")
    conf.set("MEOS", "3Dlinestyle", "-")

    # Psychr
    conf.add_section("Psychr")
    conf.set("Psychr", "chart", "True")
    conf.set("Psychr", "virial", "False")
    conf.set("Psychr", "coolprop", "False")
    conf.set("Psychr", "refprop", "False")

    conf.set("Psychr", "saturation"+"Color", "#000000")
    conf.set("Psychr", "saturation"+"alpha", "255")
    conf.set("Psychr", "saturation"+"lineWidth", "0.5")
    conf.set("Psychr", "saturation"+"lineStyle", "-")
    conf.set("Psychr", "saturation"+"marker", "None")
    conf.set("Psychr", "saturation"+"markersize", "3")
    conf.set("Psychr", "saturation"+"markerfacecolor", "#ff0000")
    conf.set("Psychr", "saturation"+"markeredgewidth", "1")
    conf.set("Psychr", "saturation"+"markeredgecolor", "#000000")

    conf.set("Psychr", "crux"+"Color", "#0000ff")
    conf.set("Psychr", "crux"+"alpha", "255")
    conf.set("Psychr", "crux"+"lineWidth", "0.5")
    conf.set("Psychr", "crux"+"lineStyle", "-")
    conf.set("Psychr", "crux"+"marker", "None")
    conf.set("Psychr", "crux"+"markersize", "3")
    conf.set("Psychr", "crux"+"markerfacecolor", "#ff0000")
    conf.set("Psychr", "crux"+"markeredgewidth", "1")
    conf.set("Psychr", "crux"+"markeredgecolor", "#000000")

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
        {"start": 0.7, "end": 1.3, "step": 0.01, "color": "#00aa00",
         "linewidth": 0.8, "linestyle": ":", "label": "False",
         "units": "False", "position": 90}]
    for linea, value in zip(lineas, values):
        conf.set("Psychr", linea+"Start", str(value["start"]))
        conf.set("Psychr", linea+"End", str(value["end"]))
        conf.set("Psychr", linea+"Step", str(value["step"]))
        conf.set("Psychr", linea+"Custom", "False")
        conf.set("Psychr", linea+"List", "")
        conf.set("Psychr", linea+"Color", str(value["color"]))
        conf.set("Psychr", linea+"alpha", "255")
        conf.set("Psychr", linea+"lineWidth", str(value["linewidth"]))
        conf.set("Psychr", linea+"lineStyle", str(value["linestyle"]))
        conf.set("Psychr", linea+"marker", "None")
        conf.set("Psychr", linea+"markersize", "3")
        conf.set("Psychr", linea+"markerfacecolor", "#ff0000")
        conf.set("Psychr", linea+"markeredgewidth", "1")
        conf.set("Psychr", linea+"markeredgecolor", "#000000")
        conf.set("Psychr", linea+"Label", str(value["label"]))
        conf.set("Psychr", linea+"Units", str(value["units"]))
        conf.set("Psychr", linea+"Position", str(value["position"]))
        conf.set("Psychr", linea+"variable", str(False))

    # Moody
    conf.add_section("Moody")
    conf.set("Moody", "fanning", "False")
    conf.set("Moody", "method", "0")
    conf.set("Moody", "ed", "0, 1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, "
             "4e-4, 6e-4, 8e-4, .001, .0015, .002, .003, .004, .006, .008, "
             ".01, .0125, .015, .0175, .02, .025, .03, .035, .04, .045, "
             ".05, .06, .07")

    conf.set("Moody", "line"+"Color", "#000000")
    conf.set("Moody", "line"+"alpha", "255")
    conf.set("Moody", "line"+"lineWidth", "0.5")
    conf.set("Moody", "line"+"lineStyle", "-")
    conf.set("Moody", "line"+"marker", "None")
    conf.set("Moody", "line"+"markersize", "3")
    conf.set("Moody", "line"+"markerfacecolor", "#ff0000")
    conf.set("Moody", "line"+"markeredgewidth", "1")
    conf.set("Moody", "line"+"markeredgecolor", "#000000")

    conf.set("Moody", "crux"+"Color", "#0000ff")
    conf.set("Moody", "crux"+"alpha", "255")
    conf.set("Moody", "crux"+"lineWidth", "0.5")
    conf.set("Moody", "crux"+"lineStyle", "-")
    conf.set("Moody", "crux"+"marker", "None")
    conf.set("Moody", "crux"+"markersize", "3")
    conf.set("Moody", "crux"+"markerfacecolor", "#ff0000")
    conf.set("Moody", "crux"+"markeredgewidth", "1")
    conf.set("Moody", "crux"+"markeredgecolor", "#000000")

    conf.set("Moody", "grid", "True")
    conf.set("Moody", "grid"+"Color", "#000000")
    conf.set("Moody", "grid"+"alpha", "255")
    conf.set("Moody", "grid"+"lineWidth", "0.5")
    conf.set("Moody", "grid"+"lineStyle", ":")
    conf.set("Moody", "grid"+"marker", "None")
    conf.set("Moody", "grid"+"markersize", "3")
    conf.set("Moody", "grid"+"markerfacecolor", "#ff0000")
    conf.set("Moody", "grid"+"markeredgewidth", "1")
    conf.set("Moody", "grid"+"markeredgecolor", "#000000")
    conf.set("Moody", "grid"+"which", "both")
    conf.set("Moody", "grid"+"axis", "both")

    # Standing-Katz
    conf.add_section("Standing_Katz")
    conf.set("Standing_Katz", "method", "0")
    conf.set("Standing_Katz", "Tr", "1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35,"
             "1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.")

    conf.set("Standing_Katz", "line"+"Color", "#000000")
    conf.set("Standing_Katz", "line"+"alpha", "255")
    conf.set("Standing_Katz", "line"+"lineWidth", "0.5")
    conf.set("Standing_Katz", "line"+"lineStyle", "-")
    conf.set("Standing_Katz", "line"+"marker", "None")
    conf.set("Standing_Katz", "line"+"markersize", "3")
    conf.set("Standing_Katz", "line"+"markerfacecolor", "#ff0000")
    conf.set("Standing_Katz", "line"+"markeredgewidth", "1")
    conf.set("Standing_Katz", "line"+"markeredgecolor", "#000000")

    conf.set("Standing_Katz", "crux"+"Color", "#0000ff")
    conf.set("Standing_Katz", "crux"+"alpha", "255")
    conf.set("Standing_Katz", "crux"+"lineWidth", "0.5")
    conf.set("Standing_Katz", "crux"+"lineStyle", "-")
    conf.set("Standing_Katz", "crux"+"marker", "None")
    conf.set("Standing_Katz", "crux"+"markersize", "3")
    conf.set("Standing_Katz", "crux"+"markerfacecolor", "#ff0000")
    conf.set("Standing_Katz", "crux"+"markeredgewidth", "1")
    conf.set("Standing_Katz", "crux"+"markeredgecolor", "#000000")

    conf.set("Standing_Katz", "grid", "True")
    conf.set("Standing_Katz", "grid"+"Color", "#000000")
    conf.set("Standing_Katz", "grid"+"alpha", "255")
    conf.set("Standing_Katz", "grid"+"lineWidth", "0.5")
    conf.set("Standing_Katz", "grid"+"lineStyle", ":")
    conf.set("Standing_Katz", "grid"+"marker", "None")
    conf.set("Standing_Katz", "grid"+"markersize", "3")
    conf.set("Standing_Katz", "grid"+"markerfacecolor", "#ff0000")
    conf.set("Standing_Katz", "grid"+"markeredgewidth", "1")
    conf.set("Standing_Katz", "grid"+"markeredgecolor", "#000000")
    conf.set("Standing_Katz", "grid"+"which", "both")
    conf.set("Standing_Katz", "grid"+"axis", "both")

    # drag sphere
    conf.add_section("drag")
    conf.set("drag", "method", "0")

    conf.set("drag", "line"+"Color", "#000000")
    conf.set("drag", "line"+"alpha", "255")
    conf.set("drag", "line"+"lineWidth", "0.5")
    conf.set("drag", "line"+"lineStyle", "-")
    conf.set("drag", "line"+"marker", "None")
    conf.set("drag", "line"+"markersize", "3")
    conf.set("drag", "line"+"markerfacecolor", "#ff0000")
    conf.set("drag", "line"+"markeredgewidth", "1")
    conf.set("drag", "line"+"markeredgecolor", "#000000")

    conf.set("drag", "crux"+"Color", "#0000ff")
    conf.set("drag", "crux"+"alpha", "255")
    conf.set("drag", "crux"+"lineWidth", "0.5")
    conf.set("drag", "crux"+"lineStyle", "-")
    conf.set("drag", "crux"+"marker", "None")
    conf.set("drag", "crux"+"markersize", "3")
    conf.set("drag", "crux"+"markerfacecolor", "#ff0000")
    conf.set("drag", "crux"+"markeredgewidth", "1")
    conf.set("drag", "crux"+"markeredgecolor", "#000000")

    conf.set("drag", "grid", "True")
    conf.set("drag", "grid"+"Color", "#000000")
    conf.set("drag", "grid"+"alpha", "255")
    conf.set("drag", "grid"+"lineWidth", "0.5")
    conf.set("drag", "grid"+"lineStyle", ":")
    conf.set("drag", "grid"+"marker", "None")
    conf.set("drag", "grid"+"markersize", "3")
    conf.set("drag", "grid"+"markerfacecolor", "#ff0000")
    conf.set("drag", "grid"+"markeredgewidth", "1")
    conf.set("drag", "grid"+"markeredgecolor", "#000000")
    conf.set("drag", "grid"+"which", "both")
    conf.set("drag", "grid"+"axis", "both")

    # Openbabel
    conf.add_section("Openbabel")
    conf.set("Openbabel", "BondColor", "#000000")
    conf.set("Openbabel", "BackgroundColor", "#ffffff")
    conf.set("Openbabel", "BackColorTransparent", "True")
    conf.set("Openbabel", "AtomsColor", "True")
    conf.set("Openbabel", "AtomsAll", "False")
    conf.set("Openbabel", "AtomsEnd", "True")
    conf.set("Openbabel", "AtomsNone", "False")
    conf.set("Openbabel", "TighBond", "False")
    conf.set("Openbabel", "AsymetricDouble", "True")

    return conf


def config():
    """Function to define a first project config file"""
    conf = ConfigParser()

    # Components
    conf.add_section("Components")
    conf.set("Components", "Components", "[]")
    conf.set("Components", "Solids", "[]")

    # Thermodynamics
    conf.add_section("Thermo")
    conf.set("Thermo", "K", "0")
    conf.set("Thermo", "Alfa", "0")
    conf.set("Thermo", "Mixing", "0")
    conf.set("Thermo", "H", "0")
    conf.set("Thermo", "Cp_ideal", "0")
    conf.set("Thermo", "MEoS", "False")
    conf.set("Thermo", "iapws", "False")
    conf.set("Thermo", "GERG", "False")
    conf.set("Thermo", "freesteam", "False")
    conf.set("Thermo", "coolProp", "False")
    conf.set("Thermo", "refprop", "False")

    # Transport
    conf.add_section("Transport")
    conf.set("Transport", "RhoL", "0")
    conf.set("Transport", "MuL", "0")
    conf.set("Transport", "MuG", "0")
    conf.set("Transport", "Tension", "0")
    conf.set("Transport", "ThCondL", "0")
    conf.set("Transport", "ThCondG", "0")
    conf.set("Transport", "Pv", "0")
    conf.set("Transport", "f_acent", "0")

    conf.set("Transport", "Corr_RhoL", "0")
    conf.set("Transport", "Corr_MuL", "0")
    conf.set("Transport", "Corr_MuG", "0")
    conf.set("Transport", "Corr_ThCondL", "0")
    conf.set("Transport", "Corr_ThCondG", "0")

    conf.set("Transport", "RhoLMix", "0")
    conf.set("Transport", "MuLMix", "0")
    conf.set("Transport", "MuGMix", "0")
    conf.set("Transport", "ThCondLMix", "0")
    conf.set("Transport", "ThCondGMix", "0")

    conf.set("Transport", "Corr_RhoLMix", "0")
    conf.set("Transport", "Corr_MuGMix", "0")
    conf.set("Transport", "Corr_ThCondGMix", "0")

    conf.set("Transport", "RhoLEoS", "False")

    # Units
    conf.add_section("Units")
    conf.set("Units", "System", "0")
    for magnitud in magnitudes[:-1]:
        conf.set("Units", magnitud, "0")

    # Resolution
    conf.add_section("PFD")
    conf.set("PFD", "x", "600")
    conf.set("PFD", "y", "480")

    return conf


def getrates(filename):
    """Procedure to update change rates"""
    rates = {}

    # Get date from old file to avoid bad use of server, only one use a day
    try:
        archivo = open(filename, "r")
        olddate = datetime.date.fromisoformat(json.load(archivo)["date"])
        date = datetime.date.today()
        if date <= olddate:
            logging.info("Currency, using saved data")
            return
    except FileNotFoundError:
        pass

    conn = http.client.HTTPSConnection("api.currencyscoop.com")
    conn.request("GET", "/v1/latest?api_key=c92991dc6b5f33389fcff081bcede004")

    # Alternate web service
    # key = "cc1eW8osv25urLFAf1HUGxZe1MotMDhS"
    # conn = http.client.HTTPSConnection("api.apilayer.com")
    # conn.request("GET", "/exchangerates_data/latest?base=USD&apikey=" + key)

    res = conn.getresponse().read().decode("utf-8")
    data = json.loads(res)["response"]["rates"]
    date = json.loads(res)["response"]["date"][:10]

    for iso, value in data.items():
        if value:
            rates[iso.lower()] = 1/value

    rates["date"] = date
    json.dump(rates, open(filename, "w"), indent=4)


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
