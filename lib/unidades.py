#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=missing-class-docstring

# The unidad sublasses are autodocumented

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.



Phisics quantities module with support for unit conversion.

:class:`unidad`: Base class with all functionality

Using the dimension symbol for define the SI base units:

    * Length: L
    * Mass: M
    * Time: T
    * Temperature: Θ
    * Amount of substance: N
    * Electric current: I

The list of supported unit are:

    * :class:`Acceleration`, LT⁻²
    * :class:`Angle`
    * :class:`Area`, L²
    * :class:`CakeResistance`, LM⁻¹
    * :class:`Currency`
    * :class:`Density`, ML⁻³
    * :class:`DensityPressure`, T²L⁻²
    * :class:`DensityTemperature`, ML⁻³Θ⁻¹
    * :class:`Diffusivity`, L²T⁻¹
    * :class:`DipoleMoment`, IL
    * :class:`PotencialElectric`, ML³T⁻³I⁻¹
    * :class:`Energy`, ML²T⁻²
    * :class:`Enthalpy`, L²T⁻²
    * :class:`EnthalpyDensity`, L⁵MT⁻²
    * :class:`EnthalpyPressure`, L³M⁻¹
    * :class:`Pressure`, ML⁻¹T⁻²
    * :class:`Entropy`, ML²T⁻²Θ
    * :class:`Force`, MLT⁻²
    * :class:`Fouling`
    * :class:`Frequency`, T⁻¹
    * :class:`V2V`
    * :class:`HeatFlux`, MT⁻³
    * :class:`HeatTransfCoef`, MT⁻⁴
    * :class:`Length`, L
    * :class:`Mass`, M
    * :class:`MassFlow`, MT⁻¹
    * :class:`Mol`, N
    * :class:`MolarDensity`, NL⁻³
    * :class:`MolarEnthalpy`
    * :class:`MolarFlow`, NT⁻¹
    * :class:`MolarSpecificHeat`
    * :class:`MolarVolume`, L⁻³N
    * :class:`PackingDP`, ML⁻²T⁻²
    * :class:`Power`, ML²T⁻³
    * :class:`Pressure`, ML⁻¹T⁻²
    * :class:`DeltaP`, ML⁻¹T⁻²
    * :class:`InvPressure`, LT²M⁻¹
    * :class:`PressureTemperature`, ML⁻¹T⁻²Θ⁻¹
    * :class:`PressureDensity`, L²T⁻²
    * :class:`SolubilityParameter`, M⁰⁵L⁻⁰⁵T⁻¹
    * :class:`SpecificHeat`, L²T⁻²Θ
    * :class:`SpecificVolume`, L³M⁻¹
    * :class:`Speed`, LT⁻¹
    * :class:`Tension`, MT⁻²
    * :class:`Temperature`, Θ
    * :class:`DeltaT`, Θ
    * :class:`InvTemperature`, Θ⁻¹
    * :class:`TemperaturePressure`, ΘLT²M⁻¹
    * :class:`ThermalConductivity`, ML¹T⁻³Θ⁻¹
    * :class:`SpecificVolume_square`, L⁶M⁻²
    * :class:`Time`, T
    * :class:`UA`, ML²T⁻³Θ⁻¹
    * :class:`Viscosity`, ML⁻¹T⁻¹
    * :class:`Volume`, L³
    * :class:`VolFlow`, L³T⁻¹
    * :class:`Dimensionless`: Null unit

'''


from configparser import ConfigParser
import json
import logging
import os

from PyQt5.QtWidgets import QApplication
from PyQt5 import QtCore
import scipy.constants as k

from lib.config import conf_dir, getMainWindowConfig
from lib.utilities import representacion
from tools.firstrun import getrates

# Defining conversion factor not available in scipy
k.tonUS = 2000 * k.lb
k.tonUK = 2240 * k.lb
k.TNT = 4.184e12
k.BarrilOil = 5.8e6 * k.Btu
k.TmOil = 4.1868e10
k.TmCoal = 2.93e10
k.slug = k.lb * k.g / k.foot
k.debye = 1. / k.c * k.zepto
k.pdl = k.foot * k.lb
k.Rankine = 1 / 1.8  # only for differences
k.Reaumur = 1 / 0.8  # only for differences
k.milla = 1609.344
k.milla_nau = 1852
k.acre = 4840 * k.yard**2
k.qtUSliq = 0.25 * k.gallon
k.qtUSdry = 67.2 * k.inch**3
k.qtUK = 0.25 * k.gallon_imp
k.cwtUS = 112 * k.lb
k.cwtUK = 100 * k.lb
k.CV = k.kgf * 75
k.ozf = k.oz * k.g
k.TonfUK = k.g * k.tonUK
k.TonfUS = k.g * k.tonUS
k.statV = 300


def C2K(C):
    """Convert Celcius to Kelvin"""
    return C + 273.15


def K2C(K):
    """Convert Kelvin to Celcius"""
    return K - 273.15


def K2R(K):
    """Convert Kelvin to Rankine"""
    return K * 1.8


def R2K(K):
    """Convert Rankine to Kelvin"""
    return K / 1.8


def F2K(F):
    """Convert Fahrenheit to Kelvin"""
    return ((F - 32) / 1.8) + 273.15


def K2F(K):
    """Convert Kelvin to Fahrenheit"""
    return 1.8 * (K - 273.15) + 32


def Re2K(Re):
    """Convert Reaumur to Kelvin"""
    return (Re * 1.25) + 273.15


def K2Re(K):
    """Convert Kelvin to Reaumur"""
    return (K - 273.15) / 1.25


class unidad(float):
    """
    Generic class to model units.

    Each child class must define the following parameters:
        * __title__: Title or name of class
        * rates: Dict with conversion rates
        * __text__: List with units title
        * __units__: List with units properties names
        * __tooltip__: List with help string for units
        * _magnitudes: Opcional to units with several magnituds. Each magnitud
        is a tuple with format (Name, title)
        * __units_set__: Dict with standart unit for units system (*altsi*,
        *si*, *metric*, *cgs*, *english*)
    """
    __title__ = ""
    rates = {}
    __text__ = []
    __units__ = []
    __tooltip__ = []
    _magnitudes = []
    __units_set__ = []

    def __init__(self, data, unit="", magnitud=""):
        """Constructor

        Parameters
        ----------
        data : float
            Value of initial entity
        unit : str
            String with unit of data value input
        magnitud : str, optional
            Name of magnitud (i.e. PipeDiameter or Head for length unit)

        Notes
        -----
        Non proportional magnitudes (Temperature, Pressure) must rewrite this
        method
        """
        if not magnitud:
            magnitud = self.__class__.__name__
        self.magnitud = magnitud

        if data is None:
            self._data = 0
            self.code = "n/a"
        else:
            self._data = float(data)
            self.code = ""

        if unit == "conf":
            Config = getMainWindowConfig()
            unit = self.__units__[Config.getint('Units', magnitud)]
        elif not unit:
            unit = self.__units__[0]
        self._data = self._getBaseValue(data, unit, magnitud)

        for key in self.__class__.rates:
            # Reject unused currencies
            if self.__class__.__name__ == "Currency" and \
                    key not in self.__units__ and key not in self._uUnused:
                continue

            self.__setattr__(key, self._data / self.__class__.rates[key])

        logging.debug("%s, %f", self.__class__.__name__, self._data)

    def __new__(cls, data, unit="", magnitud=""):
        """Constructor to let multiple parameter input in float"""
        if not magnitud:
            magnitud = cls.__name__
        if data is None:
            data = 0
        elif unit:
            data = cls._getBaseValue(float(data), unit, magnitud)

        return float.__new__(cls, data)

    @classmethod
    def _getBaseValue(cls, data, unit, magnitud):
        """Convert input data to the base unit"""
        if data is None:
            data = 0
        else:
            data = float(data)

        if unit == "conf":
            Config = getMainWindowConfig()
            unit = cls.__units__[Config.getint('Units', magnitud)]
        elif not unit:
            unit = cls.__units__[0]

        try:
            conversion = cls.rates[unit]
        except KeyError:
            raise ValueError(
                QApplication.translate("pychemqt", "Wrong input code"))

        data *= conversion
        return data

    def __add__(self, other):
        """Support for += operation"""
        return self.__class__(self._data+other)

    def __sub__(self, other):
        """Support for -= operation"""
        return self.__class__(self._data-other)

    def config(self, magnitud=""):
        """Using config file return the value in the configurated unit"""
        if not magnitud:
            magnitud = self.__class__.__name__
        Config = getMainWindowConfig()
        value = Config.getint('Units', magnitud)
        return self.__getattribute__(self.__units__[value])

    @classmethod
    def text(cls, magnitud=""):
        """Using config file return the configurated unit text"""
        if not magnitud:
            magnitud = cls.__name__
        Config = getMainWindowConfig()
        return cls.__text__[Config.getint("Units", magnitud)]

    @classmethod
    def func(cls, magnitud=""):
        """Return the configurated unit name for getattribute call"""
        if not magnitud:
            magnitud = cls.__name__
        Config = getMainWindowConfig()
        return cls.__units__[Config.getint("Units", magnitud)]

    @classmethod
    def magnitudes(cls):
        """Return the magnitudes list for unit,
        if define several magnitudes, it must fill the _magnitudes variable"""
        if cls._magnitudes:
            return cls._magnitudes
        else:
            return [(cls.__name__, cls.__title__)]

    def format(self, unit="", magnitud=""):
        """Using config file return the unit value in desired numeric format"""
        if not magnitud:
            magnitud = self.__class__.__name__
        if not unit:
            unit = self.func(magnitud)
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        kwargs = eval(Preferences.get("NumericFormat", magnitud))
        value = self.__getattribute__(unit)
        return representacion(value, **kwargs)

    def get_str(self, conf=None):
        """Return a string representation of class"""
        if self.code:
            return self.code
        else:
            if not conf:
                conf = self.func(self.magnitud)
                txt = self.text(self.magnitud)
            else:
                txt = self.__text__[self.__units__.index(conf)]
            num = self.format(conf, self.magnitud)
            return num+" "+txt
    str = property(get_str)


class Dimensionless(float):
    """Dummy class to integrate dimensionless magnitudes
with support for class unidad operations: txt, config. func."""
    __title__ = QApplication.translate("pychemqt", "Dimensionless")
    __text__ = []
    _magnitudes = []

    def __init__(self, data, txt=""):

        self.txt = txt

        if data is None:
            self._data = 0
            self.code = "n/a"
        else:
            self._data = float(data)
            self.code = ""
        logging.debug("%s, %f" % (self.__class__.__name__, self._data))

    def __new__(cls, data, txt=""):
        """Discard superfluous parameters for this class"""
        if data is None:
            data = 0
        return float.__new__(cls, data)

    @classmethod
    def text(cls):
        return ""

    @classmethod
    def func(cls):
        return ""

    def config(self):
        return self

    def format(self, unit):
        """Using config file return the unit value in desired numeric format"""
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        kwargs = eval(Preferences.get("NumericFormat", "Dimensionless"))
        return representacion(self, **kwargs)

    @property
    def str(self):
        num = self.format(self)
        if self.txt:
            num += " "+self.txt
        return num


class Temperature(unidad):
    __title__ = QApplication.translate("pychemqt", "Temperature")
    __text__ = ['K', 'ºC', 'ºR', 'ºF', 'ºRe']
    __units__ = ['K', 'C', 'R', 'F', 'Re']
    __tooltip__ = ['Kelvin', 'Celsius', 'Rankine', 'Fahrenheit', 'Reaumur']
    __units_set__ = {"altsi": "C", "si": "K", "metric": "C", "cgs": "C",
                     "english": "F"}
    __test__ = [{"input": {"value": 25, "unit": "C"},
                 "prop": {"K": 298.15, "C": 25, "F": 77}}]

    def __init__(self, data, unit="K", magnitud=""):

        if not magnitud:
            magnitud = self.__class__.__name__
        self.magnitud = magnitud

        if data is None:
            data = 0
            self.code = "n/a"
        else:
            data = float(data)
            self.code = ""

        if unit == "conf":
            Config = getMainWindowConfig()
            unit = self.__units__[Config.getint('Units', magnitud)]

        if unit == "K":
            self._data = data
        elif unit == "C":
            self._data = C2K(data)
        elif unit == "F":
            self._data = F2K(data)
        elif unit == "R":
            self._data = R2K(data)
        elif unit == "Re":
            self._data = Re2K(data)
        else:
            raise ValueError(
                QApplication.translate("pychemqt", "Wrong input code"))

        self.K = self._data
        self.C = K2C(self._data)
        self.F = K2F(self._data)
        self.R = K2R(self._data)
        self.Re = K2Re(self._data)
        logging.debug("%s, %f" % (self.__class__.__name__, self._data))

    @classmethod
    def _getBaseValue(cls, data, unit, magnitud):
        if data is None:
            data = 0
        if not magnitud:
            magnitud = cls.__name__

        if unit == "conf":
            Config = getMainWindowConfig()
            unit = cls.__units__[Config.getint('Units', magnitud)]
        elif not unit:
            unit = "K"

        if unit == "C":
            data = C2K(data)
        elif unit == "F":
            data = F2K(data)
        elif unit == "R":
            data = R2K(data)
        elif unit == "Re":
            data = Re2K(data)
        elif unit != "K":
            raise ValueError(
                QApplication.translate("pychemqt", "Wrong input code"))

        return data


class DeltaT(unidad):
    __title__ = QApplication.translate("pychemqt", "Temperature increase")
    rates = {"K": 1.,
             "C": 1.,
             "F": k.Rankine,
             "R": k.Rankine,
             "Re": k.Reaumur}
    __text__ = ['K', 'ºC', 'ºR', 'ºF', 'ºRe']
    __units__ = ['K', 'C', 'R', 'F', 'Re']
    __tooltip__ = ['Kelvin', 'Celsius', 'Rankine', 'Fahrenheit', 'Reaumur']
    __units_set__ = {"altsi": "C", "si": "K", "metric": "C", "cgs": "C",
                     "english": "F"}
    __test__ = [{"input": {"value": 25, "unit": "C"},
                 "prop": {"K": 25, "F": 45}}]


class Angle(unidad):
    __title__ = QApplication.translate("pychemqt", "Angle")
    rates = {"rad": 1.,
             "deg": 2*k.pi/360,
             "min": 2*k.pi/360/60,
             "sec": 2*k.pi/360/3600,
             "grad": 2*k.pi/400}
    __text__ = ["rad", "º deg", "'", '"', "grad"]
    __units__ = ["rad", "deg", "min", "sec", "grad"]
    __tooltip__ = [QApplication.translate("pychemqt", "Radian"),
                   QApplication.translate("pychemqt", "Degree"),
                   QApplication.translate("pychemqt", "Arcminute"),
                   QApplication.translate("pychemqt", "Arcsecond"),
                   QApplication.translate("pychemqt", "Gradian")]
    __units_set__ = {"altsi": "rad", "si": "rad", "metric": "rad",
                     "cgs": "rad", "english": "rad"}
    __test__ = [{"input": {"value": 25, "unit": "deg"},
                 "prop": {"rad": 0.436332312999}}]


class Length(unidad):
    __title__ = QApplication.translate("pychemqt", "Length")
    rates = {"m": 1.,
             "cm": k.centi,
             "mm": k.milli,
             "micra": k.micro,
             "km": k.kilo,
             "inch": k.inch,
             "ft": k.foot,
             "yd": k.yard,
             "milla": k.milla,
             "milla_nau": k.milla_nau,
             "pm": k.pico,
             "A": 1e-10}
    __text__ = ['m', 'cm', 'mm', 'µm', 'km', 'inch', 'ft', 'yard', 'milla',
                "M", "pm", "Å"]
    __units__ = ['m', 'cm', 'mm', 'micra', 'km', 'inch', 'ft', 'yd', 'milla',
                 "milla_nau", "pm", "A"]
    __tooltip__ = [QApplication.translate("pychemqt", "meter"),
                   QApplication.translate("pychemqt", "centimeter"),
                   QApplication.translate("pychemqt", "milimeter"),
                   QApplication.translate("pychemqt", "micra"),
                   QApplication.translate("pychemqt", "kilometer"),
                   QApplication.translate("pychemqt", "inch"),
                   QApplication.translate("pychemqt", "foot"),
                   QApplication.translate("pychemqt", "yard"),
                   QApplication.translate("pychemqt", "mile"),
                   QApplication.translate("pychemqt", "nautical mile"),
                   QApplication.translate("pychemqt", "picometer"),
                   "Ångström"]
    _magnitudes = [
        ("Length", QApplication.translate("pychemqt", "Length")),
        ("ParticleDiameter", QApplication.translate("pychemqt",
                                                    "Particle Diameter")),
        ("Thickness", QApplication.translate("pychemqt", "Thickness")),
        ("PipeDiameter", QApplication.translate("pychemqt", "Pipe Diameter")),
        ("Head", QApplication.translate("pychemqt", "Head"))]
    __units_set__ = {
        "Length": {"altsi": "m", "si": "m", "metric": "m", "cgs": "cm",
                   "english": "ft"},
        "ParticleDiameter": {"altsi": "mm", "si": "mm", "metric": "mm",
                             "cgs": "cm", "english": "inch"},
        "Thickness": {"altsi": "mm", "si": "mm", "metric": "mm", "cgs": "cm",
                      "english": "inch"},
        "PipeDiameter": {"altsi": "mm", "si": "mm", "metric": "cm",
                         "cgs": "cm", "english": "inch"},
        "Head": {"altsi": "m", "si": "m", "metric": "m", "cgs": "cm",
                 "english": "ft"}}
    __test__ = [{"input": {"value": 12, "unit": "inch"},
                 "prop": {"m": 0.3048, "inch": 12, "ft": 1}}]


class Area(unidad):
    __title__ = QApplication.translate("pychemqt", "Area")
    rates = {"m2": 1.,
             "cm2": k.centi**2,
             "mm2": k.milli**2,
             "km2": k.kilo**2,
             "inch2": k.inch**2,
             "ft2": k.foot**2,
             "yd2": k.yard**2,
             "ha": k.hecto**2,
             "acre": k.acre}
    __text__ = ['m²', 'cm²', 'mm²', "km2", 'inch²', 'ft²', 'yd²',
                'ha', 'acre']
    __units__ = ['m2', 'cm2', 'mm2', 'km2', 'inch2', 'ft2', 'yd2', 'ha',
                 "acre"]
    __tooltip__ = [QApplication.translate("pychemqt", "square meter"),
                   QApplication.translate("pychemqt", "square centimeter"),
                   QApplication.translate("pychemqt", "square milimeter"),
                   QApplication.translate("pychemqt", "square kilometer"),
                   QApplication.translate("pychemqt", "square inch"),
                   QApplication.translate("pychemqt", "square foot"),
                   QApplication.translate("pychemqt", "square yard"),
                   QApplication.translate("pychemqt", "hectarea"),
                   QApplication.translate("pychemqt", "acre")]
    __units_set__ = {"altsi": "m2", "si": "m2", "metric": "m2", "cgs": "cm2",
                     "english": "ft2"}
    __test__ = [{"input": {"value": 1, "unit": "ft2"},
                 "prop": {"m2": 0.09290304, "inch2": 144}}]


class Volume(unidad):
    __title__ = QApplication.translate("pychemqt", "Volume")
    rates = {"m3": 1.,
             "cc": k.centi**3,
             "l": k.deci**3,
             "ml": k.deci**3*k.milli,
             "km3": k.kilo**3,
             "inch3": k.inch**3,
             "ft3": k.foot**3,
             "yd3": k.yard**3,
             "galUS": k.gallon,
             "galUK": k.gallon_imp,
             "qtUSliq": k.qtUSliq,
             "qtUSdr": k.qtUSdry,
             "qtUK": k.qtUK,
             "bbl": k.bbl,
             "onz": k.fluid_ounce,
             "onzUK": k.fluid_ounce_imp,
             "bblUK": 36*k.gallon_imp,
             "bblUS": 31.5*k.gallon}
    __text__ = ['m³', 'cm³', 'liter', 'yd³', 'ft³', 'inch³', 'galon US',
                'galon UK', 'oil', 'bbl', 'bblUK', 'onza', 'onza UK']
    __units__ = ['m3', 'cc', 'l', 'yd3', 'ft3', 'inch3', 'galUS', 'galUK',
                 'bbl', 'bblUS', "bblUK", 'onz', 'onzUK']
    __tooltip__ = [
        'm³', 'cm³',
        QApplication.translate("pychemqt", "liter"),
        'yd³', 'ft³', 'inch³',
        QApplication.translate("pychemqt", "US liquid gallon"),
        QApplication.translate("pychemqt", "Imperial gallon"),
        QApplication.translate("pychemqt", "US fluid barrel"),
        QApplication.translate("pychemqt", "UK fluid barrel"),
        QApplication.translate("pychemqt", "Oil barrel"),
        QApplication.translate("pychemqt", "US customary fluid ounce"),
        QApplication.translate("pychemqt", "Imperial fluid ounce")]
    _magnitudes = [
        ("Volume", QApplication.translate("pychemqt", "Volume")),
        ("VolLiq", QApplication.translate("pychemqt", "Liquid Volume")),
        ("VolGas", QApplication.translate("pychemqt", "Gas Volume"))]
    __units_set__ = {
        "Volume": {"altsi": "m3", "si": "m3", "metric": "m3", "cgs": "cc",
                   "english": "ft3"},
        "VolLiq": {"altsi": "m3", "si": "m3", "metric": "m3", "cgs": "cc",
                   "english": "ft3"},
        "VolGas": {"altsi": "m3", "si": "m3", "metric": "m3", "cgs": "cc",
                   "english": "ft3"}}
    __test__ = [{"input": {"value": 1, "unit": "bbl"},
                 "prop": {"l": 158.987294928, "ft3": 5.61458333333,
                          "galUS": 42}}]


class Time(unidad):
    __title__ = QApplication.translate("pychemqt", "Time")
    rates = {"s": 1.,
             "min": k.minute,
             "h": k.hour,
             "day": k.day,
             "year": k.year}
    __text__ = ['s', 'min', 'h', 'day', 'year']
    __units__ = ['s', 'min', 'h', 'day', 'year']
    __tooltip__ = [QApplication.translate("pychemqt", "second"),
                   QApplication.translate("pychemqt", "minute"),
                   QApplication.translate("pychemqt", "hour"),
                   QApplication.translate("pychemqt", "day"),
                   QApplication.translate("pychemqt", "year")]
    __units_set__ = {"altsi": "h", "si": "h", "metric": "h", "cgs": "s",
                     "english": "h"}
    __test__ = [{"input": {"value": 1, "unit": "day"},
                 "prop": {"min": 1440, "h": 24}}]


class Frequency(unidad):
    __title__ = QApplication.translate("pychemqt", "Frequency")
    rates = {"rpm": 1.,
             "rph": 1./60,
             "rps": 60.,
             "Hz": 60.,
             "rads": 30./k.pi,
             "radmin": 30*60./k.pi,
             "radh": 30*3600./k.pi}
    __text__ = ['rpm', 'rps', 'rph', 'Hz', 'rad/s', 'rad/min', 'rad/hr']
    __units__ = ['rpm', 'rps', 'rph', 'Hz', 'rads', 'radmin', "radh"]
    __units_set__ = {"altsi": "rpm", "si": "Hz", "metric": "Hz", "cgs": "Hz",
                     "english": "rpm"}
    __test__ = [{"input": {"value": 1, "unit": "rads"},
                 "prop": {"rpm": 9.54929658551}}]


class Speed(unidad):
    __title__ = QApplication.translate("pychemqt", "Speed")
    rates = {"ms": 1.,
             "cms": k.centi,
             "mms": k.milli,
             "kms": k.kilo,
             "mmin": 1./k.minute,
             "kmmin": k.kilo/k.minute,
             "kmh": k.kilo/k.hour,
             "mday": 1./k.day,
             "kmday": k.kilo/k.day,
             "fts": k.foot,
             "ftmin": k.foot/k.minute,
             "fth": k.foot/k.hour,
             "ftday": k.foot/k.day,
             "inchs": k.inch,
             "inchmin": k.inch/k.minute,
             "mph": k.mile/k.hour,
             "kt": k.nautical_mile/k.hour}
    __text__ = ['m/s', 'cm/s', 'mm/s', 'km/s', 'ft/s', 'ft/min', 'm/min',
                'km/min', 'km/h', 'km/day', 'mph', 'nudo']
    __units__ = ['ms', 'cms', 'mms', 'kms', 'fts', 'ftmin', 'mmin', 'kmmin',
                 'kmh', 'kmday', 'mph', 'kt']
    __tooltip__ = ['m/s', 'cm/s', 'mm/s', 'km/s', 'ft/s', 'ft/min', 'm/min',
                   'km/min', 'km/h', 'km/day', 'mph',
                   QApplication.translate("pychemqt", "Knot")]
    __units_set__ = {"altsi": "ms", "si": "ms", "metric": "ms", "cgs": "cms",
                     "english": "fts"}
    __test__ = [{"input": {"value": 1, "unit": "ms"},
                 "prop": {"mmin": 60, "kmh": 3.6, "fts": 3.28083989501}}]


class Acceleration(unidad):
    __title__ = QApplication.translate("pychemqt", "Acceleration")
    rates = {"ms2": 1.,
             "cms2": k.centi,
             "fts2": k.foot,
             "inchs2": k.inch,
             "yds2": k.yard,
             "mmin2": 1./k.minute**2,
             "cmmin2": k.centi/k.minute**2,
             "ftmin2": k.foot/k.minute**2,
             "inchmin2": k.inch/k.minute**2,
             "ydmin2": k.yard/k.minute**2}
    __text__ = ["m/s²", "cm/s²", "ft/s²", "inch/s²", "yd/s²", "m/min²",
                "cm/min²", "ft/min²", "inch/min²"]
    __units__ = ["ms2", "cms2", "fts2", "inchs2", "yds2", "mmin2", "cmmin2",
                 "ftmin2", "inchmin2"]
    __units_set__ = {"altsi": "ms2", "si": "ms2", "metric": "ms2",
                     "cgs": "cms2", "english": "fts2"}
    __test__ = [{"input": {"value": 9.81, "unit": "ms2"},
                 "prop": {"fts2": 32.1850393701}}]


class Mass(unidad):
    __title__ = QApplication.translate("pychemqt", "Mass")
    rates = {"kg": 1.,
             "g": 1./k.kilo,
             "mg": 1./k.mega,
             "Ton": k.kilo,
             "lb": k.pound,
             "grain": k.grain,
             "oz": k.oz,
             "slug": k.slug,
             "cwtUK": k.cwtUK,
             "cwtUS": k.cwtUS,
             "TonUK": k.tonUK,
             "TonUS": k.tonUS}
    __text__ = ['kg', 'g', 'mg', 'Ton', 'lb', 'grano', 'onza', 'slug', 'TonUK',
                'TonUS']
    __units__ = ['kg', 'g', 'mg', 'Ton', 'lb', 'grain', 'oz', 'slug', 'TonUK',
                 'TonUS']
    __tooltip__ = [QApplication.translate("pychemqt", "kilogram"),
                   QApplication.translate("pychemqt", "gram"),
                   QApplication.translate("pychemqt", "miligram"),
                   QApplication.translate("pychemqt", "ton"),
                   QApplication.translate("pychemqt", "pound"),
                   QApplication.translate("pychemqt", "grain"),
                   QApplication.translate("pychemqt", "ounce"),
                   QApplication.translate("pychemqt", "slug"),
                   QApplication.translate("pychemqt", "long ton (UK)"),
                   QApplication.translate("pychemqt", "short ton (US)")]
    __units_set__ = {"altsi": "kg", "si": "kg", "metric": "kg", "cgs": "g",
                     "english": "lb"}
    __test__ = [{"input": {"value": 1, "unit": "lb"},
                 "prop": {"kg": 0.45359237, "g": 453.59237, "oz": 16}}]


class Mol(unidad):
    __title__ = QApplication.translate("pychemqt", "Mol")
    rates = {"kmol": 1.,
             "mol": 1./k.kilo,
             "milimol": 1./k.mega,
             "lbmol": k.pound}
    __text__ = ['kmol', 'mol', 'mmol', "lbmol"]
    __units__ = ['kmol', 'mol', 'mmol', "lbmol"]
    __units_set__ = {"altsi": "kmol", "si": "kmol", "metric": "kmol",
                     "cgs": "mol", "english": "lbmol"}
    __test__ = [{"input": {"value": 1, "unit": "kmol"},
                 "prop": {"mol": 1000, "lbmol": 2.20462262185}}]


class SpecificVolume(unidad):
    __title__ = QApplication.translate("pychemqt", "Specific Volume")
    rates = {"m3kg": 1.,
             "lg": 1.,
             "lkg": k.liter,
             "ccg": k.liter,
             "mlg": k.liter,
             "m3g": 1/k.gram,
             "cckg": k.micro,
             "ft3lb": k.foot**3/k.pound,
             "inch3lb": k.inch**3/k.pound,
             "galUKlb": k.gallon_imp/k.pound,
             "galUSlb": k.gallon/k.pound,
             "bbllb": k.bbl/k.pound,
             "ft3tonUK": k.foot**3/k.tonUK,
             "ft3tonUS": k.foot**3/k.tonUS,
             "ft3slug": k.foot**3/k.slug,
             "ft3oz": k.foot**3/k.oz,
             "in3oz": k.inch**3/k.oz,
             "galUKoz": k.gallon_imp/k.oz,
             "galUSoz": k.gallon/k.oz}
    __text__ = ['m³/kg', 'cm³/g', 'ml/g', 'm³/g', 'cm³/kg', 'ft³/lb',
                'in³/lb', 'gallon UK/lb', 'gallon US/lb', 'barril/lb',
                'ft³/ton UK', 'ft³/ton US', 'ft³/slug', 'ft³/onza',
                'in³/onza', 'gallon UK/onza', 'gallon US/onza']
    __units__ = ['m3kg', 'ccg', 'mlg', 'm3g', 'cckg', 'ft3lb', 'inch3lb',
                 'galUKlb', 'galUSlb', 'bbllb', 'ft3tonUK', 'ft3tonUS',
                 'ft3slug', 'ft3oz', 'in3oz', 'galUKoz', 'galUSoz']
    __units_set__ = {"altsi": "m3kg", "si": "m3kg", "metric": "m3kg",
                     "cgs": "ccg", "english": "ft3lb"}
    __test__ = [{"input": {"value": 50, "unit": "lkg"},
                 "prop": {"m3kg": 0.05, "ft3lb": 0.800923168698}}]


class SpecificVolume_square(unidad):
    __title__ = QApplication.translate("pychemqt", "Third virial coefficient")
    rates = {"m3kg": 1.,
             "lg": 1.,
             "lkg": k.liter**2,
             "ccg": k.liter**2,
             "mlg": k.liter**2,
             "m3g": 1/k.gram**2,
             "cckg": k.micro**2,
             "ft3lb": k.foot**6/k.pound**2,
             "inch3lb": k.inch**6/k.pound**2,
             "galUKlb": k.gallon_imp**2/k.pound**2,
             "galUSlb": k.gallon**2/k.pound**2,
             "bbllb": k.bbl**2/k.pound**2,
             "ft3tonUK": k.foot**6/k.tonUK**2,
             "ft3tonUS": k.foot**6/k.tonUS**2,
             "ft3slug": k.foot**6/k.slug**2,
             "ft3oz": k.foot**6/k.oz**2,
             "in3oz": k.inch**6/k.oz**2,
             "galUKoz": k.gallon_imp**2/k.oz**2,
             "galUSoz": k.gallon**2/k.oz**2}
    __text__ = ['m⁶/kg²', 'cm⁶/g²', 'ml²/g²', 'm⁶/g²', 'cm⁶/kg²', 'ft⁶/lb²',
                'in⁶/lb²', 'gallon UK²/lb²', 'gallon US²/lb²', 'barril²/lb²',
                'ft⁶/ton UK²', 'ft⁶/ton US²', 'ft⁶/slug²', 'ft⁶/onza²',
                'in⁶/onza²', 'gallon UK²/onza²', 'gallon US²/onza²']
    __units__ = ['m3kg', 'ccg', 'mlg', 'm3g', 'cckg', 'ft3lb', 'inch3lb',
                 'galUKlb', 'galUSlb', 'bbllb', 'ft3tonUK', 'ft3tonUS',
                 'ft3slug', 'ft3oz', 'in3oz', 'galUKoz', 'galUSoz']
    __units_set__ = {"altsi": "m3kg", "si": "m3kg", "metric": "m3kg",
                     "cgs": "ccg", "english": "ft3lb"}
    __test__ = [{"input": {"value": 50, "unit": "lkg"},
                 "prop": {"m3kg": 5e-5, "ft3lb": 0.0128295584431}}]

# TODO: Add unit for fourth virial coefficient, only useful for refprop library
# when work in qt loop


class MolarVolume(unidad):
    __title__ = QApplication.translate("pychemqt", "Molar Volume")
    rates = {"m3kmol": 1.,
             "lmol": 1.,
             "lkmol": k.liter,
             "ccmol": k.liter,
             "mlmol": k.liter,
             "m3mol": 1/k.gram,
             "cckmol": k.micro,
             "ft3lbmol": k.foot**3/k.pound,
             "inch3lbmol": k.inch**3/k.pound}
    __text__ = ['m³/kmol', 'l/mol', 'l/kmol', 'cm³/mol', 'ml/mol',
                'm³/mol', 'cm³/kmol', 'ft³/lbmol', 'in³/lbmol']
    __units__ = ['m3kmol', 'lmol', 'lkmol', 'ccmol', 'mlmol', 'm3mol',
                 'cckmol', 'ft3lbmol', 'inch3lbmol']
    __units_set__ = {"altsi": "m3kmol", "si": "m3kmol", "metric": "m3kmol",
                     "cgs": "ccmol", "english": "ft3lbmol"}
    __test__ = [{"input": {"value": 50, "unit": "lkmol"},
                 "prop": {"m3kmol": 0.05, "ft3lbmol": 0.800923168698}}]


class Density(unidad):
    __title__ = QApplication.translate("pychemqt", "Density")
    rates = {"kgm3": 1.,
             "gl": 1.,
             "kgl": 1./k.liter,
             "gcc": 1./k.liter,
             "gml": 1./k.liter,
             "gm3": k.gram,
             "kgcc": 1./k.micro,
             "lbft3": k.pound/k.foot**3,
             "lbin3": k.pound/k.inch**3,
             "lbgalUK": k.pound/k.gallon_imp,
             "lbgalUS": k.pound/k.gallon,
             "lbbbl": k.pound/k.bbl,
             "tonUKft3": k.tonUK/k.foot**3,
             "tonUSft3": k.tonUS/k.foot**3,
             "slugft3": k.slug/k.foot**3,
             "ozft3": k.oz/k.foot**3,
             "ozin3": k.oz/k.inch**3,
             "ozgalUK": k.oz/k.gallon_imp,
             "ozgalUS": k.oz/k.gallon}
    __text__ = ['kg/m³', 'g/cm³', 'g/m³', 'kg/cm³', 'lb/ft³', 'lb/inch³',
                'lb/galon UK', 'lb/galon US', 'lb/barril', 'ton UK/ft³',
                'ton US/ft³', 'slug/ft³', 'onza/ft³', 'onza/inch³',
                'onza/galon UK', 'onza/galon US']
    __units__ = ['kgm3', 'gcc', 'gm3', 'kgcc', 'lbft3', 'lbin3', 'lbgalUK',
                 'lbgalUS', 'lbbbl', 'tonUKft3', 'tonUSft3', 'slugft3',
                 'ozft3', 'ozin3', 'ozgalUK', 'ozgalUS']
    _magnitudes = [
        ("Density", QApplication.translate("pychemqt", "Density")),
        ("DenLiq", QApplication.translate("pychemqt", "Liquid Density")),
        ("DenGas", QApplication.translate("pychemqt", "Gas Density"))]
    __units_set__ = {
        "Density": {"altsi": "kgm3", "si": "kgm3", "metric": "kgm3",
                    "cgs": "gcc", "english": "lbft3"},
        "DenLiq": {"altsi": "kgm3", "si": "kgm3", "metric": "kgm3",
                   "cgs": "gcc", "english": "lbft3"},
        "DenGas": {"altsi": "kgm3", "si": "kgm3", "metric": "kgm3",
                   "cgs": "gcc", "english": "lbft3"}}
    __test__ = [{"input": {"value": 1, "unit": "kgl"},
                 "prop": {"kgm3": 1000, "lbft3": 62.4279605761}}]


class MolarDensity(unidad):
    __title__ = QApplication.translate("pychemqt", "Molar Density")
    rates = {"kmolm3": 1.,
             "moll": 1.,
             "molcc": 1./k.liter,
             "kmoll": 1./k.liter,
             "molm3": k.gram,
             "kmolcc": 1./k.micro,
             "lbmolft3": k.pound/k.foot**3,
             "lbmolin3": k.pound/k.inch**3}
    __text__ = ['kmol/m³', 'mol/cm³', 'mol/m³', 'kmol/cm³', 'lbmol/ft³',
                'lbmol/inch³']
    __units__ = ['kmolm3', 'molcc', 'molm3', 'kmolcc', 'lbmolft3', 'lbmolin3']
    __units_set__ = {"altsi": "kmolm3", "si": "kmolm3", "metric": "kmolm3",
                     "cgs": "molcc", "english": "lbmolft3"}
    __test__ = [{"input": {"value": 1, "unit": "kmolm3"},
                 "prop": {"molcc": 0.001, "lbmolft3": 0.0624279605761}}]


class Force(unidad):
    __title__ = QApplication.translate("pychemqt", "Force")
    rates = {"N": 1.,
             "kN": k.kilo,
             "dyn": k.dyn,
             "kgf": k.kgf,
             "gf": k.kgf/k.kilo,
             "lbf": k.lbf,
             "ozf": k.ozf,
             "pdl": k.pdl,
             "TonfUK": k.TonfUK,
             "TonfUS": k.TonfUS}
    __text__ = ["N", "kN", "dyn", "kgf", "gf", "lbf", "ozf", "Poundal",
                "TonfUK", "TonfUS"]
    __units__ = ["N", "kN", "dyn", "kgf", "gf", "lbf", "ozf", "pdl", "TonfUK",
                 "TonfUS"]
    __tooltip__ = ["Newton", "Kilonewton", "Dyna",
                   QApplication.translate("pychemqt", "Kilogram force"),
                   QApplication.translate("pychemqt", "Gram force"),
                   QApplication.translate("pychemqt", "Pound force"),
                   QApplication.translate("pychemqt", "Ounze force"),
                   "Poundal", "TonfUK", "TonfUS"]
    __units_set__ = {"altsi": "kN", "si": "N", "metric": "N", "cgs": "dyn",
                     "english": "lbf"}
    __test__ = [{"input": {"value": 1, "unit": "pdl"},
                 "prop": {"N": 0.138254954376, "kgf": 0.0140980818502,
                          "dyn": 13825.4954376}}]


class Pressure(unidad):
    __title__ = QApplication.translate("pychemqt", "Pressure")
    rates = {"Pa": 1.,
             "MPa": k.mega,
             "hPa": k.hecto,
             "kPa": k.kilo,
             "bar": k.bar,
             "baria": 0.1,
             "mbar": k.bar/k.kilo,
             "psi": k.psi,
             "atm": k.atm,
             "kgcm2": k.g/k.centi**2,
             "mmH2O": k.g,
             "mH2O": k.g*k.kilo,
             "cmH2O": k.g*10,
             "inH2O": k.g*k.kilo*k.inch,
             "ftH2O": k.g*k.kilo*k.foot,
             "mmHg": k.torr,
             "cmHg": k.torr*10,
             "inHg": k.torr*k.inch*k.kilo,
             "ftHg": k.torr*k.foot*k.kilo,
             "torr": k.torr,
             "lbcm2": k.g*k.pound/k.centi**2,
             "lbft2": k.g*k.pound/k.foot**2,
             "dyncm2": k.dyn/k.centi**2}
    __text__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'bar g', 'mbar', 'psi',
                'psi g', 'atm', 'kg/cm²', 'kg/cm² g', 'mmH2O', 'cmH2O',
                'mH2O', 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg',
                'lb/cm²', 'lb/ft²', 'dyn/cm²']
    __units__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'barg', 'mbar', 'psi',
                 'psig', 'atm', 'kgcm2', 'kgcm2g', 'mmH2O', 'cmH2O', 'mH2O',
                 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg', 'lbcm2',
                 'lbft2', 'dyncm2']
    __tooltip__ = ["Pascal", "Hectopascal", "Kilopascal", "Megapascal", "bar",
                   QApplication.translate("pychemqt", "Bar gauge"),
                   "Milibar",
                   QApplication.translate("pychemqt", "Pound per square inch"),
                   QApplication.translate(
                       "pychemqt", "Pound per square inch gauge"),
                   QApplication.translate("pychemqt", "Atmosphere"),
                   QApplication.translate(
                       "pychemqt", "Atmosphere technical, kg/cm²"),
                   QApplication.translate(
                       "pychemqt", "Atmosphere technical gauge, kg/cm²g"),
                   QApplication.translate(
                       "pychemqt", "Milimeter of water column"),
                   QApplication.translate(
                       "pychemqt", "Centimeter of water column"),
                   QApplication.translate("pychemqt", "Meter of water column"),
                   QApplication.translate("pychemqt", "Inch of water column"),
                   QApplication.translate("pychemqt", "Foot of water column"),
                   QApplication.translate(
                       "pychemqt", "Milimeter of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Centimeter of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Inch of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Foot of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Pound per square centimeter"),
                   QApplication.translate("pychemqt", "Pound per square foot"),
                   QApplication.translate(
                       "pychemqt", "Dyn per square centimeter")]
    __units_set__ = {"altsi": "bar", "si": "Pa", "metric": "Pa",
                     "cgs": "dyncm2", "english": "psi"}
    __test__ = [{"input": {"value": 760, "unit": "mmHg"},
                 "prop": {"bar": 1.01325, "atm": 1, "psi": 14.6959487755,
                          "kgcm2g": 0}}]

    def __init__(self, data, unit="Pa", magnitud=""):

        if data is None:
            data = 0
            self.code = "n/a"
        else:
            self.code = ""
            data = float(data)

        if not magnitud:
            magnitud = self.__class__.__name__
        self.magnitud = magnitud

        if unit == "conf":
            Config = getMainWindowConfig()
            unit = self.__units__[Config.getint('Units', magnitud)]

        if unit == "barg":
            self._data = data*k.bar+k.atm
        elif unit == "psig":
            self._data = data*k.psi+k.atm
        elif unit == "kgcm2g":
            self._data = data*k.g/k.centi**2+k.atm
        else:
            self._data = data * self.__class__.rates[unit]

        for key in self.__class__.rates:
            self.__setattr__(key, self._data/self.__class__.rates[key])

        self.barg = (self.Pa-k.atm)/k.bar
        self.psig = (self.Pa-k.atm)/k.psi
        self.kgcm2g = (self.Pa-k.atm)*k.centi**2/k.g
        logging.debug("%s, %f" % (self.__class__.__name__, self._data))

    @classmethod
    def _getBaseValue(cls, data, unit, magnitud):
        if data is None:
            data = 0
        if not magnitud:
            magnitud = cls.__name__

        if unit == "conf":
            Config = getMainWindowConfig()
            unit = cls.__units__[Config.getint('Units', magnitud)]
        elif not unit:
            unit = "Pa"

        if unit == "barg":
            data = data*k.bar+k.atm
        elif unit == "psig":
            data = data*k.psi+k.atm
        elif unit == "kgcm2g":
            data = data*k.g/k.centi**2+k.atm
        elif unit in cls.rates:
            data = data * cls.rates[unit]
        else:
            raise ValueError(
                QApplication.translate("pychemqt", "Wrong input code"))

        return data


class DeltaP(unidad):
    __title__ = QApplication.translate("pychemqt", "Pressure increase")
    rates = {"Pa": 1.,
             "MPa": k.mega,
             "hPa": k.hecto,
             "kPa": k.kilo,
             "bar": k.bar,
             "barg": k.bar,
             "baria": 0.1,
             "mbar": k.bar/k.kilo,
             "psi": k.psi,
             "psig": k.psi,
             "atm": k.atm,
             "kgcm2": k.g/k.centi**2,
             "kgcm2g": k.g/k.centi**2,
             "mmH2O": k.g,
             "mH2O": k.g*k.kilo,
             "cmH2O": k.g*10,
             "inH2O": k.g*k.kilo*k.inch,
             "ftH2O": k.g*k.kilo*k.foot,
             "mmHg": k.torr,
             "cmHg": k.torr*10,
             "inHg": k.torr*k.inch*k.kilo,
             "ftHg": k.torr*k.foot*k.kilo,
             "torr": k.torr,
             "lbcm2": k.g*k.pound/k.centi**2,
             "lbft2": k.g*k.pound/k.foot**2,
             "dyncm2": k.dyn/k.centi**2}
    __text__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'bar g', 'mbar', 'psi',
                'psi g', 'atm', 'kg/cm²', 'kg/cm² g', 'mmH2O', 'cmH2O',
                'mH2O', 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg',
                'lb/cm²', 'lb/ft²', 'dyn/cm²']
    __units__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'barg', 'mbar', 'psi',
                 'psig', 'atm', 'kgcm2', 'kgcm2g', 'mmH2O', 'cmH2O', 'mH2O',
                 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg', 'lbcm2',
                 'lbft2', 'dyncm2']
    __tooltip__ = ["Pascal", "Hectopascal", "Kilopascal", "Megapascal", "bar",
                   QApplication.translate("pychemqt", "Bar gauge"),
                   "Milibar",
                   QApplication.translate("pychemqt", "Pound per square inch"),
                   QApplication.translate(
                       "pychemqt", "Pound per square inch gauge"),
                   QApplication.translate("pychemqt", "Atmosphere"),
                   QApplication.translate(
                       "pychemqt", "Atmosphere technical, kg/cm²"),
                   QApplication.translate(
                       "pychemqt", "Atmosphere technical gauge, kg/cm²g"),
                   QApplication.translate(
                       "pychemqt", "Milimeter of water column"),
                   QApplication.translate(
                       "pychemqt", "Centimeter of water column"),
                   QApplication.translate("pychemqt", "Meter of water column"),
                   QApplication.translate("pychemqt", "Inch of water column"),
                   QApplication.translate("pychemqt", "Foot of water column"),
                   QApplication.translate(
                       "pychemqt", "Milimeter of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Centimeter of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Inch of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Foot of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Pound per square centimeter"),
                   QApplication.translate("pychemqt", "Pound per square foot"),
                   QApplication.translate(
                       "pychemqt", "Dyn per square centimeter")]
    __units_set__ = {"altsi": "bar", "si": "Pa", "metric": "Pa",
                     "cgs": "dyncm2", "english": "psi"}
    __test__ = [{"input": {"value": 760, "unit": "mmHg"},
                 "prop": {"bar": 1.01325, "atm": 1, "psi": 14.6959487755,
                          "kgcm2": 1.03323}}]


class Energy(unidad):
    __title__ = QApplication.translate("pychemqt", "Energy")
    rates = {"J": 1.,
             "kJ": k.kilo,
             "MJ": k.mega,
             "cal": k.calorie,
             "kcal": k.calorie*k.kilo,
             "cal_i": k.calorie_IT,
             "erg": k.erg,
             "Btu": k.Btu,
             "kBtu": k.Btu*k.kilo,
             "MBtu": k.Btu*k.mega,
             "Wh": k.hour,
             "kWh": k.hour*k.kilo,
             "MWh": k.hour*k.mega,
             "TNT": k.TNT,
             "HPh": k.hp*k.hour,
             "CVh": k.CV*k.hour,
             "kgfm": k.kgf,
             "lbfft": k.lbf*k.foot,
             "GeV": k.eV*k.giga,
             "oil": k.BarrilOil,
             "toe": k.TmOil,
             "tce": k.TmCoal}
    __text__ = ['J', 'kJ', 'MJ', 'cal', 'kcal', 'cal int', 'erg', 'Btu',
                'kBtu', 'MBtu', 'Wh', 'kWh', 'MWh', 'HPh', 'kgf/m', 'lbf/ft',
                'TNT', 'CVh', 'GeV', 'oil', 'toe', 'tce']
    __units__ = ['J', 'kJ', 'MJ', 'cal', 'kcal', 'cal_i', 'erg', 'Btu', 'kBtu',
                 'MBtu', 'Wh', 'kWh', 'MWh', 'HPh', 'kgfm', 'lbfft',
                 'TNT', 'CVh', 'GeV', 'oil', 'toe', 'tce']
    _magnitudes = [
        ("Energy", QApplication.translate("pychemqt", "Energy")),
        ("Work", QApplication.translate("pychemqt", "Work"))]
    __tooltip__ = [
        QApplication.translate("pychemqt", "Joule"),
        QApplication.translate("pychemqt", "Kilojoule"),
        QApplication.translate("pychemqt", "Megajoule"),
        QApplication.translate("pychemqt", "Calorie"),
        QApplication.translate("pychemqt", "Kilocalorie"),
        QApplication.translate("pychemqt", "Calorie international"),
        QApplication.translate("pychemqt", "Erg"),
        QApplication.translate("pychemqt", "Btu"),
        QApplication.translate("pychemqt", "KiloBtu"),
        QApplication.translate("pychemqt", "MegaBtu"),
        QApplication.translate("pychemqt", "Watt-hour"),
        QApplication.translate("pychemqt", "Kilowatt-hour"),
        QApplication.translate("pychemqt", "Megawatt-hour"),
        QApplication.translate("pychemqt", "Horsepower-hour"),
        'kgf/m', 'lbf/ft',
        QApplication.translate("pychemqt", "Ton TNT equivalent"),
        QApplication.translate("pychemqt", "Metric horsepower-hour"),
        QApplication.translate("pychemqt", "Gigaelectronvolt"),
        QApplication.translate("pychemqt", "Barrel petrol"),
        QApplication.translate("pychemqt", "Tonne of oil equivalent"),
        QApplication.translate("pychemqt", "Tonne of coal equivalent")]
    __units_set__ = {
        "Energy": {"altsi": "MJ", "si": "MJ", "metric": "J", "cgs": "erg",
                   "english": "MBtu"},
        "Work": {"altsi": "MJ", "si": "kWh", "metric": "J", "cgs": "erg",
                 "english": "HPh"}}
    __test__ = [{"input": {"value": 1, "unit": "kcal"},
                 "prop": {"J": 4184, "Btu": 3.96566683139, "Wh": 1.162222}}]


class Enthalpy(unidad):
    __title__ = QApplication.translate("pychemqt", "Enthalpy")
    rates = {"Jkg": 1.,
             "kJkg": k.kilo,
             "Jg": k.kilo,
             "MJkg": k.mega,
             "kJg": k.mega,
             "kWhkg": k.hour*k.kilo,
             "calkg": k.calorie,
             "kcalkg": k.calorie*k.kilo,
             "calg": k.calorie*k.kilo,
             "callb": k.calorie/k.lb,
             "kcalg": k.calorie*k.mega,
             "Btulb": k.Btu/k.lb}
    __text__ = ['J/kg', 'kJ/kg', 'MJ/kg', 'cal/kg', 'kcal/kg', 'calg',
                'cal/lb', 'Btu/lb']
    __units__ = ['Jkg', 'kJkg', 'MJkg', 'calkg', 'kcalkg', 'calg', 'callb',
                 'Btulb']
    __units_set__ = {"altsi": "kJkg", "si": "Jkg", "metric": "Jkg",
                     "cgs": "calg", "english": "Btulb"}
    __test__ = [{"input": {"value": -5, "unit": "Btulb"},
                 "prop": {"kJkg": -11.63, "kcalkg": -2.77963671128}}]


class MolarEnthalpy(unidad):
    __title__ = QApplication.translate("pychemqt", "Molar Enthalpy")
    rates = {"Jkmol": 1.,
             "kJkmol": k.kilo,
             "Jmol": k.kilo,
             "MJkmol": k.mega,
             "kJmol": k.mega,
             "kWhkmol": k.hour*k.kilo,
             "calkmol": k.calorie,
             "kcalkmol": k.calorie*k.kilo,
             "calmol": k.calorie*k.kilo,
             "callbmol": k.calorie/k.lb,
             "kcalmol": k.calorie*k.mega,
             "Btulbmol": k.Btu/k.lb}
    __text__ = ['J/kmol', 'kJ/kmol', 'MJ/kmol', 'cal/kmol', 'kcal/kmol',
                'calmol', 'cal/lbmol', 'Btu/lbmol']
    __units__ = ['Jkmol', 'kJkmol', 'MJkmol', 'calkmol', 'kcalkmol', 'calmol',
                 'callbmol', 'Btulbmol']
    __units_set__ = {"altsi": "kJkmol", "si": "Jkmol", "metric": "Jkmol",
                     "cgs": "calmol", "english": "Btulbmol"}
    __test__ = [{"input": {"value": -5, "unit": "Btulbmol"},
                 "prop": {"kJkmol": -11.63, "kcalkmol": -2.77963671128}}]


class Entropy(unidad):
    __title__ = QApplication.translate("pychemqt", "Entropy")
    rates = {"JK": 1.,
             "kJK": k.kilo,
             "MJK": k.mega,
             "calK": k.calorie,
             "kcalK": k.calorie*k.kilo,
             "McalK": k.calorie*k.mega,
             "WhK": k.hour,
             "kWhK": k.hour*k.kilo,
             "MWhK": k.hour/k.mega,
             "hphF": k.hour*k.hp/k.Rankine,
             "BtuF": k.Btu/k.Rankine,
             "kBtuF": k.Btu*k.kilo/k.Rankine,
             "MBtuF": k.Btu*k.mega/k.Rankine}
    __text__ = ['J/K', 'kJ/K', 'MJ/K', 'cal/K', 'kcal/K', 'Mcal/K', 'Wh/K',
                'kWh/K', 'MWh/K', 'hph/F', 'Btu/F', 'kBtu/F', 'MBtu/F']
    __units__ = ['JK', 'kJK', 'MJK', 'calK', 'kcalK', 'McalK', 'WhK', 'kWhK',
                 'MWhK', 'hphF', 'BtuF', 'kBtuF', 'MBtuF']
    __units_set__ = {"altsi": "kJK", "si": "JK", "metric": "JK",
                     "cgs": "calK", "english": "MBtuF"}
    __test__ = [{"input": {"value": 30, "unit": "BtuF"},
                 "prop": {"kJK": 56.9730160415, "kcalK": 13.616877639}}]


class SpecificHeat(unidad):
    __title__ = QApplication.translate("pychemqt", "Specific Heat")
    rates = {"JkgK": 1.,
             "kJkgK": k.kilo,
             "JgK": k.kilo,
             "kcalkgK": k.calorie*k.kilo,
             "calgK": k.calorie*k.kilo,
             "kcalgK": k.calorie*k.kilo**2,
             "kWhkgK": k.kilo*k.hour,
             "BtulbF": k.Btu/k.lb/k.Rankine}
    _magnitudes = [
        ("SpecificHeat", QApplication.translate("pychemqt", "Specific Heat")),
        ("SpecificEntropy", QApplication.translate(
            "pychemqt", "Specific Entropy"))]
    __text__ = ['J/kg·K', 'kJ/kg·K', 'kcal/kg·K', 'cal/g·K', 'kcal/g·K',
                'kWh/kg·K', 'Btu/lb·F']
    __units__ = ['JkgK', 'kJkgK', 'kcalkgK', 'calgK', 'kcalgK', 'kWhkgK',
                 'BtulbF']
    __units_set__ = {
        "SpecificHeat": {"altsi": "kJkgK", "si": "JkgK", "metric": "JkgK",
                         "cgs": "calgK", "english": "BtulbF"},
        "SpecificEntropy": {"altsi": "kJkgK", "si": "JkgK", "metric": "JkgK",
                            "cgs": "calgK", "english": "BtulbF"}}
    __test__ = [{"input": {"value": 1, "unit": "BtulbF"},
                 "prop": {"kJkgK": 4.1868, "kcalkgK": 1.00066921606}}]


class MolarSpecificHeat(unidad):
    __title__ = QApplication.translate("pychemqt", "Molar Specific Heat")
    rates = {"JkmolK": 1.,
             "kJkmolK": k.kilo,
             "kJmolK": k.kilo*k.kilo,
             "JmolK": k.kilo,
             "kcalkmolK": k.calorie*k.kilo,
             "calmolK": k.calorie*k.kilo,
             "kcalmolK": k.calorie*k.kilo**2,
             "kWhkmolK": k.kilo*k.hour,
             "BtulbmolF": k.Btu/k.lb/k.Rankine}
    __text__ = ['J/kmol·K', 'kJ/kmol·K', 'kJ/mol·K', 'kcal/kmol·K',
                'cal/mol·K', 'kcal/mol·K', 'kWh/kmol·K', 'Btu/lbmol·F']
    __units__ = ['JkmolK', 'kJkmolK', 'kJmolK', 'kcalkmolK', 'calmolK',
                 'kcalmolK', 'kWhkmolK', 'BtulbmolF']
    __units_set__ = {"altsi": "kJkmolK", "si": "JkmolK", "metric": "JkmolK",
                     "cgs": "calmolK", "english": "BtulbmolF"}
    __test__ = [{"input": {"value": 1, "unit": "BtulbmolF"},
                 "prop": {"kJkmolK": 4.1868, "kcalkmolK": 1.00066921606}}]


class Power(unidad):
    __title__ = QApplication.translate("pychemqt", "Power")
    rates = {"W": 1.,
             "kW": k.kilo,
             "MW": k.mega,
             "hp": k.hp,
             "CV": k.CV,
             "cals": k.calorie,
             "kcalh": k.calorie*k.kilo/k.hour,
             "Jh": 1/k.hour,
             "kJh": k.kilo/k.hour,
             "MJh": k.mega/k.hour,
             "ergs": k.erg,
             "Btus": k.Btu,
             "Btumin": k.Btu/k.minute,
             "Btuh": k.Btu/k.hour,
             "MBtuh": k.Btu/k.hour*k.mega,
             "ftlbfs": k.foot*k.lb*k.g,
             "ftlbfmin": k.foot*k.lb*k.g/k.minute,
             "ftlbfh": k.foot*k.lb*k.g/k.hour}
    __text__ = ['W', 'kW', 'MW', 'hp', 'CV', 'cal/s', 'kcal/h', 'J/h', 'kJ/h',
                'MJ/h', 'erg/s', 'Btu/s', 'Btu/min', 'Btu/h', 'MBtu/h',
                'ft/lbf·s', 'ft/lbf·min', 'ft/lbf·h']
    __units__ = ['W', 'kW', 'MW', 'hp', 'CV', 'cals', 'kcalh', 'Jh', 'kJh',
                 'MJh', 'ergs', 'Btus', 'Btumin', 'Btuh', 'MBtuh', 'ftlbfs',
                 'ftlbfmin', 'ftlbfh']
    __tooltip__ = [
        QApplication.translate("pychemqt", "Watt"),
        QApplication.translate("pychemqt", "Kilowatt"),
        QApplication.translate("pychemqt", "Megawatt"),
        QApplication.translate("pychemqt", "Horsepower"),
        QApplication.translate("pychemqt", "Metric horsepower"),
        QApplication.translate("pychemqt", "Calorie per second"),
        QApplication.translate("pychemqt", "Kilocalorie per hour"),
        QApplication.translate("pychemqt", "Joule per hour"),
        QApplication.translate("pychemqt", "Kilojoule per hour"),
        QApplication.translate("pychemqt", "Megajoule per hour"),
        QApplication.translate("pychemqt", "Erg per second"),
        QApplication.translate("pychemqt", "Btu per second"),
        QApplication.translate("pychemqt", "Btu per minute"),
        QApplication.translate("pychemqt", "Btu per hour"),
        QApplication.translate("pychemqt", "MegaBtu per hour"),
        'ft/lbf·s', 'ft/lbf·min', 'ft/lbf·h']
    _magnitudes = [
        ("EnergyFlow", QApplication.translate("pychemqt", "Energy Flow")),
        ("Power", QApplication.translate("pychemqt", "Power"))]
    __units_set__ = {
        "EnergyFlow": {"altsi": "MJh", "si": "kJh", "metric": "Jh",
                       "cgs": "ergs", "english": "MBtuh"},
        "Power": {"altsi": "hp", "si": "kW", "metric": "Jh", "cgs": "ergs",
                  "english": "hp"}}
    __test__ = [{"input": {"value": 5, "unit": "Btuh"},
                 "prop": {"kW": 0.00146535535086, "hp": 0.00196507389461,
                          "kcalh": 1.26082200361}}]


class MassFlow(unidad):
    __title__ = QApplication.translate("pychemqt", "Mass Flow")
    rates = {"kgs": 1.,
             "kgmin": 1./k.minute,
             "kgh": 1./k.hour,
             "gs": k.milli,
             "gmin": k.milli/k.minute,
             "gh": k.milli/k.hour,
             "Tons": k.kilo,
             "Tonmin": k.kilo/k.minute,
             "Tonh": k.kilo/k.hour,
             "lbs": k.lb,
             "lbmin": k.lb/k.minute,
             "lbh": k.lb/k.hour,
             "TonUKs": k.tonUK,
             "TonUSs": k.tonUS,
             "TonUKmin": k.tonUK/k.minute,
             "TonUSmin": k.tonUS/k.minute,
             "TonUKh": k.tonUK/k.hour,
             "TonUSh": k.tonUS/k.hour,
             "TonUKday": k.tonUK/k.day,
             "TonUSday": k.tonUS/k.day}
    __text__ = ['kg/s', 'kg/min', 'kg/h', 'g/s', 'g/min', 'g/h', 'Ton/s',
                'Ton/min', 'Ton/h', 'lb/s', 'lb/min', 'lb/h', 'TonUK/min',
                'TonUS/min', 'TonUK/h', 'TonUS/h', 'TonUK/day', 'TonUS/day']
    __units__ = ['kgs', 'kgmin', 'kgh', 'gs', 'gmin', 'gh', 'Tons', 'Tonmin',
                 'Tonh', 'lbs', 'lbmin', 'lbh', 'TonUKmin', 'TonUSmin',
                 'TonUKh', 'TonUSh', 'TonUKday', 'TonUSday']
    __units_set__ = {"altsi": "kgh", "si": "kgh", "metric": "kgs", "cgs": "gs",
                     "english": "lbh"}
    __test__ = [{"input": {"value": 1, "unit": "gs"},
                 "prop": {"kgh": 3.6, "lbh": 7.93664143866, "gmin": 60}}]


class MolarFlow(unidad):
    __title__ = QApplication.translate("pychemqt", "Molar Flow")
    rates = {"kmols": 1.,
             "kmolmin": 1./k.minute,
             "kmolh": 1./k.hour,
             "mols": k.milli,
             "molmin": k.milli/k.minute,
             "molh": k.milli/k.hour,
             "lbmols": k.lb,
             "lbmolmin": k.lb/k.minute,
             "lbmolh": k.lb/k.hour}
    __text__ = ['kmol/s', 'kmol/min', 'kmol/h', 'mol/s', 'mol/min', 'mol/h',
                'lbmol/s', 'lbmol/min', 'lbmol/h']
    __units__ = ['kmols', 'kmolmin', 'kmolh', 'mols', 'molmin', 'molh',
                 'lbmols', 'lbmolmin', 'lbmolh']
    __units_set__ = {"altsi": "kmolh", "si": "kmolh", "metric": "kmols",
                     "cgs": "mols", "english": "lbmolh"}
    __test__ = [{"input": {"value": 1, "unit": "mols"},
                 "prop": {"kmolh": 3.6, "lbmolh": 7.93664143866,
                          "molmin": 60}}]


class VolFlow(unidad):
    __title__ = QApplication.translate("pychemqt", "Volumetric Flow")
    rates = {"m3s": 1.,
             "m3min": 1./k.minute,
             "m3h": 1./k.hour,
             "ls": k.liter,
             "lmin": k.liter/k.minute,
             "lh": k.liter/k.hour,
             "ccs": k.micro,
             "cm3s": k.micro,
             "ccmin": k.micro/k.minute,
             "cch": k.micro/k.hour,
             "ft3s": k.foot**3,
             "ft3min": k.foot**3/k.minute,
             "kft3min": k.foot**3*k.kilo/k.minute,
             "ft3h": k.foot**3/k.hour,
             "mft3day": k.foot**3*k.mega/k.day,
             "galUKh": k.gallon_imp/k.hour,
             "galUSh": k.gallon/k.hour,
             "galUKmin": k.gallon_imp/k.minute,
             "galUSmin": k.gallon/k.minute,
             "galUKs": k.gallon_imp,
             "galUSs": k.gallon,
             "bbls": k.bbl,
             "bblmin": k.bbl/k.minute,
             "bblh": k.bbl/k.hour,
             "bblday": k.bbl/k.day}
    __text__ = ['m³/s', 'm³/min', 'm³/h', 'l/s', 'l/min', 'l/h', 'cm³/s',
                'cm³/min', 'cm³/h', 'ft³/s', 'ft³/min', 'kft³/min',
                'ft³/h', 'Mft³/day', 'galon UK/h', 'galon US/h',
                'galon UK/min', 'galon US/min', 'galon UK/s', 'galon US/s',
                'barril/s', 'barril/min', 'barril/h', 'barril/day']
    __units__ = ['m3s', 'm3min', 'm3h', 'ls', 'lmin', 'lh', 'ccs', 'ccmin',
                 'cch', 'ft3s', 'ft3min', 'kft3min', 'ft3h', 'mft3day',
                 'galUKh', 'galUSh', 'galUKmin', 'galUSmin', 'galUKs',
                 'galUSs', 'bbls', 'bblmin', 'bblh', 'bblday']
    _magnitudes = [
        ("VolFlow", QApplication.translate("pychemqt", "Volumetric Flow")),
        ("QLiq", QApplication.translate("pychemqt", "Liquid Flow")),
        ("QGas", QApplication.translate("pychemqt", "Gas Flow"))]
    __units_set__ = {
        "VolFlow": {"altsi": "m3h", "si": "m3h", "metric": "m3s", "cgs": "ccs",
                    "english": "ft3h"},
        "QLiq": {"altsi": "m3h", "si": "m3h", "metric": "m3s", "cgs": "ccs",
                 "english": "ft3h"},
        "QGas": {"altsi": "m3h", "si": "m3h", "metric": "m3s", "cgs": "ccs",
                 "english": "ft3h"}}
    __test__ = [{"input": {"value": 1, "unit": "lmin"},
                 "prop": {"m3h": 0.06, "ft3min": 0.0353146667215,
                          "ccs": 16.6666666667}}]


class Diffusivity(unidad):
    __title__ = QApplication.translate("pychemqt", "Diffusivity")
    rates = {"m2s": 1.,
             "cm2s": k.centi**2,
             "mm2s": k.milli**2,
             "ft2s": k.foot**2,
             "inch2s": k.inch**2,
             "m2h": 1./k.hour,
             "ft2h": k.foot**2/k.hour,
             "inch2h": k.inch**2/k.hour,
             "St": k.centi**2,
             "cSt": k.milli**2}
    __text__ = ["m²/s", "cm²/s", "mm²/s", "ft²/s", "inch²/s", "m²/h",
                "ft²/h", "inch²/h", "St", "cSt"]
    __units__ = ["m2s", "cm2s", "mm2s", "ft2s", "inch2s", "m2h", "ft2h",
                 "inch2h", "St", "cSt"]
    _magnitudes = [
        ("Diffusivity", QApplication.translate("pychemqt", "Diffusivity")),
        ("KViscosity", QApplication.translate(
            "pychemqt", "Kinematic viscosity"))]
    __tooltip__ = ["m²/s", "cm²/s", "mm²/s", "ft²/s", "inch²/s", "m²/h",
                   "ft²/h", "inch²/h", "Stokes", "Centistokes"]
    __units_set__ = {
        "Diffusivity": {"altsi": "m2s", "si": "m2s", "metric": "m2s",
                        "cgs": "cm2s", "english": "ft2s"},
        "KViscosity": {"altsi": "m2s", "si": "m2s", "metric": "m2s",
                       "cgs": "cm2s", "english": "ft2s"}}
    __test__ = [{"input": {"value": 5, "unit": "St"},
                 "prop": {"m2s": 0.0005, "ft2s": 0.00538195520835}}]


class HeatFlux(unidad):
    __title__ = QApplication.translate("pychemqt", "Heat Flux")
    rates = {"Wm2": 1.,
             "kWm2": k.kilo,
             "calhm2": k.calorie/k.hour,
             "calsm2": k.calorie,
             "calscm2": k.calorie/k.centi**2,
             "kcalhm2": k.kilo*k.calorie/k.hour,
             "Btuhft2": k.Btu/k.hour/k.foot**2, "Btusft2": k.Btu/k.foot**2}
    __text__ = ["W/m²", "kW/m²", "cal/hm²", "cal/sm²", "cal/scm²",
                "kcal/hm²", "Btu/hft²", "Btu/sft²"]
    __units__ = ["Wm2", "kWm2", "calhm2", "calsm2", "calscm2", "kcalhm2",
                 "Btuhft2", "Btusft2"]
    __units_set__ = {"altsi": "Wm2", "si": "Wm2", "metric": "Wm2",
                     "cgs": "calscm2", "english": "Btuhft2"}
    __test__ = [{"input": {"value": 1, "unit": "Btuhft2"},
                 "prop": {"Wm2": 3.15459074506, "kcalhm2": 2.71427501965}}]


class ThermalConductivity(unidad):
    __title__ = QApplication.translate("pychemqt", "Thermal Conductivity")
    rates = {"WmK": 1.,
             "mWmK": 1./k.kilo,
             "kWmK": k.kilo,
             "JhmK": 1./k.hour,
             "kJhmK": k.kilo/k.hour,
             "calscmK": k.calorie/k.centi,
             "calhcmK": k.calorie/k.centi/k.hour,
             "calhmmK": k.calorie/k.milli/k.hour,
             "kcalhmK": k.calorie*k.kilo/k.hour,
             "lbfsF": k.lbf/k.Rankine,
             "lbfts3F": k.lb*k.foot/k.Rankine,
             "BtuhftF": k.Btu/k.hour/k.foot/k.Rankine}
    __text__ = ['W/m·K', 'mW/m·K', "kW/m·K", 'J/h·m·K', 'cal/s·cm·K',
                'cal/h·cm·K', 'kcal/h·m·K', 'lbf/s·F', 'lb/ft·s³·F',
                'Btu/h·ft·F']
    __units__ = ['WmK', 'mWmK', "kWmK", 'JhmK', 'calscmK', 'calhcmK',
                 'kcalhmK', 'lbfsF', 'lbfts3F', 'BtuhftF']
    __units_set__ = {"altsi": "mWmK", "si": "WmK", "metric": "WmK",
                     "cgs": "calscmK", "english": "lbfts3F"}
    __test__ = [{"input": {"value": 50, "unit": "WmK"},
                 "prop": {"BtuhftF": 28.8894658271, "kcalhmK": 43.0210325048}}]


class UA(unidad):
    __title__ = QApplication.translate("pychemqt", "UA")
    rates = {"WK": 1.,
             "kWK": k.kilo,
             "mWK": k.milli,
             "JhK": 1./k.hour,
             "kJhK": k.kilo/k.hour,
             "calhK": k.calorie/k.hour,
             "kcalhK": k.calorie*k.kilo/k.hour,
             "calsK": k.calorie,
             "kcalsK": k.calorie*k.kilo,
             "BtuhF": k.Btu/k.hour/k.foot**2/k.Rankine,
             "BtusF": k.Btu/k.foot**2/k.Rankine}
    __text__ = ['W/K', 'kW/K', 'mW/K', 'J/h·K', 'kJ/h·K', 'cal/h·K',
                'kcal/h·K', 'cal/s·K', 'kcal/s·K', 'Btu/h·F', 'Btu/s·F']
    __units__ = ['WK', 'kWK', 'mWK', 'JhK', 'kJhK', 'calhK', 'kcalhK',
                 'calsK', 'kcalsK', 'BtuhF', 'BtusF']
    __units_set__ = {"altsi": "mWK", "si": "WK", "metric": "WK",
                     "cgs": "calsK", "english": "BtuhF"}
    __test__ = [{"input": {"value": 1, "unit": "BtuhF"},
                 "prop": {"WK": 5.67826334111, "kcalhK": 4.88569503537}}]


class HeatTransfCoef(unidad):
    __title__ = QApplication.translate("pychemqt", "Heat Transfer Coefficient")
    rates = {"Wm2K": 1.,
             "kWm2K": k.kilo,
             "Jhm2K": 1./k.hour,
             "kJhm2K": k.kilo/k.hour,
             "calhm2K": k.calorie/k.hour,
             "kcalhm2K": k.calorie*k.kilo/k.hour,
             "calsm2K": k.calorie,
             "kcalsm2K": k.calorie*k.kilo,
             "calscm2K": k.lbf/k.Rankine,
             "kcalscm2K": k.calorie*k.kilo/k.centi**2,
             "Btuhft2F": k.Btu/k.hour/k.foot**2/k.Rankine,
             "Btusft2F": k.Btu/k.foot**2/k.Rankine}
    __text__ = ['W/m²·K', 'kW/m²·K', 'J/h·m²·K', 'kJ/h·m²·K', 'cal/h·m³·K',
                'kcal/h·m²·K', 'cal/s·m²·K', 'kcal/s·m²·K', 'cal/s·cm²·K',
                'kcal/s·cm²·K', 'Btu/h·ft²·F', 'Btu/s·ft²·F']
    __units__ = ['Wm2K', 'kWm2K', 'Jhm2K', 'kJhm2K', 'calhm2K', 'kcalhm2K',
                 'calsm2K', 'kcalsm2K', 'calscm2K', 'kcalscm2K', 'Btuhft2F',
                 'Btusft2F']
    __units_set__ = {"altsi": "Wm2K", "si": "Wm2K", "metric": "Wm2K",
                     "cgs": "calscm2K", "english": "Btuhft2F"}
    __test__ = [{"input": {"value": 1, "unit": "Btuhft2F"},
                 "prop": {"Wm2K": 5.67826334111, "kcalhm2K": 4.88569503537}}]


class Fouling(unidad):
    __title__ = QApplication.translate("pychemqt", "Fouling Factor")
    rates = {"m2KW": 1.,
             "m2KkW": 1./k.kilo,
             "hm2KJ": k.hour,
             "hm2KkJ": k.hour/k.kilo,
             "hm2Kcal": k.hour/k.calorie,
             "hm2Kkcal": k.hour/k.kilo/k.calorie,
             "sm2Kcal": 1./k.calorie,
             "sm2Kkcal": 1./k.calorie/k.kilo,
             "scm2Kcal": k.Rankine/k.lbf,
             "scm2Kkcal": k.centi**2/k.calorie/k.kilo,
             "hft2FBtu": k.hour*k.foot**2*k.Rankine/k.Btu,
             "sft2FBtu": k.foot**2*k.Rankine/k.Btu}
    __text__ = ['m²·K/W', 'm²·K/kW', 'h·m²·K/J', 'h·m²·K/kJ',
                'h·m³·K/cal', 'h·m²·K/kcal', 's·m²·K/cal', 's·m²·K/kcal',
                's·cm²·K/cal', 's·cm²·K/kcal', 'h·ft²·F/Btu', 's·ft²·F/Btu']
    __units__ = ['m2KW', 'm2KkW', 'hm2KJ', 'hm2KkJ', 'hm2Kcal', 'hm2Kkcal',
                 'sm2Kcal', 'sm2Kkcal', 'scm2Kcal', 'scm2Kkcal', 'hft2FBtu',
                 'sft2FBtu']
    __units_set__ = {"altsi": "m2KW", "si": "m2KW", "metric": "m2KW",
                     "cgs": "scm2Kkcal", "english": "hft2FBtu"}
    __test__ = [{"input": {"value": 1, "unit": "hft2FBtu"},
                 "prop": {"m2KW": 0.176110183682, "hm2Kkcal": 0.204679169035}}]


class Tension(unidad):
    __title__ = QApplication.translate("pychemqt", "Surface Tension")
    rates = {"Nm": 1.,
             "mNm": k.milli,
             "dyncm": k.dyn/k.centi,
             "lbfft": k.lbf/k.foot}
    __text__ = ['N/m', 'mN/m', 'dyn/cm', 'lbf/ft']
    __units__ = ['Nm', 'mNm', 'dyncm', 'lbfft']
    __units_set__ = {"altsi": "mNm", "si": "Nm", "metric": "Nm",
                     "cgs": "dyncm", "english": "lbfft"}
    __test__ = [{"input": {"value": 1, "unit": "lbfft"},
                 "prop": {"Nm": 14.5939029372, "dyncm": 14593.9029372}}]


class Viscosity(unidad):
    __title__ = QApplication.translate("pychemqt", "Viscosity")
    rates = {"Pas": 1.,
             "mPas": k.milli,
             "muPas": k.micro,
             "P": 0.1,
             "cP": k.milli,
             "microP": 0.1*k.micro,
             "dynscm2": k.milli,
             "reyn": k.g*k.pound/k.inch**2,
             "lbfts": k.lbf,
             "lbfft2": k.pound/k.foot,
             "lbfinch2": k.g*k.pound/k.inch**2,
             "lbfth": k.pound/k.foot/k.hour}
    __text__ = ['Pa·s', 'mPa·s', 'µPa·s', 'P', 'cP', 'dyn/s·cm²', 'µP',
                'reyn', 'lb/ft·s', 'lbf/ft²', 'lbf/in²', 'lb/ft·h']
    __units__ = ['Pas', 'mPas', 'muPas', 'P', 'cP', 'dynscm2', 'microP',
                 'reyn', 'lbfts', 'lbfft2', 'lbfinch2', 'lbfth']
    __tooltip__ = [QApplication.translate("pychemqt", "Pascal per second"),
                   QApplication.translate("pychemqt", "Milipascal per second"),
                   QApplication.translate(
                       "pychemqt", "Micropascal per second"),
                   "Poise", "Centipoise", 'dyn/s·cm²', 'microPoise',
                   "Reyn", 'lb/ft·s', 'lbf/ft²', 'lbf/in²', 'lb/ft·h']
    __units_set__ = {"altsi": "muPas", "si": "Pas", "metric": "Pas",
                     "cgs": "dynscm2", "english": "cP"}
    __test__ = [{"input": {"value": 1, "unit": "P"},
                 "prop": {"cP": 100, "Pas": 0.1, "lbfth": 241.90883105}}]


class SolubilityParameter(unidad):
    __title__ = QApplication.translate("pychemqt", "Solubility Parameter")
    rates = {"Jm3": 1.,
             "calcc": (k.calorie*k.mega)**0.5,
             "Btuft3": (k.Btu*k.foot**-3)**0.5}
    __text__ = ["(J/m³)^0.5", "(cal/cm³)^0.5", "(Btu/ft³)^0.5"]
    __units__ = ["Jm3", "calcc", "Btuft3"]
    __units_set__ = {"altsi": "Jm3", "si": "Jm3", "metric": "Jm3",
                     "cgs": "calcc", "english": "Btuft3"}
    __test__ = [{"input": {"value": 1, "unit": "Btuft3"},
                 "prop": {"Jm3": 193.025764622, "calcc": 0.0943668467764}}]


class PotencialElectric(unidad):
    __title__ = QApplication.translate("pychemqt", "Electric Potencial")
    rates = {"Vm": 1.,
             "kVm": k.kilo,
             "MVm": k.mega,
             "Vcm": k.centi,
             "kVcm": k.kilo/k.centi,
             "MVcm": k.mega/k.centi,
             "statVm": k.statV,
             "statVcm": k.statV/k.centi}
    __text__ = ["V/m", "kV/m", "MV/m", "V/cm", "kV/cm", "MV/cm", "statV/m",
                "statV/cm"]
    __units__ = ["Vm", "kVm", "MVm", "Vcm", "kVcm", "MVcm", "statVm",
                 "statVcm"]
    __units_set__ = {"altsi": "Vm", "si": "Vm", "metric": "Vm", "cgs": "Vcm",
                     "english": "statVm"}
    __test__ = [{"input": {"value": 3, "unit": "statVcm"},
                 "prop": {"Vm": 90000}}]


class DipoleMoment(unidad):
    __title__ = QApplication.translate("pychemqt", "Dipole Moment")
    rates = {"Cm": 1.,
             "Debye": k.debye}
    __text__ = ['C·m', 'Debye']
    __units__ = ['Cm', 'Debye']
    __tooltip__ = [QApplication.translate("pychemqt", "Coulomb per meter"),
                   "Debye"]
    __units_set__ = {"altsi": "Cm", "si": "Cm", "metric": "Cm", "cgs": "Cm",
                     "english": "Debye"}
    __test__ = [{"input": {"value": 1, "unit": "Debye"},
                 "prop": {"Cm": 3.33564095198e-30, "Debye": 1}}]


class CakeResistance(unidad):
    __title__ = QApplication.translate("pychemqt", "Cake Resistance")
    rates = {"mkg": 1.,
             "cmg": k.centi/k.kilo,
             "ftlb": k.foot/k.pound}
    __text__ = ['m/kg', 'c/gr', "ft/lb"]
    __units__ = ['mkg', 'cmg', "ftlb"]
    __units_set__ = {"altsi": "mkg", "si": "mkg", "metric": "mkg",
                     "cgs": "cmg", "english": "ftlb"}
    __test__ = [{"input": {"value": 1, "unit": "ftlb"},
                 "prop": {"mkg": 0.67196897514}}]


class PackingDP(unidad):
    __title__ = QApplication.translate("pychemqt", "Packing Pressure drop")
    rates = {"mmH2Om": 1.,
             "inH2Oft": k.inch/k.milli/k.foot}
    __text__ = ['mmH2O/m', 'inH2O/ft']
    __units__ = ['mmH2Om', 'inH2Oft']
    __units_set__ = {"altsi": "mmH2Om", "si": "mmH2Om", "metric": "mmH2Om",
                     "cgs": "mmH2Om", "english": "inH2Oft"}
    __test__ = [{"input": {"value": 1, "unit": "inH2Oft"},
                 "prop": {"mmH2Om": 83.3333333333}}]


class V2V(unidad):
    __title__ = QApplication.translate("pychemqt", "Gas-Oil ratio")
    rates = {"m3m3": 1.,
             "ft3ft3": 1.,
             "ll": 1.,
             "ft3bbl": k.foot**3/k.bbl}
    __text__ = ["m³m³", "ft/bbl"]
    __units__ = ["m3m3", "ft3bbl"]
    __tooltip__ = ['m³/m³', 'cubic foot/oil barrel']
    __units_set__ = {"altsi": "m3m3", "si": "m3m3", "metric": "m3m3",
                     "cgs": "m3m3", "english": "ft3bbl"}
    __test__ = [{"input": {"value": 1, "unit": "ft3bbl"},
                 "prop": {"m3m3": 0.178107606679}}]


class InvTemperature(unidad):
    __title__ = QApplication.translate("pychemqt", "Temperature inverse")
    rates = {"K": 1.,
             "C": 1.,
             "F": 1./k.Rankine,
             "R": 1./k.Rankine,
             "Re": 1./k.Reaumur}
    __text__ = ['1/K', '1/ºC', '1/ºR', '1/ºF', '1/ºRe']
    __units__ = ['K', 'C', 'R', 'F', 'Re']
    __tooltip__ = ['1/Kelvin', '1/Celsius', '1/Rankine', '1/Fahrenheit',
                   '1/Reaumur']
    __units_set__ = {"altsi": "C", "si": "K", "metric": "C", "cgs": "C",
                     "english": "F"}
    __test__ = [{"input": {"value": 25, "unit": "C"},
                 "prop": {"K": 25, "F": 13.888888889}}]


class InvPressure(unidad):
    __title__ = QApplication.translate("pychemqt", "Pressure inverse")
    rates = {"Pa": 1.,
             "MPa": 1./k.mega,
             "hPa": 1./k.hecto,
             "kPa": 1./k.kilo,
             "bar": 1./k.bar,
             "barg": 1./k.bar,
             "baria": 10.,
             "mbar": k.kilo/k.bar,
             "psi": 1./k.psi,
             "psig": 1./k.psi,
             "atm": 1./k.atm,
             "kgcm2": 1./k.g*k.centi**2,
             "kgcm2g": 1./k.g*k.centi**2,
             "mmH2O": 1./k.g,
             "mH2O": 1./k.g/k.kilo,
             "cmH2O": 1./k.g/10,
             "inH2O": 1./k.g/k.kilo/k.inch,
             "ftH2O": 1/k.g/k.kilo/k.foot,
             "mmHg": 1./k.torr,
             "cmHg": 1./k.torr/10,
             "inHg": 1./k.torr/k.inch/k.kilo,
             "ftHg": 1./k.torr/k.foot/k.kilo,
             "torr": 1./k.torr,
             "lbcm2": 1./k.g/k.pound*k.centi**2,
             "lbft2": k.foot**2/k.g/k.pound,
             "dyncm2": k.centi**2/k.dyn}
    __text__ = ['1/Pa', '1/hPa', '1/kPa', '1/MPa', '1/bar', '1/bar g',
                '1/mbar', '1/psi', '1/psi g', '1/atm', '1/kg/cm²',
                '1/kg/cm² g', '1/mmH2O', '1/cmH2O', '1/mH2O', '1/inH2O',
                '1/ftH2O', '1/mmHg', '1/cmHg', '1/inHg', '1/ftHg', '1/lb/cm²',
                '1/lb/ft²', '1/dyn/cm²']
    __units__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'barg', 'mbar', 'psi',
                 'psig', 'atm', 'kgcm2', 'kgcm2g', 'mmH2O', 'cmH2O', 'mH2O',
                 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg', 'lbcm2',
                 'lbft2', 'dyncm2']
    __tooltip__ = ["Pascal", "Hectopascal", "Kilopascal", "Megapascal", "bar",
                   QApplication.translate("pychemqt", "Bar gauge"),
                   "Milibar",
                   QApplication.translate("pychemqt", "Pound per square inch"),
                   QApplication.translate(
                       "pychemqt", "Pound per square inch gauge"),
                   QApplication.translate("pychemqt", "Atmosphere"),
                   QApplication.translate(
                       "pychemqt", "Atmosphere technical, kg/cm²"),
                   QApplication.translate(
                       "pychemqt", "Atmosphere technical gauge, kg/cm²g"),
                   QApplication.translate(
                       "pychemqt", "Milimeter of water column"),
                   QApplication.translate(
                       "pychemqt", "Centimeter of water column"),
                   QApplication.translate("pychemqt", "Meter of water column"),
                   QApplication.translate("pychemqt", "Inch of water column"),
                   QApplication.translate("pychemqt", "Foot of water column"),
                   QApplication.translate(
                       "pychemqt", "Milimeter of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Centimeter of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Inch of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Foot of mercury column"),
                   QApplication.translate(
                       "pychemqt", "Pound per square centimeter"),
                   QApplication.translate("pychemqt", "Pound per square foot"),
                   QApplication.translate(
                       "pychemqt", "Dyn per square centimeter")]
    __units_set__ = {"altsi": "bar", "si": "Pa", "metric": "Pa",
                     "cgs": "dyncm2", "english": "psi"}
    __test__ = [{"input": {"value": 1, "unit": "mmHg"},
                 "prop": {"bar": 750.062, "atm": 760, "psi": 51.7149}}]


class EnthalpyPressure(unidad):
    __title__ = QApplication.translate("pychemqt", "Enthalpy per pressure")
    rates = {"JkgPa": 1.,
             "kJkgkPa": 1.,
             "kJkgMPa": k.milli,
             "Jkgatm": 1./101325,
             "kJkgatm": 1./101.325,
             "Btulbpsi": k.Btu/k.lb/k.psi}
    __text__ = ['J/kgPa', 'kJ/kgkPa', 'kJ/kgMPa', "J/kgatm", "kJ/kgatm",
                "Btu/lb psi"]
    __units__ = ['JkgPa', 'kJkgkPa', 'kJkgMPa', "Jkgatm", "kJkgatm",
                 "Btulbpsi"]
    __units_set__ = {"altsi": "kJkgkPa", "si": "JkgPa", "metric": "JkgPa",
                     "cgs": "kJkgkPa", "english": "Btulbpsi"}
    __test__ = [{"input": {"value": 5, "unit": "JkgPa"},
                 "prop": {"JkgPa": 5, "kJkgMPa": 5000}}]


class EnthalpyDensity(unidad):
    __title__ = QApplication.translate("pychemqt", "Enthalpy per density")
    rates = {"Jkgkgm3": 1.,
             "kJkgkgm3": k.kilo,
             "Btulb2ft3": k.Btu/k.pound**2*k.foot**3}
    __text__ = ['J/kgkgm³', 'kJ/kgkgm³', "Btulb/lbft³"]
    __units__ = ['Jkgkgm3', 'kJkgkgm3', "Btulb2ft3"]
    __units_set__ = {"altsi": "kJkgkgm3", "si": "Jkgkgm3", "metric": "Jkgkgm3",
                     "cgs": "kJkgkgm3", "english": "Btulb2ft3"}
    __test__ = [{"input": {"value": 5, "unit": "Jkgkgm3"},
                 "prop": {"kJkgkgm3": 0.005}}]


class TemperaturePressure(unidad):
    __title__ = QApplication.translate("pychemqt", "Temperature per pressure")
    rates = {"KPa": 1.,
             "KkPa": k.milli,
             "Kbar": 1e-5,
             "KMPa": k.micro,
             "Katm": 1/101325.,
             "Fpsi": k.Rankine/k.psi}
    __text__ = ['K/Pa', 'K/kPa', "K/bar", 'K/MPa', "K/atm", "F/psi"]
    __units__ = ['KPa', 'KkPa', "Kbar", 'KMPa', "Katm", "Fpsi"]
    __units_set__ = {"altsi": "KkPa", "si": "KPa", "metric": "KPa",
                     "cgs": "KPa", "english": "Fpsi"}
    __test__ = [{"input": {"value": 1, "unit": "KPa"},
                 "prop": {"KPa": 1, "KkPa": 1000}}]


class PressureTemperature(unidad):
    __title__ = QApplication.translate("pychemqt", "Pressure per Temperature")
    rates = {"PaK": 1.,
             "kPaK": k.kilo,
             "barK": 1e5,
             "MPaK": k.mega,
             "atmK": 101325.,
             "psiF": k.psi/k.Rankine}
    __text__ = ['Pa/K', 'kPa/K', 'bar/K', 'MPa/K', "atm/K", "psi/F"]
    __units__ = ['PaK', 'kPaK', 'barK', 'MPaK', "atmK", "psiF"]
    __units_set__ = {"altsi": "kPaK", "si": "PaK", "metric": "PaK",
                     "cgs": "PaK", "english": "psiF"}
    __test__ = [{"input": {"value": 1000, "unit": "PaK"},
                 "prop": {"kPaK": 1, "atmK": 0.00986923266716}}]


class PressureDensity(unidad):
    __title__ = QApplication.translate("pychemqt", "Pressure per density")
    rates = {"Pakgm3": 1.,
             "kPakgm3": k.kilo,
             "barkgm3": 1e5,
             "Pagcc": k.liter,
             "MPakgm3": k.mega,
             "atmkgm3": 101325.,
             "psilbft3": k.psi/k.pound*k.foot**3}
    __text__ = ['Pa/kgm³', 'kPa/kgm³', 'bar/kgm³', 'MPa/kgm³', "atm/kgm³",
                "Pa/gcm³", "psi/lbft³"]
    __units__ = ['Pakgm3', 'kPakgm3', 'barkgm3', 'MPakgm3', "atmkgm3",
                 "Pagcc", "psilbft3"]
    __units_set__ = {"altsi": "kPakgm3", "si": "Pakgm3", "metric": "Pakgm3",
                     "cgs": "Pagcc", "english": "psilbft3"}
    __test__ = [{"input": {"value": 1e3, "unit": "Pakgm3"},
                 "prop": {"kPakgm3": 1, "atmkgm3": 0.00986923266716}}]


class DensityPressure(unidad):
    __title__ = QApplication.translate("pychemqt", "Density per pressure")
    rates = {"kgm3Pa": 1.,
             "kgm3kPa": k.milli,
             "kgm3bar": 1/1e5,
             "gccPa": 1./k.liter,
             "kgm3MPa": k.micro,
             "kgm3atm": 1/101325.,
             "lbft3psi": k.pound/k.foot**3/k.psi}
    __text__ = ['kg/m³Pa', 'kg/m³kPa', 'kg/m³MPa', "kg/m³bar", 'kg/m³atm',
                "gcm³/Pa", "lb/ft³psi"]
    __units__ = ['kgm3Pa', 'kgm3kPa', 'kgm3MPa', "kgm3bar", "kgm3atm", "gccPa",
                 "lbft3psi"]
    __units_set__ = {"altsi": "kgm3kPa", "si": "kgm3Pa", "metric": "kgm3kPa",
                     "cgs": "gccPa", "english": "lbft3psi"}
    __test__ = [{"input": {"value": 5, "unit": "lbft3psi"},
                 "prop": {"kgm3Pa": 0.0116164084484,
                          "kgm3atm": 1177.03258603}}]


class DensityTemperature(unidad):
    __title__ = QApplication.translate("pychemqt", "Density per temperature")
    rates = {"kgm3K": 1.,
             "gccK": 1./k.liter,
             "lbft3F": k.pound/k.foot**3/k.Rankine}
    __text__ = ['kg/m³K', 'g/cm³K', "lb/ft³F"]
    __units__ = ['kgm3K', 'gccK', "lbft3F"]
    __units_set__ = {"altsi": "kgm3K", "si": "kgm3K", "metric": "kgm3K",
                     "cgs": "gccK", "english": "lbft3F"}
    __test__ = [{"input": {"value": 1, "unit": "gccK"},
                 "prop": {"kgm3K": 1000, "lbft3F": 34.6822003201}}]


class Currency(unidad):
    """Class that models a currency rate
    Supported many currency codes from ISO 4217 using the currencyscoop web
    service

    >>> S=Currency(5, "eur")
    """
    filename = conf_dir+"moneda.dat"
    try:
        archivo = open(filename, "r")
        rates = json.load(archivo)
    except (FileNotFoundError, TypeError):
        getrates(filename)
        archivo = open(filename, "r")
        rates = json.load(archivo)
    archivo.close
    date = rates.pop("date")
    __title__ = QApplication.translate("pychemqt", "Currency")

    # Main currencies
    _uMain = [
      ("usd", QApplication.translate("pychemqt", "United States dollar"), "$"),
      ("eur", QApplication.translate("pychemqt", "Euro"), "€"),
      ("gbp", QApplication.translate("pychemqt", "Pound sterling"), "£"),
      ("jpy", QApplication.translate("pychemqt", "Japanese yen"), "¥"),
      ("cny", QApplication.translate("pychemqt", "Chinese yuan"), "¥"),
      ("rub", QApplication.translate("pychemqt", "Russian rouble"), "руб"),
      ("aud", QApplication.translate("pychemqt", "Australian dollar"), "A$"),
      ("brl", QApplication.translate("pychemqt", "Brazilian real"), "R$"),
      ("cad", QApplication.translate("pychemqt", "Canadian dollar"), "C$"),
      ("chf", QApplication.translate("pychemqt", "Swiss franc"), "Fr.")]

    # Europe
    _uEurope = [
      ("dkk", QApplication.translate("pychemqt", "Danish krone"), "kr"),
      ("isk", QApplication.translate("pychemqt", "Icelandic króna"), "Íkr"),
      ("nok", QApplication.translate("pychemqt", "Norwegian krone"), "kr"),
      ("sek", QApplication.translate("pychemqt", "Swedish krona"), "kr"),
      ("all", QApplication.translate("pychemqt", "Albanian lek"), "L"),
      ("bgn", QApplication.translate("pychemqt", "Bulgarian lev"), "лв"),
      ("czk", QApplication.translate("pychemqt", "Czech koruna"), "Kč"),
      ("huf", QApplication.translate("pychemqt", "Hungarian forint"), "Ft"),
      ("pln", QApplication.translate("pychemqt", "Polish złoty"), "zł"),
      ("ron", QApplication.translate("pychemqt", "Romanian new leu"), "RON"),
      ("bam", QApplication.translate(
          "pychemqt", "Bosnia and Herzebgovina convertible mark"), "KM"),
      ("hrk", QApplication.translate("pychemqt", "Croatian kuna"), "kn"),
      ("mkd", QApplication.translate("pychemqt", "Macedonian denar"), "ден"),
      ("mdl", QApplication.translate("pychemqt", "Moldovan leu"), "lei"),
      ("rsd", QApplication.translate("pychemqt", "Serbian dinar"), "дин."),
      ("byn", QApplication.translate("pychemqt", "Belarusian ruble"), "p."),
      ("uah", QApplication.translate("pychemqt", "Ukrainian hryvnia"), "₴"),
      ("try", QApplication.translate("pychemqt", "Turkish lira"), "TL")]

    # America
    _uAmerica = [
      ("ars", QApplication.translate("pychemqt", "Argentine peso"), "$"),
      ("bob", QApplication.translate("pychemqt", "Bolivian boliviano"), "Bs"),
      ("clf", QApplication.translate("pychemqt", "Chilean Unit of Account"),
         "UF"),
      ("clp", QApplication.translate("pychemqt", "Chilean peso"), "$"),
      ("cop", QApplication.translate("pychemqt", "Colombian peso"), "$"),
      ("crc", QApplication.translate("pychemqt", "Costa Rican colon"), "₡"),
      ("cuc", QApplication.translate("pychemqt", "Cuban convertible peso"),
       "CUC$"),
      ("cup", QApplication.translate("pychemqt", "Cuban peso"), "₱"),
      ("dop", QApplication.translate("pychemqt", "Dominican peso"), "RD$"),
      ("gtq", QApplication.translate("pychemqt", "Guatemalan quetzal"), "Q"),
      ("hnl", QApplication.translate("pychemqt", "Honduran lempira"), "L"),
      ("mxn", QApplication.translate("pychemqt", "Mexican peso"), "$"),
      ("mxv", QApplication.translate("pychemqt", "Mexican UDI"), "$"),
      ("nio", QApplication.translate("pychemqt", "Nicaraguan córdoba"), "C$"),
      ("pab", QApplication.translate("pychemqt", "Panamanian balboa"), "฿"),
      ("pyg", QApplication.translate("pychemqt", "Paraguayan guaraní"), "₲"),
      ("pen", QApplication.translate("pychemqt", "Peruvian nuevo sol"), "S/."),
      ("svc", QApplication.translate("pychemqt", "Salvadoran colón"), "₡"),
      ("uyu", QApplication.translate("pychemqt", "Uruguayan peso"), "$U"),
      # ("vef", QApplication.translate("pychemqt", "Venezuelan bolívar"), "Bs"),
      ("ved", QApplication.translate("pychemqt", "Venezuelan digital bolívar"),
          "Bs.D"),
      ("ves", QApplication.translate(
          "pychemqt", "Venezuelan sovereign bolívar"), "Bs.S"),

      ("awg", QApplication.translate("pychemqt", "Aruban florin"), "Afl."),
      ("bsd", QApplication.translate("pychemqt", "Bahamian dollar"), "B$"),
      ("bbd", QApplication.translate("pychemqt", "Barbados dollar"), "Bds$"),
      ("bzd", QApplication.translate("pychemqt", "Belize dollar"), "BZ$"),
      ("bmd", QApplication.translate("pychemqt", "Bermudean dollar"), "BD$"),
      ("kyd", QApplication.translate("pychemqt", "Cayman Islands dollar"),
          "CI$"),
      ("xcd", QApplication.translate("pychemqt", "East Caribbean dollar"),
          "EC$"),
      ("gyd", QApplication.translate("pychemqt", "Guyanese dollar"), "GY$"),
      ("htg", QApplication.translate("pychemqt", "Haitian gourde"), "G"),
      ("jmd", QApplication.translate("pychemqt", "Jamaican dollar"), "J$"),
      ("ang", QApplication.translate(
          "pychemqt", "Netherlands Antillean guilder"), "f"),
      ("srd", QApplication.translate("pychemqt", "Surinamese dollar"), "$"),
      ("ttd", QApplication.translate("pychemqt", "Trinidad and Tobago dollar"),
          "TT$")]

    # Asia
    _uAsia = [
      ("afn", QApplication.translate("pychemqt", "Afghan afghani"), "Af"),
      ("bhd", QApplication.translate("pychemqt", "Bahraini dinar"), "BD"),
      ("bnd", QApplication.translate("pychemqt", "Brunei dollar"), "B$"),
      ("ils", QApplication.translate("pychemqt", "Israeli new shekel"), "₪"),
      ("iqd", QApplication.translate("pychemqt", "Iraqi dinar"), "د.ع"),
      ("irr", QApplication.translate("pychemqt", "Iranian rial"), "﷼"),
      ("jod", QApplication.translate("pychemqt", "Jordanian dinar"), "JD"),
      ("kwd", QApplication.translate("pychemqt", "Kuwaiti dinar"), "د.ك"),
      ("lbp", QApplication.translate("pychemqt", "Lebanese pound"), "ل.ل."),
      ("omr", QApplication.translate("pychemqt", "Omani rial"), "﷼"),
      ("pkr", QApplication.translate("pychemqt", "Pakistani rupee"), "Rs"),
      ("qar", QApplication.translate("pychemqt", "Qatari riyal"), "QR"),
      ("sar", QApplication.translate("pychemqt", "Saudi riyal"), "ر.س"),
      ("syp", QApplication.translate("pychemqt", "Syrian pound"), "£S"),
      ("aed", QApplication.translate(
          "pychemqt", "United Arab Emirates dirham"), "إ.د"),
      ("yer", QApplication.translate("pychemqt", "Yemeni rial"), "﷼"),

      ("amd", QApplication.translate("pychemqt", "Armenian dram"), "֏"),
      ("azn", QApplication.translate("pychemqt", "Azerbaijan manat"), "₼"),
      ("gel", QApplication.translate("pychemqt", "Georgian lari"), "ლ"),
      ("kzt", QApplication.translate("pychemqt", "Kazakhstani tenge"), "₸"),
      ("kgs", QApplication.translate("pychemqt", "Kyrgyzstani som"), "som"),
      ("tjs", QApplication.translate("pychemqt", "Tajikistani somoni"), "som"),
      ("tmt", QApplication.translate("pychemqt", "Turkmenistan manat"), "T"),
      ("uzs", QApplication.translate("pychemqt", "Uzbekistan som"), "som"),

      ("bdt", QApplication.translate("pychemqt", "Bangladeshi taka"), "৳"),
      ("btn", QApplication.translate("pychemqt", "Bhutanese ngultrum"), "Nu."),
      ("cnh", QApplication.translate("pychemqt", "Renminbi"), "¥"),
      ("khr", QApplication.translate("pychemqt", "Cambodian riel"), "៛"),
      ("kpw", QApplication.translate("pychemqt", "North Korean won"), "₩"),
      ("hkd", QApplication.translate("pychemqt", "Hong Kong dollar"), "HK$"),
      ("inr", QApplication.translate("pychemqt", "Indian rupee"), "₨"),
      ("idr", QApplication.translate("pychemqt", "Indonesian rupiah"), "Rp"),
      ("lak", QApplication.translate("pychemqt", "Lao kip"), "₭"),
      ("mop", QApplication.translate("pychemqt", "Macanese pataca"), "MOP$"),
      ("myr", QApplication.translate("pychemqt", "Malaysian ringgit"), "RM"),
      ("mnt", QApplication.translate("pychemqt", "Mongolian tögrög"), "₮"),
      ("mmk", QApplication.translate("pychemqt", "Myanmar kyat"), "K"),
      ("npr", QApplication.translate("pychemqt", "Nepalese rupee"), "रु"),
      ("twd", QApplication.translate("pychemqt", "New Taiwan dollar"), "NT$"),
      ("php", QApplication.translate("pychemqt", "Philippine peso"), "PhP"),
      ("sgd", QApplication.translate("pychemqt", "Singapore dollar"), "S$"),
      ("krw", QApplication.translate("pychemqt", "South Korean won"), "₩"),
      ("lkr", QApplication.translate("pychemqt", "Sri Lankan rupee"), "₨"),
      ("thb", QApplication.translate("pychemqt", "Thai baht"), "฿"),
      ("vnd", QApplication.translate("pychemqt", "Vietnamese dong"), "₫")]

    # Africa
    _uAfrica = [
      ("dzd", QApplication.translate("pychemqt", "Algerian dinar"), "دج"),
      ("aoa", QApplication.translate("pychemqt", "Angolan kwanza"), "Kz"),
      ("bwp", QApplication.translate("pychemqt", "Botswana pula"), "P"),
      ("bif", QApplication.translate("pychemqt", "Burundian franc"), "FBu"),
      ("cve", QApplication.translate("pychemqt", "Cape Verde escudo"), "$"),
      ("kmf", QApplication.translate("pychemqt", "Comoro franc"), "CF"),
      ("cdf", QApplication.translate("pychemqt", "Congolese franc"), "FC"),
      ("djf", QApplication.translate("pychemqt", "Djiboutian franc"), "Fdj"),
      ("egp", QApplication.translate("pychemqt", "Egyptian pound"), "E£"),
      ("ern", QApplication.translate("pychemqt", "Eritrean nakfa"), "Nfk"),
      ("etb", QApplication.translate("pychemqt", "Ethiopian birr"), "Br"),
      ("gmd", QApplication.translate("pychemqt", "Ghambian dalasi"), "D"),
      ("ghs", QApplication.translate("pychemqt", "Ghanaian cedi"), "GH₵"),
      ("gnf", QApplication.translate("pychemqt", "Guinean franc"), "GFr"),
      ("kes", QApplication.translate("pychemqt", "Kenyan shilling"), "KSh"),
      ("lsl", QApplication.translate("pychemqt", "Lesotho loti"), "L"),
      ("lrd", QApplication.translate("pychemqt", "Liberian dollar"), "L$"),
      ("lyd", QApplication.translate("pychemqt", "Libyan dinar"), "ل.د"),
      ("mru", QApplication.translate("pychemqt", "Mauritanian ouguiya"), "UM"),
      ("mur", QApplication.translate("pychemqt", "Mauritian rupee"), "₨"),
      ("mga", QApplication.translate("pychemqt", "Malagasy ariary"), "Ar"),
      ("mwk", QApplication.translate("pychemqt", "Malawian kwacha"), "MK"),
      ("mvr", QApplication.translate("pychemqt", "Maldivian rufiyaa"), "MRf"),
      ("mad", QApplication.translate("pychemqt", "Moroccan dirham"), "درهم"),
      ("mzn", QApplication.translate("pychemqt", "Mozambican metical"), "MT"),
      ("nad", QApplication.translate("pychemqt", "Namibian dollar"), "N$"),
      ("ngn", QApplication.translate("pychemqt", "Nigerian naira"), "₦"),
      ("rwf", QApplication.translate("pychemqt", "Rwandan franc"), "FRw"),
      ("stn", QApplication.translate(
          "pychemqt", "São Tomé and Príncipe dobra"), "Db"),
      ("scr", QApplication.translate("pychemqt", "Seychelles rupee"), "SR"),
      ("sll", QApplication.translate(
          "pychemqt", "Sierra Leonean leone"), "Le"),
      ("sos", QApplication.translate("pychemqt", "Somali shilling"), "Sh.So."),
      ("zar", QApplication.translate("pychemqt", "South African rand"), "R"),
      ("sdg", QApplication.translate("pychemqt", "Sudanese pound"), "ج.س"),
      ("szl", QApplication.translate("pychemqt", "Swazi lilangeni"), "L"),
      ("tzs", QApplication.translate("pychemqt", "Tanzanian shilling"), "TSh"),
      ("tnd", QApplication.translate("pychemqt", "Tunisian dinar"), "د.ت"),
      ("ugx", QApplication.translate("pychemqt", "Ugandan shilling"), "USh"),
      ("zmw", QApplication.translate("pychemqt", "Zambian kwacha"), "ZK"),
      ("zwl", QApplication.translate("pychemqt", "Zimbabwean dollar"), "Z$"),
      ("xaf", QApplication.translate(
          "pychemqt", "Central AFrican CFA franc"), "FCFA"),
      ("xof", QApplication.translate(
          "pychemqt", "West African CFA franc"), "CFA")]

    # Oceania
    _uOceania = [
      ("fjd", QApplication.translate("pychemqt", "Fiji dollar"), "FJ$"),
      ("nzd", QApplication.translate("pychemqt", "New Zealand dollar"), "NZ$"),
      ("pgk", QApplication.translate("pychemqt", "Papua New Guinean kina"),
          "K"),
      ("sbd", QApplication.translate("pychemqt", "Salomon Islands dollar"),
          "SI$"),
      ("wst", QApplication.translate("pychemqt", "Samoan tala"), "WS$"),
      ("top", QApplication.translate("pychemqt", "Tongan pa'anga"), "T$"),
      ("tvd", QApplication.translate("pychemqt", "Tuvalu dollar"), "TV$"),
      ("vuv", QApplication.translate("pychemqt", "Vanuatu vatu"), "VT"),
      ("xpf", QApplication.translate("pychemqt", "CFP franc"), "F")]

    # Crypto
    _uCrypto = [
        ("ada", "Cardano", "₳"),
        ("bch", "Bitcoin Cash", ""),
        ("xbt", "Bitcoin", "₿"),
        ("ltc", "Litecoin", "Ł"),
        ("doge", "Dogecoin", "Ð"),
        ("xrp", "Ripple", ""),
        ("xlm", "Stellar", ""),
        ("eth", "Ethereum", "Ξ"),
        ("dot", "Polkadot", ""),
        ("uni", "Uniswap", ""),
        ("link", "Chainlink", "")]


    # Unused
    _uUnused = [
      "fkp", # Falkland Islands pound, same to sterling pound
      "ggp", # Guernsey pound, same to sterling pound
      "gip", # Gibraltar pound, same to sterling pound
      "imp", # Manx pound, same to sterling pound
      "jep", # Jersey pound, same to sterling pound
      "shp", # Saint Helena pound, same to sterling pound
      "ats", # Old austrian schilling
      "azm", # Old azerbaijani manat
      "bef", # Old belgian franc
      "btc", # Synonims for bitcoin
      "byr", # Old belarusian ruble
      "cyp", # Old cypriot pound
      "dem", # Old german mark
      "eek", # Old estonian kroon
      "esp", # Old spanish peseta
      "fim", # Old finnish markka
      "frf", # Old french franc
      "ghc", # Old ghanaian cedi
      "grd", # Old greek drachma
      "iep", # Old irish pound
      "itl", # Old italian lira
      "ltl", # Old Lithuanian litas
      "luf", # Old Luxembourg franc
      "lvl", # Old latvian lats
      "mgf", # Old malagasy franc
      "mtl", # Old maltese lira
      "mzm", # Old mozambican metical
      "mro", # Old mauritanian ouguiya
      "nlg", # Old dutch guilder
      "pte", # Old portuguese escudo
      "rol", # Old romanian leu
      "sdd", # Old Sudanese dinar
      "sit", # Old Slovenian tolar
      "skk", # Old Czechoslovak koruna
      "sle", # Old Sierra Leonean leone
      "spl", # Seborgan Luigino
      "srg", # Old Surinamese guilder
      "std", # Old São Tomé and Príncipe dobra
      "tmm", # Old turkmenistani manat
      "trl", # Old Turkish lira
      "val", # Old vatican lira
      "veb", # Old venezuelan bolivar
      "vef", # Old venezuelan bolivar fuerte
      "xag", # Silver (one troy ounce)
      "xau", # Gold (one troy ounce)
      "xdr", # Special drawing rights (from FMI)
      "xpd", # Palladium (one troy ounce)
      "xpt", # Platinum (one troy ounce)
      "zmk", # Old Zambian kwacha
      "zwd" # Old Zimbabwean dollar
      ]

    _uTotal = _uMain + _uEurope + _uAmerica + _uAfrica + _uAsia + _uOceania \
        + _uCrypto
    __text__ = []
    __units__ = []
    __tooltip__ = []
    for unit, tip, txt in _uTotal:
        __units__.append(unit)
        __tooltip__.append(tip)
        __text__.append(txt)

    __units_set__ = {"altsi": "usd", "si": "usd", "metric": "usd",
                     "cgs": "usd", "english": "usd"}

    @property
    def str(self):
        if self.code:
            return self.code

        conf = self.func(self.magnitud)
        num = self.format(conf, self.magnitud)
        txt = self.text(self.magnitud)
        return " "+txt+num


if os.environ["icu"] == "True":
    import icu
    locale = QtCore.QLocale.system().name()

    subclasses = unidad.__subclasses__()
    names = [unit.__title__ for unit in subclasses]
    collator = icu.Collator.createInstance(icu.Locale(locale))
    sortfunc = collator.getSortKey
    title_sorted = sorted(names, key=sortfunc)
    _all = [0]*len(names)
    for _unit in subclasses:
        i = title_sorted.index(_unit.__title__)
        _all[i] = _unit
else:
    _all = sorted(unidad.__subclasses__(), key=lambda item: item.__title__)


# Auto documenting unidad subclasses
for _clas in _all:
    if _clas.__doc__:
        continue
    doc = QApplication.translate(
        "pychemqt", "Class to model a %s measure" % _clas.__title__)
    doc += os.linesep + os.linesep
    doc += QApplication.translate("pychemqt", "Supported units") + "::"
    doc += os.linesep
    default = True
    for i, key in enumerate(_clas.__units__):
        # Add list of supported unit with name and symbol
        if _clas.__tooltip__:
            name = _clas.__tooltip__[i]
        else:
            name = _clas.__text__[i]
        doc += "    * %s (%s)" % (name, key)

        # Mark the default unit
        if default:
            doc += " (%s)" % QApplication.translate("pychemqt", "default")
            default = False
        doc += os.linesep + os.linesep

    # Add doctest example
    doc += "Examples" + os.linesep
    doc += "--------" + os.linesep + os.linesep
    title = _clas.__name__
    for test in _clas.__test__:
        doc += ">>> %s = %s(%g, '%s')" % (
            title[0], title, test["input"]["value"], test["input"]["unit"])
        doc += os.linesep + ">>> "
        template = []
        values = []
        for key in test["prop"].keys():
            template.append("%g")
            values.append("%s.%s" % (title[0], key))
        doc += '"%s"' % " ".join(template)
        doc += ' % (' + ", ".join(values)
        doc += ")" + os.linesep
        values = []
        for value in test["prop"].values():
            values.append("%g" % value)
        doc += "'" + " ".join(values) + "'" + os.linesep
    doc += os.linesep + os.linesep
    _clas.__doc__ = doc


_magnitudes = []
for unit in _all:
    for magnitud in unit.magnitudes():
        _magnitudes.append(magnitud+(unit, ))
_magnitudes.append(("Dimensionless",
                    QApplication.translate("pychemqt", "Dimensionless"),
                    Dimensionless))

unit_set = {}
for unit in _all:
    if unit._magnitudes:
        unit_set.update(unit.__units_set__)
    else:
        unit_set[unit.__name__] = unit.__units_set__

units_set = {}
for _set in ("altsi", "si", "metric", "cgs", "english"):
    units_set[_set] = []
    for magnitud, titulo, unit in _magnitudes[:-1]:
        units_set[_set].append(unit.__units__.index(unit_set[magnitud][_set]))
