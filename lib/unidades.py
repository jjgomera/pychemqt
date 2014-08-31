#!/usr/bin/python
# -*- coding: utf-8 -*-

import cPickle
from ConfigParser import ConfigParser

from PyQt4.QtGui import QApplication
import scipy.constants as k

from config import conf_dir, getMainWindowConfig, representacion
from firstrun import getrates

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
    Generic class to model units
    Each child class must define the following parameters:
        __title__: Title or name of class
        rates: Dict with conversion rates
        __text__: List with units title
        __units__: List with units properties names
        __tooltip__: List with help string for units
        _magnitudes: Opcional to units with several magnituds
            Each magnitud is a tuple with format (Name, title)
        __units_set__: Dict with standart unit for units system,
            altsi, si, metric, cgs, english
    """
    __title__ = ""
    rates = {}
    __text__ = []
    __units__ = []
    __tooltip__ = []
    _magnitudes = []
    __units_set__ = []

    def __init__(self, data, unit="", magnitud=""):
        """Non proportional magnitudes (Temperature, Pressure)
        must rewrite this method"""
        if not magnitud:
            magnitud = self.__class__.__name__
        self.magnitud = magnitud

        if data is None:
            self._data = 0
            self.code = "n/a"
        else:
            self._data = data
            self.code = ""

        if unit == "conf":
            unit = self.__units__[self.Config.getint('Units', magnitud)]
        try:
            conversion = self.__class__.rates[unit]
        except:
            raise ValueError(
                QApplication.translate("pychemqt", "Wrong input code"))

        self._data *= conversion
        for key in self.__class__.rates:
            self.__setattr__(key, self._data / self.__class__.rates[key])

    def __new__(cls, data, unit="", magnitud=""):
        """Non proportional magnitudes (Temperature, Pressure)
        must rewrite this method"""
        cls.Config = getMainWindowConfig()

        if data:
            data = float(data)
        else:
            data = 0.

        if not magnitud:
            magnitud = cls.__name__

        if not unit:
            unit = cls.__units__[0]
        elif unit == "conf":
            unit = cls.__units__[cls.Config.getint('Units', magnitud)]

        return float.__new__(cls, data * cls.rates[unit])

    def config(self, magnitud=""):
        """Using config file return the value in the configurated unit"""
        if not magnitud:
            magnitud = self.__class__.__name__
        value = self.Config.getint('Units', magnitud)
        return self.__getattribute__(self.__units__[value])

    @classmethod
    def text(cls, magnitud=""):
        """Using config file return the configurated unit text"""
        if not magnitud:
            magnitud = cls.__name__
        return cls.__text__[cls.Config.getint("Units", magnitud)]

    @classmethod
    def func(cls, magnitud=""):
        """Return the configurated unit name for getattribute call"""
        if not magnitud:
            magnitud = cls.__name__
        return cls.__units__[cls.Config.getint("Units", magnitud)]

    @classmethod
    def magnitudes(cls):
        """Return the magnitudes list for unit,
if a unit define several magnitudes, must be fill the _magnitudes variable"""
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
                txt=self.__text__[self.__units__.index(conf)]
            num = self.format(conf, self.magnitud)
            return num+" "+txt
    str = property(get_str)


class Dimensionless(float):
    """Dummy class to integrate dimensionless magnitudes
with support for class unidad operations: txt, config. func."""
    __title__ = QApplication.translate("pychemqt", "Dimensionless")
    __text__ = [""]
    _magnitudes = []

    def __new__(cls, *args, **kwargs):
        """Discard superfluous parameters for this class"""
        if args[0] is None:
            val = 0
            cls.code = "n/a"
        else:
            val = args[0]
            cls.code = ""

        return float.__new__(cls, val)

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
        return num


class Temperature(unidad):
    """Class that models a temperature measure
    Supported units:

    * Kelvin (K) default
    * Celsius (C)
    * Fahrenheit (F)
    * Rankine (R)
    * Reaumur (Re)

    >>> T=Temperature(25, "C")
    >>> print T.K, T.C, T.F
    298.15 25.0 77.0
    """
    __title__ = QApplication.translate("pychemqt", "Temperature")
    __text__ = ['K', u'ºC', u'ºR', u'ºF', u'ºRe']
    __units__ = ['K', 'C', 'R', 'F', 'Re']
    __tooltip__ = ['Kelvin', 'Celsius', 'Rankine', 'Fahrenheit', 'Reaumur']
    __units_set__ = {"altsi": "C", "si": "K", "metric": "C", "cgs": "C",
                     "english": "F"}

    def __init__(self, data, unit="K", magnitud=""):

        if not magnitud:
            magnitud = self.__class__.__name__
        self.magnitud = magnitud

        if data is None:
            data = 0
            self.code = "n/a"
        else:
            data = data
            self.code = ""

        if unit == "conf":
            unit = self.__units__[self.Config.getint('Units', magnitud)]

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

    def __new__(cls, data, unit="K", magnitud=""):
        cls.Config = getMainWindowConfig()

        if data is None:
            cls.code = "n/a"
            data = 0.
        else:
            cls.code = ""
            data = float(data)

        if not magnitud:
            magnitud = cls.__name__

        if unit == "conf":
            unit = cls.__units__[cls.Config.getint('Units', magnitud)]

        if unit == "K":
            cls._data = data
        elif unit == "C":
            cls._data = C2K(data)
        elif unit == "F":
            cls._data = F2K(data)
        elif unit == "R":
            cls._data = R2K(data)
        elif unit == "Re":
            cls._data = Re2K(data)
        else:
            raise ValueError(
                QApplication.translate("pychemqt", "Wrong input code"))

        return float.__new__(cls, cls._data)


class DeltaT(unidad):
    """Class that models a delta temperature measure
    Supported units:

    * Kelvin (default)
    * Celsius
    * Fahrenheit
    * Rankine
    * Reaumur

    >>> T=DeltaT(25, "C")
    >>> print T.K, T.F
    25.0 45.0
    """
    __title__ = QApplication.translate("pychemqt", "Temperature increase")
    rates = {"K": 1.,
             "C": 1.,
             "F": k.Rankine,
             "R": k.Rankine,
             "Re": k.Reaumur}
    __text__ = ['K', u'ºC', u'ºR', u'ºF', u'ºRe']
    __units__ = ['K', 'C', 'R', 'F', 'Re']
    __tooltip__ = ['Kelvin', 'Celsius', 'Rankine', 'Fahrenheit', 'Reaumur']
    __units_set__ = {"altsi": "C", "si": "K", "metric": "C", "cgs": "C",
                     "english": "F"}

    def __init__(self, data, unit="K", magnitud=""):
        super(DeltaT, self).__init__(data, unit, magnitud)


class Angle(unidad):
    """Class that models a angle measure
    Supported units:

    * radian (rad) default
    * grade (deg)
    * minute (min)
    * second (sec)
    * gradian (grad)

    >>> angle=Angle(25, "deg")
    >>> print angle.rad
    0.436332312999
    """
    __title__ = QApplication.translate("pychemqt", "Angle")
    rates = {"rad": 1.,
             "deg": 2*k.pi/360,
             "min": 2*k.pi/360/60,
             "sec": 2*k.pi/360/3600,
             "grad": 2*k.pi/400}
    __text__ = ["rad", u"º deg", "'", '"', "grad"]
    __units__ = ["rad", "deg", "min", "sec", "grad"]
    __tooltip__ = [QApplication.translate("pychemqt", "Radian"),
                   QApplication.translate("pychemqt", "Degree"),
                   QApplication.translate("pychemqt", "Arcminute"),
                   QApplication.translate("pychemqt", "Arcsecond"),
                   QApplication.translate("pychemqt", "Gradian")]
    __units_set__ = {"altsi": "rad", "si": "rad", "metric": "rad",
                     "cgs": "rad", "english": "rad"}

    def __init__(self, data, unit="rad", magnitud=""):
        super(Angle, self).__init__(data, unit, magnitud)


class Length(unidad):
    """Class that models a length measure
    Supported units:

    * meter (m) default
    * millimeter (mm)
    * centimeter (cm)
    * micrometer (micra)
    * kilometer (km)
    * inch (inch)
    * foot (ft)
    * yard (yd)
    * milla (milla)
    * milla nautica (milla_nau)
    * amstrong (A)

    >>> L=Length(12, "inch")
    >>> print L.m, L.inch, L.ft
    0.3048 12.0 1.0
    """
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
    __text__ = ['m', 'cm', 'mm', u'µm', 'km', 'inch', 'ft', 'yard', 'milla',
                "M", "pm", u"Å"]
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
                   QApplication.translate("pychemqt", "icometer"), u"Ångström"]
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

    def __init__(self, data, unit="m", magnitud=""):
        super(Length, self).__init__(data, unit, magnitud)


class Area(unidad):
    """Class that models a area measure
    Supported units:

    * square meter (m2) default
    * square centimeter (cm2)
    * square milimeter (mm2)
    * square foot (ft2)
    * square inch (in2)
    * square yard (yd2)
    * hectare (ha)
    * acre (acre)

    >>> S=Area(1, "ft2")
    >>> print S.m2, S.inch2
    0.09290304 144.0
    """
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
    __text__ = [u'm²', u'cm²', u'mm²', u"km2", u'inch²', u'ft²', u'yd²',
                'ha', u'acre']
    __units__ = ['m2', 'cm2', 'mm2', 'km2', 'inch2', 'ft2', 'yd2', 'ha', "acre"]
    __units_set__ = {"altsi": "m2", "si": "m2", "metric": "m2", "cgs": "cm2",
                     "english": "ft2"}

    def __init__(self, data, unit="m2", magnitud=""):
        super(Area, self).__init__(data, unit, magnitud)


class Volume(unidad):
    """Class that models a volume measure
    Supported units:

    * cubic meter (m3) default
    * cubic centimeter (cc)
    * liter (l)
    * mililiter (ml)
    * cubic yard (yd3)
    * cubic foot (ft3)
    * cubic inch (inch3)
    * US gallon (galUS)
    * british gallon (galUK)
    * US liquid quart (qtUSliq)
    * US dry quart (qtUSdry)
    * british quart (qtUK)
    * barrel of petrolium (bbl)
    * ounce (onz)
    * british ounce (onzUK)

    >>> V=Volume(1, "bbl")
    >>> print V.l, V.ft3, V.galUS
    158.987294928 5.61458333333 42.0
    """
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
    __text__ = [u'm³', u'cm³', 'liter', u'yd³', u'ft³', u'inch³', 'galon US',
                'galon UK', 'oil', 'bbl', 'bblUK', 'onza', 'onza UK']
    __units__ = ['m3', 'cc', 'l', 'yd3', 'ft3', 'inch3', 'galUS', 'galUK',
                 'bbl', 'bblUS', "bblUK", 'onz', 'onzUK']
    __tooltip__ = [u'm³', u'cm³',
                   QApplication.translate("pychemqt", "liter"),
                   u'yd³', u'ft³', u'inch³',
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

    def __init__(self, data, unit="m3", magnitud=""):
        super(Volume, self).__init__(data, unit, magnitud)


class Time(unidad):
    """Class that models a time measure
    Supported units:

    * second (s) default
    * minute (min)
    * hour (h)
    * day (day)
    * year (year)

    >>> t=Time(1, "day")
    >>> print t.min, t.h
    1440.0 24.0
    """
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

    def __init__(self, data, unit="s", magnitud=""):
        super(Time, self).__init__(data, unit, magnitud)


class Frequency(unidad):
    """Class that models a frequency measure
    Supported units:

    * rpm default
    * rps
    * rph
    * Hz
    * rad/s (rads)
    * rad/min (radmin)
    * rad/hour (radhr)

    >>> t=Frequency(1, "rads")
    >>> print t.rpm
    9.54929658551
    """
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

    def __init__(self, data, unit="rpm", magnitud=""):
        super(Frequency, self).__init__(data, unit, magnitud)


class Speed(unidad):
    """Class that models a speed measure
    Supported units:

    * meter per second (ms) default
    * centimeter per second (cms)
    * milimeter per second (mms)
    * kilometer per second (kms)
    * meter per minute (mmin)
    * kilometer per minute (kmmin)
    * kilometer per hour (kmh)
    * meter per day (mday)
    * kilometer per day (kmday)
    * foot per second (fts)
    * foot per minute (ftmin)
    * foot per hour (fth)
    * foot per day (ftday)
    * inch per second (inchs)
    * inch per minute (inchmin)
    * mille per hour (mph)
    * knot (kt)

    >>> V=Speed(1, "ms")
    >>> print V.mmin, V.kmh, V.fts
    60.0 3.6 3.28083989501
    """
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
                'km/min', 'km/h',  'km/day', 'mph', 'nudo']
    __units__ = ['ms', 'cms', 'mms', 'kms', 'fts', 'ftmin',  'mmin', 'kmmin',
                 'kmh', 'kmday', 'mph', 'kt']
    __tooltip__ = ['m/s', 'cm/s', 'mm/s', 'km/s', 'ft/s', 'ft/min', 'm/min',
                   'km/min', 'km/h',  'km/day', 'mph',
                   QApplication.translate("pychemqt", "Knot")]
    __units_set__ = {"altsi": "ms", "si": "ms", "metric": "ms", "cgs": "cms",
                     "english": "fts"}

    def __init__(self, data, unit="ms", magnitud=""):
        super(Speed, self).__init__(data, unit, magnitud)


class Acceleration(unidad):
    """Class that models a acceleration measure
    Supported units:

    * meter per square second (ms2) default
    * centimeter per square second (cms2)
    * foot per square second (fts2)
    * inch per square second (inchs2)

    >>> g=Acceleration(9.81)
    >>> print g.fts2
    32.1850393701
    """
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
    __text__ = [u"m/s²", u"cm/s²", u"ft/s²", u"inch/s²", u"yd/s²", u"m/min²",
                u"cm/min²", u"ft/min²", u"inch/min²"]
    __units__ = ["ms2", "cms2", "fts2", "inchs2", "yds2", "mmin2", "cmmin2",
                 "ftmin2", "inchmin2"]
    __units_set__ = {"altsi": "ms2", "si": "ms2", "metric": "ms2",
                     "cgs": "cms2", "english": "fts2"}

    def __init__(self, data, unit="ms2", magnitud=""):
        super(Acceleration, self).__init__(data, unit, magnitud)


class Mass(unidad):
    """Class that models a mass measure
    Supported units:

    * kilogram (kg) default
    * gram (g)
    * miligram (mg)
    * tonne (Ton)
    * pound (lb)
    * grain (gr)
    * ounce avoirdupois (oz)
    * slug (slug)
    * hundredweight US (cwtUS)
    * hundredweight UK (cwtUK)
    * short tonne (TonUK)
    * long tonne (TonUS)

    >>> M=Mass(1, "lb")
    >>> print M.kg, M.g, M.oz
    0.45359237 453.59237 16.0
    """
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

    def __init__(self, data, unit="kg", magnitud=""):
        super(Mass, self).__init__(data, unit, magnitud)


class Mol(unidad):
    """Class that models a mol measure
    Supported units:

    * kilomol (kmol) default
    * mol (mol)
    * milimol (mmol)

    >>> M=Mol(1, "kmol")
    >>> print M.mol, M.lbmol
    1000.0 2.20462262185
    """
    __title__ = QApplication.translate("pychemqt", "Mol")
    rates = {"kmol": 1.,
             "mol": 1./k.kilo,
             "milimol": 1./k.mega,
             "lbmol": k.pound}
    __text__ = ['kmol', 'mol', 'mmol', "lbmol"]
    __units__ = ['kmol', 'mol', 'mmol', "lbmol"]
    __units_set__ = {"altsi": "kmol", "si": "kmol", "metric": "kmol",
                     "cgs": "mol", "english": "lbmol"}

    def __init__(self, data, unit="kmol", magnitud=""):
        super(Mol, self).__init__(data, unit, magnitud)


class SpecificVolume(unidad):
    """Class that models a specific volume measure
    Supported units:

    * cubic meter per kilogram (m3kg) default
    * liter per kilogram (lkg) (same as cc/g, and  ml/g)
    * cubic meter per gram (m3g)
    * cubic centimeter per kilogram (cckg)
    * cubic foot per pound (ft3lb)
    * cubic inch  per pound (in3lb)
    * british gallon per pound (galUKlb)
    * gallon per pound(galUSlb)
    * barrel per pound (bbllb)
    * cubic foot per slug (ft3slug)
    * cubic foot per ounce (ft3oz)
    * cubic inch per ounce (inch3oz)
    * british gallon  per ounce (galUKoz)
    * gallon per ounce (galUSoz)

    >>> R=SpecificVolume(50, "lkg")
    >>> print  R.m3kg, R.ft3lb
    0.05 0.800923168698
    """
    __title__ = QApplication.translate("pychemqt", "Specific Volume")
    rates = {"m3kg": 1.,
             "lg": 1.,
             "lkg": k.liter,
             "ccg": k.liter,
             "mlg": k.liter,
             "m3g": k.gram,
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
    __text__ = [u'm³/kg', u'cm³/g', u'ml/g', u'm³/g', u'cm³/kg', u'ft³/lb',
                u'in³/lb', 'gallon UK/lb', 'gallon US/lb', 'barril/lb',
                u'ft³/ton UK', u'ft³/ton US', u'ft³/slug', u'ft³/onza',
                u'in³/onza', 'gallon UK/onza', 'gallon US/onza']
    __units__ = ['m3kg', 'lkg', 'ccg', 'mlg', 'm3g', 'cckg', 'ft3lb',
                 'inch3lb', 'galUKlb', 'galUSlb', 'bbllb', 'ft3tonUK',
                 'ft3tonUS', 'ft3slug',  'ft3oz', 'in3oz', 'galUKoz', 'galUSoz']
    __units_set__ = {"altsi": "m3kg", "si": "m3kg", "metric": "m3kg",
                     "cgs": "ccg", "english": "ft3lb"}

    def __init__(self, data, unit="m3kg", magnitud=""):
        super(SpecificVolume, self).__init__(data, unit, magnitud)


class SpecificVolume_square(unidad):
    """Class that models a specific volume squarre measure (useful too for
    third virial coefficient
    Supported units:

    * cubic meter per kilogram (m3kg) default
    * liter per kilogram (lkg) (same as cc/g, and  ml/g)
    * cubic meter per gram (m3g)
    * cubic centimeter per kilogram (cckg)
    * cubic foot per pound (ft3lb)
    * cubic inch  per pound (in3lb)
    * british gallon per pound (galUKlb)
    * gallon per pound(galUSlb)
    * barrel per pound (bbllb)
    * cubic foot per slug (ft3slug)
    * cubic foot per ounce (ft3oz)
    * cubic inch per ounce (inch3oz)
    * british gallon  per ounce (galUKoz)
    * gallon per ounce (galUSoz)

    >>> R=SpecificVolume_square(50, "lkg")
    >>> print  R.m3kg, R.ft3lb
    5e-05 0.0128295584431
    """
    __title__ = QApplication.translate("pychemqt", "Specific Volume")
    rates = {"m3kg": 1.,
             "lg": 1.,
             "lkg": k.liter**2,
             "ccg": k.liter**2,
             "mlg": k.liter**2,
             "m3g": k.gram**2,
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
    __text__ = [u'm³/kg', u'cm³/g', u'ml/g', u'm³/g', u'cm³/kg', u'ft³/lb',
                u'in³/lb', 'gallon UK/lb', 'gallon US/lb', 'barril/lb',
                u'ft³/ton UK', u'ft³/ton US', u'ft³/slug', u'ft³/onza',
                u'in³/onza', 'gallon UK/onza', 'gallon US/onza']
    __units__ = ['m3kg', 'lkg', 'ccg', 'mlg', 'm3g', 'cckg', 'ft3lb',
                 'inch3lb', 'galUKlb', 'galUSlb', 'bbllb', 'ft3tonUK',
                 'ft3tonUS', 'ft3slug',  'ft3oz', 'in3oz', 'galUKoz', 'galUSoz']
    __units_set__ = {"altsi": "m3kg", "si": "m3kg", "metric": "m3kg",
                     "cgs": "ccg", "english": "ft3lb"}

    def __init__(self, data, unit="m3kg", magnitud=""):
        super(SpecificVolume_square, self).__init__(data, unit, magnitud)


class MolarVolume(unidad):
    """Class that models a specific molar volume measure
    Supported units:

    * cubic meter per kilomol (m3kmol) default
    * liter per kilomol (lkmol) (same as cc/mol, and  ml/mol)
    * cubic meter per mol (m3mol)
    * cubic centimeter per kilomol (cckmol)
    * cubic foot per poundmol (ft3lbmol)
    * cubic inch per poundmol (in3lbmol)

    >>> R=MolarVolume(50, "lkmol")
    >>> print  R.m3kmol, R.ft3lbmol
    0.05 0.800923168698
    """
    __title__ = QApplication.translate("pychemqt", "Molar Volume")
    rates = {"m3kmol": 1.,
             "lmol": 1.,
             "lkmol": k.liter,
             "ccmol": k.liter,
             "mlmol": k.liter,
             "m3mol": k.gram,
             "cckmol": k.micro,
             "ft3lbmol": k.foot**3/k.pound,
             "inch3lbmol": k.inch**3/k.pound}
    __text__ = [u'm³/kmol', u'l/mol', u'l/kmol', u'cm³/mol', u'ml/mol',
                u'm³/mol', u'cm³/kmol', u'ft³/lbmol', u'in³/lbmol']
    __units__ = ['m3kmol', 'lmol', 'lkmol', 'ccmol', 'mlmol', 'm3mol',
                 'cckmol', 'ft3lbmol', 'inch3lbmol']
    __units_set__ = {"altsi": "m3kmol", "si": "m3kmol", "metric": "m3kmol",
                     "cgs": "ccmol", "english": "ft3lbmol"}

    def __init__(self, data, unit="m3kmol", magnitud=""):
        super(MolarVolume, self).__init__(data, unit, magnitud)


class Density(unidad):
    """Class that models a density measure
    Supported units:

    * kilogram per cubic foot (kgm3) default (same as gr/l)
    * kilogram per liter (kgl) (same as g/cc, g/ml)
    * gram per cubic meter (gm3)
    * kilogram per cubic centimeter (kgcc)
    * pound per cubic foot (lbft3)
    * pound per cubic inch (lbin3)
    * pound per british gallon (lbgalUK)
    * pound per gallon (lbgalUS)
    * pound per barril (lbbbl)
    * short tonne per cubic foot (tonUKft3)
    * long tonne per cubic foot (tonUSft3)
    * slug per cubic foot (slugft3)
    * ounce per cubic foot (ozft3)
    * ounce per cubic inch (ozin3)
    * ounce per british gallon (ozgalUK)
    * ounce per gallon (ozgalUS)

    >>> R=Density(1, "kgl")
    >>> print R.kgm3, R.lbft3
    1000.0 62.4279605761
    """
    __title__ = QApplication.translate("pychemqt", "Density")
    rates = {"kgm3": 1.,
             "gl": 1.,
             "kgl": 1./k.liter,
             "gcc": 1./k.liter,
             "gml": 1./k.liter,
             "gm3": 1./k.gram,
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
    __text__ = [u'kg/m³', u'g/cm³', u'g/m³', u'kg/cm³', u'lb/ft³', u'lb/inch³',
                'lb/galon UK', 'lb/galon US', 'lb/barril', u'ton UK/ft³',
                u'ton US/ft³', u'slug/ft³', u'onza/ft³', u'onza/inch³',
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

    def __init__(self, data, unit="kgm3", magnitud=""):
        super(Density, self).__init__(data, unit, magnitud)


class MolarDensity(unidad):
    """Class that models a molar density measure
    Supported units:

    * kilomol per cubic foot (kmolm3) default (same as mol/l)
    * kilomol per liter (kmoll) (same as mol/cc, mol/ml)
    * mol per cubic meter (molm3)
    * kilomol per cubic centimeter (kmolcc)
    * lbmol per cubic foot (lbft3)
    * lbmol per cubic inch (lbin3)

    >>> R=MolarDensity(1, "kmolm3")
    >>> print R.molcc, R.lbmolft3
    1.0 0.0624279605761
    """
    __title__ = QApplication.translate("pychemqt", "Molar Density")
    rates = {"kmolm3": 1.,
             "molcc": 1.,
             "kmoll": 1./k.liter,
             "molm3": 1./k.gram,
             "kmolcc": 1./k.micro,
             "lbmolft3": k.pound/k.foot**3,
             "lbmolin3": k.pound/k.inch**3}
    __text__ = [u'kmol/m³', u'mol/cm³', u'mol/m³', u'kmol/cm³', u'lbmol/ft³',
                u'lbmol/inch³']
    __units__ = ['kmolm3', 'molcc', 'molm3', 'kmolcc', 'lbmolft3', 'lbmolin3']
    __units_set__ = {"altsi": "kmolm3", "si": "kmolm3", "metric": "kmolm3",
                     "cgs": "molcc", "english": "lbmolft3"}

    def __init__(self, data, unit="kmolm3", magnitud=""):
        super(MolarDensity, self).__init__(data, unit, magnitud)


class Force(unidad):
    """Class that models a force measure
    Supported units:

    * Newton (N) default
    * Kilonewton (kN)
    * dyn (dyn)
    * kilogram force (kgf)
    * gram force (gf)
    * pound force (lbf)
    * ounce force (ozf)
    * poundal (pdl)
    * Long tonne force (TonfUS)
    * Short tonne force (TonfUK)

    >>> F=Force(1, "pdl")
    >>> print F.N, F.kgf, F.dyn
    0.138254954376 0.0140980818502 13825.4954376
    """
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
    __units_set__ = {"altsi": "kN", "si": "N", "metric": "N", "cgs": "dyn",
                     "english": "lbf"}

    def __init__(self, data, unit="N", magnitud=""):
        super(Force, self).__init__(data, unit, magnitud)


class Pressure(unidad):
    """Class that models a pressure measure
    Supported units:

    * Pascal (Pa) default
    * Megapascal (MPa)
    * Hectopascal (hPa)
    * Kilopascal (kPa)
    * Bar (bar)
    * Bar gauge (barg)
    * Milibar (mbar)
    * Pound per square inch (psi)
    * Pound per square inch gauge (psig)
    * Atmosphere (atm)
    * Atmosphere technical, kg/cm2 (kgcm2)
    * Atmosphere technical gauge (kgcm2g)
    * Milimeter of water column (mmH2O)
    * Meter of water column (mH2O)
    * Centimeter of water column (cmH2O)
    * Inch of water column (inH2O)
    * Foot of water column (ftH2O)
    * Milimeter of mercury column (mmHg)
    * Torricelli (torr)
    * Centimeter of mercury column (cmHg)
    * Inch of mercury column (inHg)
    * Foot of mercury column (ftHg)
    * Pound per cubic curadrado (lbcm2)
    * Pound per cubic foot (lbft2)
    * Dyn per cubic centimeter (dyncm2)

    >>> P=Pressure(760, "mmHg")
    >>> print P.bar, P.atm, P.psi, P.kgcm2g
    1.01325 1.0 14.6959487755 0.0
    """
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
                'psi g', 'atm', u'kg/cm²', u'kg/cm² g', 'mmH2O', 'cmH2O',
                'mH2O', 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg',
                u'lb/cm²', u'lb/ft²', u'dyn/cm²']
    __units__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'barg', 'mbar', 'psi',
                 'psig', 'atm', 'kgcm2', 'kgcm2g', 'mmH2O', 'cmH2O', 'mH2O',
                 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg', 'lbcm2',
                 'lbft2', 'dyncm2']
    __units_set__ = {"altsi": "bar", "si": "Pa", "metric": "Pa",
                     "cgs": "dyncm2", "english": "psi"}

    def __init__(self, data, unit="Pa", magnitud=""):
        if not magnitud:
            magnitud = self.__class__.__name__
        self.magnitud = magnitud

        if unit == "conf":
            unit = self.__units__[self.Config.getint('Units', magnitud)]

        for key in self.__class__.rates:
            self.__setattr__(key, self._data/self.__class__.rates[key])

        self.barg = (self.Pa-k.atm)/k.bar
        self.psig = (self.Pa-k.atm)/k.psi
        self.kgcm2g = (self.Pa-k.atm)*k.centi**2/k.g

    def __new__(cls, data, unit="Pa", magnitud=""):
        cls.Config = getMainWindowConfig()

        if data is None:
            data = 0
            cls.code = "n/a"
        else:
            cls.code = ""
            data = float(data)

        if unit == "conf":
            unit = cls.__units__[cls.Config.getint('Units', 'Pressure')]

        if unit == "barg":
            cls._data = data*k.bar+k.atm
        elif unit == "psig":
            cls._data = data*k.psi+k.atm
        elif unit == "kgcm2g":
            cls._data = data*k.g/k.centi**2+k.atm
        else:
            cls._data = data * cls.rates[unit]

        return float.__new__(cls, cls._data)


class DeltaP(unidad):
    """Class that models a delta pressure measure
    Supported units:

    * Pascal (Pa) default
    * Megapascal (MPa)
    * Hectopascal (hPa)
    * Kilopascal (kPa)
    * Bar (bar)
    * Bar gauge (barg)
    * Milibar (mbar)
    * Pound per square inch (psi)
    * Pound per square inch gauge (psig)
    * Atmosphere (atm)
    * Atmosphere technical, kg/cm2 (kgcm2)
    * Atmosphere technical gauge (kgcm2g)
    * Milimeter of water column (mmH2O)
    * Meter of water column (mH2O)
    * Centimeter of water column (cmH2O)
    * Inch of water column (inH2O)
    * Foot of water column (ftH2O)
    * Milimeter of mercury column (mmHg)
    * Torricelli (torr)
    * Centimeter of mercury column (cmHg)
    * Inch of mercury column (inHg)
    * Foot of mercury column (ftHg)
    * Pound per cubic curadrado (lbcm2)
    * Pound per cubic foot (lbft2)
    * Dyn per cubic centimeter (dyncm2)

    >>> P=Pressure(760, "mmHg")
    >>> print P.bar, P.atm, P.psi, P.kgcm2g
    1.01325 1.0 14.6959487755 0.0
    """
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
                'psi g', 'atm', u'kg/cm²', u'kg/cm² g', 'mmH2O', 'cmH2O',
                'mH2O', 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg',
                u'lb/cm²', u'lb/ft²', u'dyn/cm²']
    __units__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'barg', 'mbar', 'psi',
                 'psig', 'atm', 'kgcm2', 'kgcm2g', 'mmH2O', 'cmH2O', 'mH2O',
                 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg', 'lbcm2',
                 'lbft2', 'dyncm2']
    __units_set__ = {"altsi": "bar", "si": "Pa", "metric": "Pa",
                     "cgs": "dyncm2", "english": "psi"}

    def __init__(self, data, unit="Pa", magnitud=""):
        super(DeltaP, self).__init__(data, unit, magnitud)


class Energy(unidad):
    """Class that models a energy measure
    Supported units:

    * Joule (J) default
    * Kilojoule (kJ)
    * Megajoule (MJ)
    * Calorie (cal)
    * Kilocalorie (kcal)
    * Calorie international (cal_i)
    * Erg (erg)
    * Btu (Btu)
    * kiloBtu (kBtu)
    * MegaBtu (MBtu)
    * Watt-hour (Wh)
    * Kilowatt-hour (KWh)
    * Megawatt-hour (MWh)
    * TonTNT (TNT)
    * Horsepower-hour (HPh)
    * Metric horsepower·hour (CVh)
    * Kilogram force per meter (kgfm)
    * Pound force per foot (lbfft)
    * Gigaelectronvolt (GeV)
    * Barrel petrol (oil)
    * Tonne of oil equivalent (toe)
    * Tonne of coal equivalent (tce)

    >>> E=Energy(1, "kcal")
    >>> print E.J, E.Btu, E.Wh
    4184.0 3.96566683139 1.16222222222
    """
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
    __units_set__ = {
        "Energy": {"altsi": "MJ", "si": "MJ", "metric": "J", "cgs": "erg",
                   "english": "MBtu"},
        "Work": {"altsi": "MJ", "si": "kWh", "metric": "J", "cgs": "erg",
                 "english": "HPh"}}

    def __init__(self, data, unit="J", magnitud=""):
        super(Energy, self).__init__(data, unit, magnitud)


class Enthalpy(unidad):
    """Class that models a enthalpy measure
    Supported units:

    * Joule per kilogram (Jkg) default
    * Kilojoule per kilogram (kJkg)
    * Megajoule per kilogram (MJkg)
    * Kilowatt hour per kilogram (kWhkg)
    * Calorie per kilogram (calkg)
    * Kilocalorie per kilogram (kcalkg)
    * Calorie per gram (calg)
    * Calorie per pound (callb)
    * Kilocalorie per gram (kcalg)
    * Btu per pound (Btulb)

    >>> H=Enthalpy(-5, "Btulb")
    >>> print H.kJkg, H.kcalkg
    -11.63 -2.77963671128
    """
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

    def __init__(self, data, unit="Jkg", magnitud=""):
        super(Enthalpy, self).__init__(data, unit, magnitud)


class MolarEnthalpy(unidad):
    """Class that models a enthalpy measure in molar base
    Supported units:

    * Joule per kilogram (Jkmol) default
    * Kilojoule per kilogram (kJkmol)
    * Megajoule per kilogram (MJkmol)
    * Kilowatt hour per kilogram (kWhkmol)
    * Calorie per kilogram (calkmol)
    * Kilocalorie per kilogram (kcalkmol)
    * Calorie per gram (calmol)
    * Calorie per pound (callbmol)
    * Kilocalorie per gram (kcalmol)
    * Btu per pound (Btulbmol)

    >>> H=MolarEnthalpy(-5, "Btulbmol")
    >>> print H.kJkmol, H.kcalkmol
    -11.63 -2.77963671128
    """
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

    def __init__(self, data, unit="Jkmol", magnitud=""):
        super(MolarEnthalpy, self).__init__(data, unit, magnitud)


class Entropy(unidad):
    """Class that models a entropy measure
    Supported units:

    * Joule per kelvin (JK) default
    * Kilojoule per kelvin(kJK)
    * Megajoule per kelvin (MJK)
    * Calorie per kelvin (calK)
    * Kilocalorie per kelvin (kcalK)
    * Megacalorie per kelvin (McalK)
    * Watt hour per kelvin (WhK)
    * Kilowatt hour per kelvin (kWhK)
    * Megawatt hour per kelvin (MWhK)
    * Horsepower hour per fahrenheit (HPhF)
    * Btu per fahrenheit (BtuF)
    * KiloBtu per fahrenheit (kBtuF)
    * MegaBtu per fahrenheit (MBtuF)

    >>> S=Entropy(30, "BtuF")
    >>> print S.kJK, S.kcalK
    56.9730160415 13.616877639
    """
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

    def __init__(self, data, unit="JK", magnitud=""):
        super(Entropy, self).__init__(data, unit, magnitud)


class SpecificHeat(unidad):
    """Class that models a specific heat measure
    Supported units:

    * Joule per kilogram kelvin (JkgK) default
    * Kilojoule per kilogram kelvin (kJkgK)
    * Joule per gram kelvin (JgK)
    * Calorie per kilogram kelvin (kcalkgK)
    * Calorie per gram kelvin (calgK)
    * Kilocalorie per gram kelvin (kcalgK)
    * Kilowatt hour per kilogram kelvin (kWhkgK)
    * Btu per pound fahrenheit (BtulbF)

    >>> C=SpecificHeat(1, "BtulbF")
    >>> print C.kJkgK, C.kcalkgK
    4.1868 1.00066921606
    """
    __title__ = QApplication.translate("pychemqt", "Specific Heat")
    rates = {"JkgK": 1.,
             "kJkgK": k.kilo,
             "JgK": k.kilo,
             "kcalkgK": k.calorie*k.kilo,
             "calgK": k.calorie*k.kilo,
             "kcalgK": k.calorie*k.kilo**2,
             "kWhkgK": k.kilo*k.hour,
             "BtulbF": k.Btu/k.lb/k.Rankine}
    __text__ = [u'J/kg·K', u'kJ/kg·K', u'kcal/kg·K', u'cal/g·K', u'kcal/g·K',
                u'kWh/kg·K', u'Btu/lb·F']
    __units__ = ['JkgK', 'kJkgK', 'kcalkgK', 'calgK', 'kcalgK', 'kWhkgK',
                 'BtulbF']
    __units_set__ = {"altsi": "kJkgK", "si": "JkgK", "metric": "JkgK",
                     "cgs": "calgK", "english": "BtulbF"}

    def __init__(self, data, unit="JkgK", magnitud=""):
        super(SpecificHeat, self).__init__(data, unit, magnitud)


class Power(unidad):
    """Class that models a power measure
    Supported units:

    * Watt (W) default
    * Kilowattt (kW)
    * Megawatt (MW)
    * Horsepower (hp)
    * Metric horsepower (CV)
    * Calorie per second (cals)
    * Kilocalorie per hour (kcalh)
    * Joule per hour (Jh)
    * Erg per second (ergs)
    * Btu per hour (Btuh)
    * Megabtu per hour (MBtuh)
    * Btu per minute (Btumin)
    * Btu per second (Btus)
    * Foot pound force per hour (ftlbfh)
    * Foot pound force per minute (ftlbfmin)
    * Foot pound force per second (ftlbfs)

    >>> W=Power(5, "Btuh")
    >>> print W.kW, W.hp, W.kcalh
    0.00146535535086 0.00196507389461 1.26082200361
    """
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
                u'ft/lbf·s', u'ft/lbf·min', u'ft/lbf·h']
    __units__ = ['W', 'kW', 'MW', 'hp', 'CV', 'cals', 'kcalh', 'Jh', 'kJh',
                 'MJh', 'ergs', 'Btus', 'Btumin', 'Btuh', 'MBtuh', 'ftlbfs',
                 'ftlbfmin', 'ftlbfh']
    _magnitudes = [
        ("EnergyFlow", QApplication.translate("pychemqt", "Energy Flow")),
        ("Power", QApplication.translate("pychemqt", "Power"))]
    __units_set__ = {
        "EnergyFlow": {"altsi": "MJh", "si": "kJh", "metric": "Jh",
                       "cgs": "ergs", "english": "MBtuh"},
        "Power": {"altsi": "hp", "si": "kW", "metric": "Jh", "cgs": "ergs",
                  "english": "hp"}}

    def __init__(self, data, unit="W", magnitud=""):
        super(Power, self).__init__(data, unit, magnitud)


class MassFlow(unidad):
    """Class that models a mass flow measure
    Supported units:

    * kg per second (kgs) default
    * kg per minute (kgmin)
    * kg per hour (kgh)
    * g per second (gs)
    * g per minute (gmin)
    * g per hour (gh)
    * Tonne per second (Tons)
    * Tonne per minute (Tonmin)
    * Tonne per hour (Tonh)
    * pound per second (lbs)
    * pound per minute (lbmin)
    * pound per hour (lbh)
    * Short Tonne per second (TonUKs)
    * Long Tonne per second (TonUSs)
    * Short Tonne per minute (TonUKmin)
    * Long Tonne per minute (TonUSmin)
    * Short Tonne per hour (TonUKh)
    * Long Tonne per hour (TonUSh)
    * Short Tonne per day (TonUKday)
    * Long Tonne per day (TonUSday)
    * Short Tonne per year (TonUKyear)
    * Long Tonne per year (TonUSyear)

    >>> G=MassFlow(1, "gs")
    >>> print G.kgh, G.lbh, G.gmin
    3.6 7.93664143866 60.0
    """
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

    def __init__(self, data, unit="kgs", magnitud=""):
        super(MassFlow, self).__init__(data, unit, magnitud)


class MolarFlow(unidad):
    """Class that models a molar flow measure
    Supported units:

    * kmol per second (kmols) default
    * kmol per minute (kmolmin)
    * kmol per hour (kmolh)
    * gmol per second (mols)
    * gmol per minute (molmin)
    * gmol per hour (gmolh)
    * lbmol per second (lbmols)
    * lbmol per minute (lbmolmin)
    * lbmol per hour (lbmolh)

    >>> G=MolarFlow(1, "mols")
    >>> print G.kmolh, G.lbmolh, G.molmin
    3.6 7.93664143866 60.0
    """
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

    def __init__(self, data, unit="kmols", magnitud=""):
        super(MolarFlow, self).__init__(data, unit, magnitud)


class VolFlow(unidad):
    """Class that models a volumetric flow measure
    Supported units:

    * cubic meters per second (m3s) default
    * cubic meters per minute (m3min)
    * cubic meters per hour (m3h)
    * liters per second (ls)
    * liters per minute (lmin)
    * liters per hour (lh)
    * cubic centimeters per second (ccs, cm3s)
    * cubic centimeters per minute (ccmin)
    * cubic centimeters per hour (cch)
    * cubic foot per second (ft3s)
    * cubic foot per minute (ft3min)
    * kilo cubic foot per minute (kft3min)
    * cubic foot per hour (ft3h)
    * mega cubic foot per day (mft3day)
    * British gallon per hour (galUKh)
    * Gallon per hour (galUSh)
    * British gallon per minute (galUKmin)
    * Gallon per minute (galUSmin)
    * British gallon per second (galUKs)
    * Gallon per second (galUSs)
    * Barrel per day (bblday)
    * Barrel per hour (bblh)
    * Barrel per minute (bblmin)
    * Barrel per second (bbls)

    >>> V=VolFlow(1, "lmin")
    >>> print V.m3h, V.ft3min, V.ccs
    0.06 0.0353146667215 16.6666666667
    """
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
    __text__ = [u'm³/s', u'm³/min', u'm³/h', 'l/s', 'l/min', 'l/h', u'cm³/s',
                u'cm³/min', u'cm³/h', u'ft³/s', u'ft³/min', u'kft³/min',
                u'ft³/h', u'Mft³/day', 'galon UK/h', 'galon US/h',
                'galon UK/min', 'galon US/min', 'galon UK/s', 'galon US/s',
                'barril/s', 'barril/min', 'barril/h', 'barril/day']
    __units__ = ['m3s', 'm3min', 'm3h', 'ls', 'lmin', 'lh', 'ccs', 'ccmin',
                 'cch', 'ft3s', 'ft3min', 'kft3min',  'ft3h', 'mft3day',
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

    def __init__(self, data, unit="m3s", magnitud=""):
        super(VolFlow, self).__init__(data, unit, magnitud)


class Diffusivity(unidad):
    """Class that models a diffusivity measure (useful for kinematic viscosity)
    Supported units:

    * Square meter per second (m2s) default
    * Square centimeter per second (cm2s)
    * Square milimeter per second (mm2s)
    * Square foot per second (ft2s)
    * Square inch per second (inch2s)
    * Square meter per hour (m2h)
    * Square foot per hour (ft2h)
    * Square inch per hour (inch2h)
    * Stokes (St)
    * Centistokes (cSt)

    >>> k=Diffusivity(5, "St")
    >>> print k.m2s, k.ft2s
    0.0005 0.00538195520835
    """
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
    __text__ = [u"m²/s", u"cm²/s", u"mm²/s", u"ft²/s", u"inch²/s", u"m²/h",
                u"ft²/h", u"inch²/h", "St", "cSt"]
    __units__ = ["m2s", "cm2s", "mm2s", "ft2s", "inch2s", "m2h", "ft2h",
                 "inch2h", "St", "cSt"]
    _magnitudes = [
        ("Diffusivity", QApplication.translate("pychemqt", "Diffusivity")),
        ("KViscosity", QApplication.translate("pychemqt", "Kinematic viscosity"))]
    __units_set__ = {
        "Diffusivity": {"altsi": "m2s", "si": "m2s", "metric": "m2s",
                        "cgs": "cm2s", "english": "ft2s"},
        "KViscosity": {"altsi": "m2s", "si": "m2s", "metric": "m2s",
                       "cgs": "cm2s", "english": "ft2s"}}

    def __init__(self, data, unit="m2s", magnitud=""):
        super(Diffusivity, self).__init__(data, unit, magnitud)


class HeatFlux(unidad):
    """Class that models a heat flux measure
    Supported units:

    * W/m2 (Wm2) default
    * kW/m2 (kWm2)
    * cal/hm2 (calhm2)
    * cal/sm2 (calsm2)
    * cal/scm2 (calscm2)
    * kcal/m2h (kcalhm2)
    * Btu/ft2h (Btuhft2)
    * Btu/ft2s (Btusft2)

    >>> H=HeatFlux(1, "Btuhft2")
    >>> print H.Wm2, H.kcalhm2
    3.15459074506 2.71427501965
    """
    __title__ = QApplication.translate("pychemqt", "Heat Flux")
    rates = {"Wm2": 1.,
             "kWm2": k.kilo,
             "calhm2": k.calorie/k.hour,
             "calsm2": k.calorie,
             "calscm2": k.calorie/k.centi**2,
             "kcalhm2": k.kilo*k.calorie/k.hour,
             "Btuhft2": k.Btu/k.hour/k.foot**2, "Btusft2": k.Btu/k.foot**2}
    __text__ = [u"W/m²", u"kW/m²", u"cal/hm²", u"cal/sm²", u"cal/scm²",
                u"kcal/hm²", u"Btu/hft²", u"Btu/sft²"]
    __units__ = ["Wm2", "kWm2", "calhm2", "calsm2", "calscm2", "kcalhm2",
                 "Btuhft2", "Btusft2"]
    __units_set__ = {"altsi": "Wm2", "si": "Wm2", "metric": "Wm2",
                     "cgs": "calscm2", "english": "Btuhft2"}

    def __init__(self, data, unit="Wm2", magnitud=""):
        super(HeatFlux, self).__init__(data, unit, magnitud)


class ThermalConductivity(unidad):
    """Class that models a thermal conductivity measure
    Supported units:

    * Watt per meter Kelvin (WmK) default
    * Joule per hour meter Kelvin (JhmK)
    * Kilojoule per hour meter Kelvin (kJhmK)
    * Calorie per second, centimeter kelvin (calscmK)
    * Calorie per hour, centimeter kelvin (calhcmK)
    * Calorie per hour, milimeter kelvin (calhmmK)
    * Kilocalorie per hour meter kelvin (kcalhmK)
    * Pound force per second fahrenheit (lbfsF)
    * Pound pie per square second fahrenheit (lbfts3F)
    * Btu per hour foot fahrenheit (BtuhftF)

    >>> k=ThermalConductivity(50)
    >>> print k.WmK, k.BtuhftF, k.kcalhmK
    50.0 28.8894658271 43.0210325048
    """
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
    __text__ = [u'W/m·K', u'mW/m·K', u"kW/m·K", u'J/h·m·K', u'cal/s·cm·K',
                u'cal/h·cm·K', u'kcal/h·m·K', u'lbf/s·F', u'lb/ft·s³·F',
                u'Btu/h·ft·F']
    __units__ = ['WmK', 'mWmK', "kWmK", 'JhmK', 'calscmK', 'calhcmK',
                 'kcalhmK', 'lbfsF', 'lbfts3F', 'BtuhftF']
    __units_set__ = {"altsi": "mWmK", "si": "WmK", "metric": "WmK",
                     "cgs": "calscmK", "english": "lbfts3F"}

    def __init__(self, data, unit="WmK", magnitud=""):
        super(ThermalConductivity, self).__init__(data, unit, magnitud)


class UA(unidad):
    """Class that models a UA measure
    Supported units:

    * Watt per Kelvin (WK) default
    * kilowatt per Kelvin (kWK)
    * miliwatt per kelvin (mWK)
    * Joule per hour Kelvin (JhK)
    * Kilojoule per hour Kelvin (kJhK)
    * Calorie per hour Kelvin (calhK)
    * Kilocalorie per hour Kelvin (kcalhK)
    * Calorie per second Kelvin (calsK)
    * Kilocalorie per second Kelvin (kcalsK)
    * Btu per hour fahrenheit (BtuhF)
    * Btu per second fahrenheit (BtusF)

    >>> h=UA(1, "BtuhF")
    >>> print h.WK, h.kcalhK
    5.67826334111 4.88569503537
    """
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
    __text__ = [u'W/K', u'kW/K', u'mW/K', u'J/h·K', u'kJ/h·K', u'cal/h·K',
                u'kcal/h·K', u'cal/s·K', u'kcal/s·K', u'Btu/h·F', u'Btu/s·F']
    __units__ = ['WK', 'kWK', 'mWK', 'JhK', 'kJhK', 'calhK', 'kcalhK',
                 'calsK', 'kcalsK', 'BtuhF', 'BtusF']
    __units_set__ = {"altsi": "mWK", "si": "WK", "metric": "WK",
                     "cgs": "calsK", "english": "BtuhF"}

    def __init__(self, data, unit="WK", magnitud=""):
        super(UA, self).__init__(data, unit, magnitud)


class HeatTransfCoef(unidad):
    """Class that models a heat transfer coefficient measure
    Supported units:

    * Watt per m2 Kelvin (Wm2K) default
    * kilowatt per m2 Kelvin (kWm2K)
    * Joule per hour m2 Kelvin (Jhm2K)
    * Kilojoule per hour m2 Kelvin (kJhm2K)
    * Calorie per hour m2 Kelvin (calhm2K)
    * Kilocalorie per hour m2 Kelvin (kcalhm2K)
    * Calorie per second m2 Kelvin (calsm2K)
    * Kilocalorie per second m2 Kelvin (kcalsm2K)
    * Calorie per second cm2 Kelvin (calscm2K)
    * Kilocalorie per second cm2 Kelvin (kcalscm2K)
    * Btu per hour foot2 fahrenheit (Btuhft2F)
    * Btu per second foot2 fahrenheit (Btusft2F)

    >>> h=HeatTransfCoef(1, "Btuhft2F")
    >>> print h.Wm2K, h.kcalhm2K
    5.67826334111 4.88569503537
    """
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
    __text__ = [u'W/m²·K', u'kW/m²·K', u'J/h·m²·K', u'kJ/h·m²·K', u'cal/h·m³·K',
                u'kcal/h·m²·K', u'cal/s·m²·K', u'kcal/s·m²·K', u'cal/s·cm²·K',
                u'kcal/s·cm²·K' , u'Btu/h·ft²·F',u'Btu/s·ft²·F']
    __units__ = ['Wm2K', 'kWm2K', 'Jhm2K', 'kJhm2K', 'calhm2K', 'kcalhm2K',
                 'calsm2K', 'kcalsm2K', 'calscm2K', 'kcalscm2K' , 'Btuhft2F',
                 'Btusft2F']
    __units_set__ = {"altsi": "Wm2K", "si": "Wm2K", "metric": "Wm2K",
                     "cgs": "calscm2K", "english": "Btuhft2F"}

    def __init__(self, data, unit="Wm2K", magnitud=""):
        super(HeatTransfCoef, self).__init__(data, unit, magnitud)


class Fouling(unidad):
    """Class that models a fouling factor resistence, inverse of heat
    transmision coefficient
    Supported units:

    * m2 Kelvin per Watt (m2KW) default
    * m2 Kelvin per kilowatt (m2KkW)
    * hour m2 Kelvin per Joule (hm2KJ)
    * hour m2 Kelvin per Kilojoule (hm2KkJ)
    * hour m2 Kelvin per Calorie (hm2Kcal)
    * hour m2 Kelvin per Kilocalorie (hm2Kkcal)
    * second m2 Kelvin per Calorie (sm2Kcal)
    * second m2 Kelvin per Kilocalorie (sm2Kkcal)
    * second cm2 Kelvin per Calorie (scm2Kcal)
    * second cm2 Kelvin per Kilocalorie (scm2Kkcal)
    * hour foot2 fahrenheit per Btu (hft2FBtu)
    * second foot2 fahrenheit per Btu (sft2FBtu)

    >>> h=Fouling(1, "hft2FBtu")
    >>> print h.m2KW, h.hm2Kkcal
    0.176110183682 0.204679169035
    """
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
    __text__ = [u'm²·K/W', u'm²·K/kW', u'h·m²·K/J', u'h·m²·K/kJ',
                u'h·m³·K/cal', u'h·m²·K/kcal', u's·m²·K/cal', u's·m²·K/kcal',
                u's·cm²·K/cal', u's·cm²·K/kcal', u'h·ft²·F/Btu', u's·ft²·F/Btu']
    __units__ = ['m2KW', 'm2KkW', 'hm2KJ', 'hm2KkJ', 'hm2Kcal', 'hm2Kkcal',
                 'sm2Kcal', 'sm2Kkcal', 'scm2Kcal', 'scm2Kkcal', 'hft2FBtu',
                 'sft2FBtu']
    __units_set__ = {"altsi": "m2KW", "si": "m2KW", "metric": "m2KW",
                     "cgs": "scm2Kkcal", "english": "hft2FBtu"}

    def __init__(self, data, unit="m2KW", magnitud=""):
        super(Fouling, self).__init__(data, unit, magnitud)


class Tension(unidad):
    """Class that models a surface tension measure
    Supported units:

    * Newton per meter (Nm) default
    * Dyn per centimeter (dyncm)
    * Pound force per foot (lbfft)

    >>> s=Tension(1, "lbfft")
    >>> print s.Nm, s.dyncm
    14.5939029372 14593.9029372
    """
    __title__ = QApplication.translate("pychemqt", "Surface Tension")
    rates = {"Nm": 1.,
             "mNm": k.milli,
             "dyncm": k.dyn/k.centi,
             "lbfft": k.lbf/k.foot}
    __text__ = ['N/m', 'mN/m', 'dyn/cm', 'lbf/ft']
    __units__ = ['Nm', 'mNm', 'dyncm', 'lbfft']
    __units_set__ = {"altsi": "mNm", "si": "Nm", "metric": "Nm",
                     "cgs": "dyncm", "english": "lbfft"}

    def __init__(self, data, unit="Nm", magnitud=""):
        super(Tension, self).__init__(data, unit, magnitud)


class Viscosity(unidad):
    """Class that models a viscosity measure
    Supported units:

    * Pascal per second (Pas) default
    * Milipascal per second (mPas)
    * Micropascal per second (muPas)
    * Dinas second per square centimeter (dynscm2)
    * Poise (P)
    * Centipoise (cP)
    * Reyn (reyn)
    * Pound per foot second (lbfts)
    * Pound force second per square foot (lbfft2)
    * Pound force second per square inch (lbfinch2)
    * Pound per foot hour (lbfth)

    >>> m=Viscosity(1, "P")
    >>> print m.cP, m.Pas, m.lbfth
    100.0 0.1 241.90883105
    """
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
    __text__ = [u'Pa·s', u'mPa·s', u'µPa·s', 'P', 'cP', u'dyn/s·cm²', u'µP',
                'reyn', u'lb/ft·s', u'lbf/ft²', u'lbf/in²', u'lb/ft·h']
    __units__ = ['Pas', 'mPas', 'muPas', 'P', 'cP', 'dynscm2', 'microP',
                 'reyn', 'lbfts', 'lbfft2', 'lbfinch2', 'lbfth']
    __units_set__ = {"altsi": "muPas", "si": "Pas", "metric": "Pas",
                     "cgs": "dynscm2", "english": "cP"}

    def __init__(self, data, unit="Pas", magnitud=""):
        super(Viscosity, self).__init__(data, unit, magnitud)


class SolubilityParameter(unidad):
    """Class that models a solubility parameter measure
    Supported units:

    * raiz(J/m3) (Jm3) default
    * raiz(cal/cc) (calcc)
    * raiz(Btuft3) (Btuft3)

    >>> S=SolubilityParameter(1, "Btuft3")
    >>> print S.Jm3, S.calcc
    193.025764622 0.0943668467764
    """
    __title__ = QApplication.translate("pychemqt", "Solubility Parameter")
    rates = {"Jm3": 1.,
             "calcc": (k.calorie*k.mega)**0.5,
             "Btuft3": (k.Btu*k.foot**-3)**0.5}
    __text__ = [u"(J/m³)^0.5", u"(cal/cm³)^0.5", u"(Btu/ft³)^0.5"]
    __units__ = ["Jm3", "calcc", "Btuft3"]
    __units_set__ = {"altsi": "Jm3", "si": "Jm3", "metric": "Jm3",
                     "cgs": "calcc", "english": "Btuft3"}

    def __init__(self, data, unit="Jm3", magnitud=""):
        super(SolubilityParameter, self).__init__(data, unit, magnitud)


class PotencialElectric(unidad):
    """Class that models a potencia electric measure
    Supported units:
    * Volt per meter (Vm) default
    * Kilovolt per meter (kVm)
    * Megavolt per meter (MVm)
    * Volt per centimeter (Vcm)
    * Kilovolt per centimeter (kVcm)
    * Megavolt per centimeter (MVcm)
    * Statvolt per meter (statVm)
    * statvolt per centimeter (statVcm)

    >>> e=PotencialElectric(3, "statVcm")
    >>> print e.Vm
    90000.0
    """
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
    __units__ = ["Vm", "kVm", "MVm", "Vcm", "kVcm", "MVcm", "statVm", "statVcm"]
    __units_set__ = {"altsi": "Vm", "si": "Vm", "metric": "Vm", "cgs": "Vcm",
                     "english": "statVm"}

    def __init__(self, data, unit="Vm", magnitud=""):
        super(PotencialElectric, self).__init__(data, unit, magnitud)


class DipoleMoment(unidad):
    """Class that models a solubility parameter measure
    Supported units:
    * Coulomb per meter (Cm) default
    * debyes (Debye)

    >>> dp=DipoleMoment(1, "Debye")
    >>> print dp.Cm, dp.Debye
    3.33564095198e-30 1.0
    """
    __title__ = QApplication.translate("pychemqt", "Dipole Moment")
    rates = {"Cm": 1.,
             "Debye": k.debye}
    __text__ = [u'C·m', 'Debye']
    __units__ = ['Cm', 'Debye']
    __units_set__ = {"altsi": "Cm", "si": "Cm", "metric": "Cm", "cgs": "Cm",
                     "english": "Debye"}

    def __init__(self, data, unit="Cm", magnitud=""):
        super(DipoleMoment, self).__init__(data, unit, magnitud)


class CakeResistance(unidad):
    """Class that models a a cake resistence measure
    Supported units:
    * meter per kilogram (mkg) default
    * cm per gram (cmg)
    * foot per pound (ftlb)

    >>> dp=CakeResistance(1, "ftlb")
    >>> print dp.mkg
    0.67196897514
    """
    __title__ = QApplication.translate("pychemqt", "Cake Resistance")
    rates = {"mkg": 1.,
             "cmg": k.centi/k.kilo,
             "ftlb": k.foot/k.pound}
    __text__ = ['m/kg', 'c/gr', "ft/lb"]
    __units__ = ['mkg', 'cmg', "ftlb"]
    __units_set__ = {"altsi": "mkg", "si": "mkg", "metric": "mkg",
                     "cgs": "cmg", "english": "ftlb"}

    def __init__(self, data, unit="mkg", magnitud=""):
        super(CakeResistance, self).__init__(data, unit, magnitud)


class PackingDP(unidad):
    """Class that models a packing drop pressure measure
    Supported units:
    * mmH2O per meter (mmH2Om) default
    * inH2O per foot (inH2Oft)

    >>> dp=PackingDP(1, "inH2Oft")
    >>> print dp.mmH2Om
    83.3333333333
    """
    __title__ = QApplication.translate("pychemqt", "Packing Pressure drop")
    rates = {"mmH2Om": 1.,
             "inH2Oft": k.inch/k.milli/k.foot}
    __text__ = ['mmH2O/m', 'inH2O/ft']
    __units__ = ['mmH2Om', 'inH2Oft']
    __units_set__ = {"altsi": "mmH2Om", "si": "mmH2Om", "metric": "mmH2Om",
                     "cgs": "mmH2Om", "english": "inH2Oft"}

    def __init__(self, data, unit="mmH2Om", magnitud=""):
        super(PackingDP, self).__init__(data, unit, magnitud)


class V2V(unidad):
    """Class that models a volume ratio (Ratio gas-oil)
    Supported units:

    * cubic meter/cubic meter (m3m3) default
    * cubic foot/cubic foot (ft3ft3)
    * liter/liter (ll)
    * cubic foot/oil (ft3bbl)

    >>> V=V2V(1, "ft3bbl")
    >>> print V.m3m3
    0.178107606679
    """
    __title__ = QApplication.translate("pychemqt", "Gas-Oil ratio")
    rates = {"m3m3": 1.,
             "ft3ft3": 1.,
             "ll": 1.,
             "ft3bbl": k.foot**3/k.bbl}
    __text__ = [u"m³m³", u"ft/bbl"]
    __units__ = ["m3m3", "ft3bbl"]
    __units_set__ = {"altsi": "m3m3", "si": "m3m3", "metric": "m3m3",
                     "cgs": "m3m3", "english": "ft3bbl"}

    def __init__(self, data, unit="m3m3", magnitud=""):
        super(V2V, self).__init__(data, unit, magnitud)


class InvTemperature(unidad):
    """Class that models a  inverse temperature measure
    Supported units:

    * Kelvin (default)
    * Celsius
    * Fahrenheit
    * Rankine
    * Reaumur

    >>> T=InvTemperature(25, "C")
    >>> print T.K, T.F
    25.0 13.8888888889
    """
    __title__ = QApplication.translate("pychemqt", "Temperature inverse")
    rates = {"K": 1.,
             "C": 1.,
             "F": 1./k.Rankine,
             "R": 1./k.Rankine,
             "Re": 1./k.Reaumur}
    __text__ = ['1/K', u'1/ºC', u'1/ºR', u'1/ºF', u'1/ºRe']
    __units__ = ['K', 'C', 'R', 'F', 'Re']
    __units_set__ = {"altsi": "C", "si": "K", "metric": "C", "cgs": "C",
                     "english": "F"}

    def __init__(self, data, unit="K", magnitud=""):
        super(InvTemperature, self).__init__(data, unit, magnitud)


class InvPressure(unidad):
    """Class that models a inverse pressure measure
    Supported units:

    * Pascal (Pa) default
    * Megapascal (MPa)
    * Hectopascal (hPa)
    * Kilopascal (kPa)
    * Bar (bar)
    * Bar gauge (barg)
    * Milibar (mbar)
    * Pound per square inch (psi)
    * Pound per square inch gauge (psig)
    * Atmosphere (atm)
    * Atmosphere technical, kg/cm2 (kgcm2)
    * Atmosphere technical gauge (kgcm2g)
    * Milimeter of water column (mmH2O)
    * Meter of water column (mH2O)
    * Centimeter of water column (cmH2O)
    * Inch of water column (inH2O)
    * Foot of water column (ftH2O)
    * Milimeter of mercury column (mmHg)
    * Torricelli (torr)
    * Centimeter of mercury column (cmHg)
    * Inch of mercury column (inHg)
    * Foot of mercury column (ftHg)
    * Pound per cubic curadrado (lbcm2)
    * Pound per cubic foot (lbft2)
    * Dyn per cubic centimeter (dyncm2)

    >>> P=Pressure(760, "mmHg")
    >>> print P.bar, P.atm, P.psi, P.kgcm2g
    1.01325 1.0 14.6959487755 0.0
    """
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
                '1/mbar', '1/psi', '1/psi g', '1/atm', u'1/kg/cm²',
                u'kg/cm² g', 'mmH2O', 'cmH2O', 'mH2O', 'inH2O', 'ftH2O',
                'mmHg', 'cmHg', 'inHg', 'ftHg', u'lb/cm²', u'lb/ft²', u'dyn/cm²']
    __units__ = ['Pa', 'hPa', 'kPa', 'MPa', 'bar', 'barg', 'mbar', 'psi',
                 'psig', 'atm', 'kgcm2', 'kgcm2g', 'mmH2O', 'cmH2O', 'mH2O',
                 'inH2O', 'ftH2O', 'mmHg', 'cmHg', 'inHg', 'ftHg', 'lbcm2',
                 'lbft2', 'dyncm2']
    __units_set__ = {"altsi": "bar", "si": "Pa", "metric": "Pa",
                     "cgs": "dyncm2", "english": "psi"}

    def __init__(self, data, unit="Pa", magnitud=""):
        super(InvPressure, self).__init__(data, unit, magnitud)


class EnthalpyPressure(unidad):
    """Class that models a enthalpy per pressure measure
    Supported units:

    * Joule per kilogram pascal(JkgPa) default
    * Kilojoule per kilogram kilopascal (kJkgkPa)
    * Kilojoule per kilogram megapascal (kJkgMPa)

    >>> H=EnthalpyPressure(5, "JkgPa")
    >>> print H.JkgPa, H.kJkgMPa
    5.0 5000.0
    """
    __title__ = QApplication.translate("pychemqt", "Enthalpy per pressure")
    rates = {"JkgPa": 1.,
             "kJkgkPa": 1.,
             "kJkgMPa": k.milli}
    __text__ = ['J/kgPa', 'kJ/kgkPa', 'kJ/kgMPa']
    __units__ = ['JkgPa', 'kJkgkPa', 'kJkgMPa']
    __units_set__ = {"altsi": "kJkgkPa", "si": "JkgPa", "metric": "JkgPa",
                     "cgs": "kJkgkPa", "english": "kJkgkPa"}

    def __init__(self, data, unit="JkgPa", magnitud=""):
        super(EnthalpyPressure, self).__init__(data, unit, magnitud)


class TemperaturePressure(unidad):
    """Class that models a Temperature/Pressure measure
    Supported units:

    * Kelvin per pascal(KPa) default
    * Kelvin per kilopascal (KkPa)
    * Kelvin per megapascal (KMPa)

    >>> H=TemperaturePressure(1, "KPa")
    >>> print H.KPa, H.KkPa
    1.0 1000.0
    """
    __title__ = QApplication.translate("pychemqt", "Temperature per pressure")
    rates = {"KPa": 1.,
             "KkPa": k.milli,
             "KMPa": k.micro}
    __text__ = ['K/Pa', 'K/kPa', 'K/MPa']
    __units__ = ['KPa', 'KkPa', 'KMPa']
    __units_set__ = {"altsi": "KkPa", "si": "KPa", "metric": "KPa",
                     "cgs": "KPa", "english": "KPa"}

    def __init__(self, data, unit="KPa", magnitud=""):
        super(TemperaturePressure, self).__init__(data, unit, magnitud)


class PressureTemperature(unidad):
    """Class that models a Pressure/Temperature measure
    Supported units:

    * pascal per Kelvin (PaK) default
    * kilopascal per Kelvin (kPaK)
    * megapascal per Kelvin (MPaK)

    >>> H=PressureTemperature(1000, "PaK")
    >>> print H.PaK, H.atmK
    1000.0 0.00101325
    """
    __title__ = QApplication.translate("pychemqt", "Pressure per Temperature")
    rates = {"PaK": 1.,
             "kPaK": k.kilo,
             "barK": 1e5,
             "MPaK": k.mega,
             "atmK": 1e6/1.01325}
    __text__ = ['Pa/K', 'kPa/K', 'bar/K', 'MPa/K', "atm/K"]
    __units__ = ['PaK', 'kPaK', 'barK', 'MPaK',  "atmK"]
    __units_set__ = {"altsi": "kPaK", "si": "PaK", "metric": "PaK",
                     "cgs": "PaK", "english": "PaK"}

    def __init__(self, data, unit="PaK", magnitud=""):
        super(PressureTemperature, self).__init__(data, unit, magnitud)


class Currency(unidad):
    """Class that models a currency rate
    Supported currency codes from ISO 4217: 'pkr', 'ars', 'xcd', 'myr', 'inr',
    'hnl', 'jpy', 'czk', 'brl', 'lkr', 'sek', 'sgd', 'ttd', 'isk', 'usd',
    'aud', 'chf', 'zar', 'xpf', 'cny', 'vef', 'gtq', 'pen', 'hkd', 'hrk',
    'ang', 'xaf', 'eur', 'huf', 'vnd', 'nok', 'rub', 'pab', 'mxn', 'mad',
    'mmk', 'pln', 'php', 'jmd', 'rsd', 'cop', 'ils', 'twd', 'ghs', 'clp',
    'idr', 'krw', 'fjd', 'try', 'tnd', 'dkk', 'bsd', 'aed', 'gbp', 'nzd',
    'ron', 'thb', 'cad'

    >>> S=Currency(5, "eur")
    """
    try:
        archivo = open(conf_dir+"moneda.dat", "r")
    except:
        getrates()
        archivo = open(conf_dir+"moneda.dat", "r")
    rates = cPickle.load(archivo)
    archivo.close
    fecha = rates.pop("date")
    __title__ = QApplication.translate("pychemqt", "Currency")
    __text__ = ['$', u'€', u'£', u'¥', u'¥', u'руб', 'A$', 'R$', 'C$', 'Fr.',
                'kr', 'HK$', u'₨', u'₩', u'₨', 'RM', 'NZ$', 'S$', 'NT$',
                'R', u'฿', 'kr', 'kr', '$', u'Kč', u'Ft', u'zł', 'RON', u'Íkr',
                'kn',  'TL', 'PhP.', '$', '$', '$', u'د.ت', u'درهم', '$', 'L',
                u'฿', 'S/.', 'Rs', 'EC$', 'Dhs',  u'NAƒ', '$', 'F', 'Bs.', 'Q',
                'FCFA', u'₫',  'K', 'B$', u'дин.', u'GH₵', 'Rp', 'FJ$', u'₪']
    __units__ = ['usd', 'eur', 'gbp', 'jpy', 'cny', 'rub', 'aud', 'brl', 'cad',
                 'chf', 'dkk', 'hkd', 'inr', 'krw', 'lkr', 'myr', 'nzd', 'sgd',
                 'twd', 'zar', 'thb', 'sek', 'nok', 'mxn', 'czk', 'huf', 'pln',
                 'ron', 'isk', 'hrk', 'try', 'php', 'cop', 'ars', 'clp', 'tnd',
                 'mad', 'jmd', 'hnl', 'pab', 'pen', 'pkr', 'xcd', 'aed', 'ang',
                 'ttd', 'xpf', 'vef', 'gtq', 'xaf', 'vnd', 'mmk', 'bsd', 'rsd',
                 'ghs', 'idr', 'fjd', 'ils']
    __tooltip__ = [QApplication.translate("pychemqt", "United States dollar"),
                   QApplication.translate("pychemqt", "Euro"),
                   QApplication.translate("pychemqt", "Pound sterling"),
                   QApplication.translate("pychemqt", "Japanese yen"),
                   QApplication.translate("pychemqt", "Chinese yuan"),
                   QApplication.translate("pychemqt", "Russian rouble"),
                   QApplication.translate("pychemqt", "Australian dollar"),
                   QApplication.translate("pychemqt", "Brazilian real"),
                   QApplication.translate("pychemqt", "Canadian dollar"),
                   QApplication.translate("pychemqt", "Swiss franc"),
                   QApplication.translate("pychemqt", "Danish krone"),
                   QApplication.translate("pychemqt", "Hong Kong dollar"),
                   QApplication.translate("pychemqt", "Indian rupee"),
                   QApplication.translate("pychemqt", "South Korean won"),
                   QApplication.translate("pychemqt", "Sri Lankan rupee"),
                   QApplication.translate("pychemqt", "Malaysian ringgit"),
                   QApplication.translate("pychemqt", "New Zealand dollar"),
                   QApplication.translate("pychemqt", "Singapore dollar"),
                   QApplication.translate("pychemqt", "New Taiwan dollar"),
                   QApplication.translate("pychemqt", "South African rand"),
                   QApplication.translate("pychemqt", "Thai baht"),
                   QApplication.translate("pychemqt", "Swedish krona"),
                   QApplication.translate("pychemqt", "Norwegian krone"),
                   QApplication.translate("pychemqt", "Mexican peso"),
                   QApplication.translate("pychemqt", "Czech koruna"),
                   QApplication.translate("pychemqt", "Hungarian forint"),
                   QApplication.translate("pychemqt", "Polish złoty"),
                   QApplication.translate("pychemqt", "Romanian new leu"),
                   QApplication.translate("pychemqt", "Icelandic króna"),
                   QApplication.translate("pychemqt", "Croatian kuna"),
                   QApplication.translate("pychemqt", "Turkish lira"),
                   QApplication.translate("pychemqt", "Philippine peso"),
                   QApplication.translate("pychemqt", "Colombian peso"),
                   QApplication.translate("pychemqt", "Argentine peso"),
                   QApplication.translate("pychemqt", "Chilean peso"),
                   QApplication.translate("pychemqt", "Tunisian dinar"),
                   QApplication.translate("pychemqt", "Moroccan dirham"),
                   QApplication.translate("pychemqt", "Jamaican dollar"),
                   QApplication.translate("pychemqt", "Honduran lempira"),
                   QApplication.translate("pychemqt", "Panamanian balboa"),
                   QApplication.translate("pychemqt", "Peruvian nuevo sol"),
                   QApplication.translate("pychemqt", "Pakistani rupee"),
                   QApplication.translate("pychemqt", "East Caribbean dollar"),
                   QApplication.translate("pychemqt", "United Arab Emirates dirham"),
                   QApplication.translate("pychemqt", "Netherlands Antillean guilder"),
                   QApplication.translate("pychemqt", "Trinidad and Tobago dollar"),
                   QApplication.translate("pychemqt", "CFP franc"),
                   QApplication.translate("pychemqt", "Venezuelan bolívar fuerte"),
                   QApplication.translate("pychemqt", "Guatemalan quetzal"),
                   QApplication.translate("pychemqt", "CFA franc"),
                   QApplication.translate("pychemqt", "Vietnamese dong"),
                   QApplication.translate("pychemqt", "Myanma kyat"),
                   QApplication.translate("pychemqt", "Bahamian dollar"),
                   QApplication.translate("pychemqt", "Serbian dinar"),
                   QApplication.translate("pychemqt", "Ghanaian cedi"),
                   QApplication.translate("pychemqt", "Indonesian rupiah"),
                   QApplication.translate("pychemqt", "Fiji dollar"),
                   QApplication.translate("pychemqt", "Israeli new shekel")]
    __units_set__ = {"altsi": "usd", "si": "usd", "metric": "usd",
                     "cgs": "usd", "english": "usd"}

    def __init__(self, data=None, unit='usd', magnitud=""):
        super(Currency, self).__init__(data, unit, magnitud)

    @property
    def str(self):
        if self.code:
            return self.code
        else:
            conf = self.func(self.magnitud)
            num = self.format(conf, self.magnitud)
            txt = self.text(self.magnitud)
            return " "+txt+num


_all = unidad.__subclasses__()

_magnitudes = []
for unit in _all:
    for magnitud in unit.magnitudes():
        _magnitudes.append(magnitud+(unit, ))
_magnitudes.append(("Dimensionless",
                    QApplication.translate("pychemqt", "Dimensionless"),
                    Dimensionless))

# For get a fresh new list of magnitudes when we add some new, the list can be
# add start of lib/firstrun.py file:
# magnitudes=[]
# for magnitud, title, unit in _magnitudes:
#     magnitudes.append(magnitud)
# print magnitudes


# Run this when add some new magnitude to rebuild units_set
# unit_set={}
# for unidad in _all:
#     if unidad._magnitudes:
#         unit_set.update(unidad.__units_set__)
#     else:
#         unit_set[unidad.__name__]=unidad.__units_set__
#
# sets={}
# for set in ("altsi", "si", "metric", "cgs", "english"):
#     sets[set]=[]
#     for magnitud, titulo, unidad in _magnitudes:
#         sets[set].append(unidad.__units__.index(unit_set[magnitud][set]))
# print sets

units_set = {'cgs': [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 3, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 2, 23, 23, 6, 6, 5, 5, 3, 3, 10, 10, 3, 3, 6, 6, 6, 1, 1, 4, 3, 7, 8, 9, 2, 5, 1, 3, 0, 1, 0, 0, 1, 23, 1, 0, 0, 0],
             'si': [0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 11, 0, 0, 0, 0, 8, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             'altsi': [1, 1, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 4, 2, 2, 1, 1, 1, 1, 9, 3, 2, 2, 2, 2, 2, 0, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 4, 1, 1, 1, 0],
             'metric': [1, 1, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
             'english': [3, 3, 0, 6, 5, 5, 5, 6, 5, 4, 4, 4, 2, 0, 4, 2, 4, 3, 6, 7, 4, 4, 4, 4, 5, 7, 7, 9, 13, 7, 7, 12, 6, 14, 3, 11, 8, 12, 12, 12, 3, 3, 6, 7, 9, 10, 10, 3, 4, 2, 6, 1, 2, 1, 1, 3, 7, 1, 0, 0, 0]}


if __name__ == "__main__":
    import doctest
    doctest.testmod()

#    T=Temperature(5, "C")
#    print T.str.encode("utf-8")

    l=Length(5, "ft")
    print l.m, l.inch
#    print l.str("A")
