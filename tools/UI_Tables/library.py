#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Library function for plugin
#   - getMethod: Return the thermo method name to use
#   - getClassFluid: Return the thermo class to calculate
#   - calcPoint: Calculate point state and check state in P-T range of eq
#   - get_propiedades: Get the properties to show in tables
#   - _getData: Get values of properties in fluid
###############################################################################


from ast import literal_eval
from configparser import ConfigParser

from numpy import nan

from lib import mEoS, coolProp, refProp, config, unidades
from lib.thermo import ThermoAdvanced


N_PROP = len(ThermoAdvanced.properties())
KEYS = ThermoAdvanced.propertiesKey()
UNITS = ThermoAdvanced.propertiesUnit()


def getMethod(pref=None):
    """Return the thermo method name to use"""
    if pref is None:
        pref = ConfigParser()
        pref.read(config.conf_dir + "pychemqtrc")

        if pref.getboolean("MEOS", 'coolprop') and \
                pref.getboolean("MEOS", 'refprop'):
            txt = "refprop"
        elif pref.getboolean("MEOS", 'coolprop'):
            txt = "coolprop"
        else:
            txt = "meos"
    else:
        txt = pref["method"]
    return txt


def getClassFluid(method, fluid):
    """Return the thermo class to calculate
    Really return the base instance to add kwargs to calculate"""

    if method == "refprop":
        # RefProp case, the base instance with the ids kwargs to define the
        # defined compount
        index = mEoS.__all__[fluid].id
        fluid = refProp.RefProp(ids=[index])
        fluid._fixed()

    elif method == "coolprop":
        # CoolProp case, the base instance with the ids kwargs to define the
        # defined compount
        index = mEoS.__all__[fluid].id
        fluid = coolProp.CoolProp(ids=[index])

    else:
        # MEOS case, the instance of specified mEoS subclass
        fluid = mEoS.__all__[fluid]()

    return fluid


def getLimit(fluid, conf):
    method = getMethod()
    if method == "meos":
        if isinstance(conf, dict):
            option = conf
        else:
            option = {}
            option["eq"] = conf.getint("MEoS", "eq")
            option["visco"] = conf.getint("MEoS", "visco")
            option["thermal"] = conf.getint("MEoS", "thermal")
        Tmin = fluid.eq[option["eq"]]["Tmin"]
        Tmax = fluid.eq[option["eq"]]["Tmax"]
        Pmin = fluid._new(T=fluid.eq[option["eq"]]["Tmin"], x=1, **option).P
        Pmax = fluid.eq[option["eq"]]["Pmax"]*1000
    elif method == "coolprop":
        Tmin = fluid.eq["Tmin"]
        Tmax = fluid.eq["Tmax"]
        Pmin = fluid.eq["Pmin"]
        Pmax = fluid.eq["Pmax"]
    elif method == "refprop":
        import refprop
        refprop.setup("def", fluid.name)
        limit = refprop.limitx([1], t=-1)
        # Using the tiple point temperature the Tmin value returned here can be
        # lower if define en thermal or viscosity procedures
        Tmin = fluid.Tt
        try:
            Pmin = fluid(T=fluid.Tt, x=1).P
        except:
            Pmin = 100
        Tmax = limit["tmax"]
        Pmax = limit["pmax"]*1000

    return Tmin, Tmax, Pmin, Pmax


def calcPoint(fluid, conf, **kwargs):
    """Procedure to calculate point state and check state in P-T range of eq"""
    method = getMethod()
    Tmin, Tmax, Pmin, Pmax = getLimit(fluid, conf)

    if "T" in kwargs:
        if kwargs["T"] < Tmin or kwargs["T"] > Tmax:
            return None
    if "P" in kwargs:
        if kwargs["P"] < Pmin-1 or kwargs["P"] > Pmax+1:
            return None

    fluido = fluid._new(**kwargs)

    if fluido.status not in [1, 3]:
        return None

    # Discard any point below the melting line, in solid state
    if method == "meos":
        if fluido._melting and fluido._melting["Tmin"] <= fluido.T \
                <= fluido._melting["Tmax"]:
            Pmel = fluido._Melting_Pressure(fluido.T)
            Pmax = min(Pmax, Pmel)

    # Discard any point out of limit of equation
    if fluido.P < Pmin-1 or fluido.P > Pmax+1 or fluido.T < Tmin \
            or fluido.T > Tmax:
        return None

    return fluido


def get_propiedades(conf):
    """Procedure to get the properties to show in tables
    Input:
        conf: configparser instance with mainwindow preferences
    eanOutput:
        array with properties, key and units
    """
    booleanos = conf.get("MEoS", "properties")
    order = conf.get("MEoS", "propertiesOrder")
    if isinstance(booleanos, str):
        booleanos = literal_eval(booleanos)
    if isinstance(order, str):
        order = literal_eval(order)

    propiedades = []
    keys = []
    units = []
    for indice, boolean in zip(order, booleanos):
        if boolean:
            name, key, unit = ThermoAdvanced.properties()[indice]
            propiedades.append(name)
            keys.append(key)
            units.append(unit)
    return propiedades, keys, units


def _getData(fluid, keys, phase=True, unit=None):
    """Procedure to get values of properties in fluid
    Input:
        fluid: fluid instance to get values
        keys: array with desired parameter to get
        phase: boolean to get the properties values for both phases
        unit: unidades subclass
    """
    fila = []
    for i, key in enumerate(keys):
        if not key:
            continue
        p = fluid.__getattribute__(key)

        if isinstance(p, list):
            p = p[0]

        if isinstance(p, str):
            txt = p
        # elif isinstance(p, list):
            # txt = repr(p)
        else:
            if unit and unit[i]:
                txt = p.__getattribute__(unit[i])
            else:
                txt = p.config()
        fila.append(txt)

        # Add two phases properties is requested
        if phase and key in ThermoAdvanced.propertiesPhase():
            # Liquid
            p = fluid.Liquido.__getattribute__(key)
            if isinstance(p, str):
                txt = p
            elif isinstance(p, unidades.unidad):
                if unit and unit[i]:
                    txt = p.__getattribute__(unit[i])
                else:
                    txt = p.config()
            else:
                txt = p
            fila.append(txt)
            # Gas
            p = fluid.Gas.__getattribute__(key)
            if isinstance(p, str):
                txt = p
            elif isinstance(p, unidades.unidad):
                if unit and unit[i]:
                    txt = p.__getattribute__(unit[i])
                else:
                    txt = p.config()
            else:
                txt = p
            fila.append(txt)
    return fila


def saveProperties(fluids):
    """Save all available properties of a list of fluids"""
    dat = {}
    for key in ThermoAdvanced.propertiesKey():
        prop = []
        for fluid in fluids:
            if isinstance(fluid, list):
                # Special case when the fluids is a 2D array
                prop_i = []
                for fluidi in fluid:
                    if fluidi is not None and fluidi.status in [1, 3]:
                        p = fluidi.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        prop_i.append(p)
                    else:
                        prop_i.append(nan)

                prop.append(prop_i)

            else:
                if fluid is not None and fluid.status in [1, 3]:
                    p = fluid.__getattribute__(key)
                    if key in ["fi", "f"]:
                        p = p[0]
                    prop.append(p)
        dat[key] = prop
    return dat
