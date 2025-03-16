#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=wrong-import-position, assignment-from-no-return


"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Module with global configuration of pychemqt

Variables:

  * :const:`conf_dir`: User configuration path
  * :const:`Preferences`: ConfigParser instance with pychemqt preferences
  * :const:`currentConfig`: ConfigParser instance with the configuration of
    current pychemqt project open or the last open project

Configuration tools

  * :func:`getComponents`: Get component list from project
  * :func:`getMainWindowConfig`: Return config of current project
  * :func:`setMainWindowConfig`: Update currentconfig variable
  * :class:`Entity`: General class for model object
"""

from configparser import ConfigParser
import os

from tools.qt import QtWidgets

# TODO: Delete when it isn´t necessary debug
# os.environ["pychemqt"] = "/home/jjgomera/Programacion/pychemqt/"
# os.environ["freesteam"] = "False"
# os.environ["openbabel"] = "True"
# os.environ["CoolProp"] = "True"
# os.environ["refprop"] = "False"
# os.environ["ezodf"] = "True"
# os.environ["openpyxl"] = "True"
# os.environ["xlwt"] = "True"
# os.environ["icu"] = "True"
# os.environ["reportlab"] = "True"
# os.environ["Qsci"] = "False"


from lib.sql import databank


conf_dir = os.path.expanduser('~') + os.sep + ".pychemqt" + os.sep
IMAGE_PATH = os.path.join(os.environ["pychemqt"], "images") + os.sep

Preferences = ConfigParser()
Preferences.read(conf_dir + "pychemqtrc")
# FIXME: This instance is not updated when preferences are changed

global currentConfig
currentConfig = ConfigParser()
currentConfig.read(conf_dir + "pychemqtrc_temporal")


def getComponents(solidos=False, config=None, name=True):
    """
    Procedure to get index of component in current project

    Parameters
    ----------
    solidos : bool
        Return too solids components
    config : ConfigParser
        It's possible use a custom config instance
    name : bool
        Return too the name and molecular weight of components

    Returns
    -------
    id : list
        List of index of components
    name : list
        List of name of components
    M : list
        List of molecular weight of components
    """
    if not config:
        config = getMainWindowConfig()
    if solidos:
        indices = config.get("Components", "Solids")
        if not isinstance(indices, list):
            indices = eval(indices)
    else:
        indices = config.get("Components", "Components")
        if not isinstance(indices, list):
            indices = eval(indices)

    if name:
        nombres = []
        M = []
        for id in indices:
            query = "select name, M from compuestos where id == %s" % str(id)
            databank.execute(query)
            texto = databank.fetchone()
            nombres.append(texto[0])
            M.append(texto[1])
        return indices, nombres, M
    else:
        return indices


def getMainWindowConfig():
    """Return config of current project"""
    return currentConfig


def setMainWindowConfig(config=None):
    """Set config as current project"""
    global currentConfig
    if config:
        currentConfig = config
        return
    else:
        widget = QtWidgets.QApplication.activeWindow()
        if isinstance(widget, QtWidgets.QMainWindow) and \
           widget.__class__.__name__ == "UI_pychemqt":
            currentConfig = widget.currentConfig
        else:
            lista = QtWidgets.QApplication.topLevelWidgets()
            for widget in lista:
                if isinstance(widget, QtWidgets.QMainWindow) and \
                   widget.__class__.__name__ == "UI_pychemqt":
                    currentConfig = widget.currentConfig
                    break


class Entity(object):
    """General class for model object, with basic functionality:

        * clear object
        * definition of boolean characteristic of object
        * input used for definition
        * note properties with description
        * save/load from file
        * properties for report, tooltip

    Child class include:

        * Corriente, Mezcla, Solids
        * equipment
    """
    _bool = False
    kwargs_forbidden = ["entrada"]
    notas = ""
    notasPlain = ""
    _dependence = ""

    def __init__(self, **kwargs):
        """Class constructor, copy kwargs for child class, it can be customize
        for child class to add functionality"""
        self.kwargs = self.__class__.kwargs.copy()

        # Values defined as integer Entrada_con_unidades return as float
        # and it must be corrected
        self.kwargsInteger = []
        for key, value in self.kwargs.items():
            if isinstance(value, int):
                self.kwargsInteger.append(key)

        if kwargs:
            self.__call__(**kwargs)

    def __call__(self, **kwargs):
        """Add callable functionality, so it can be possible add kwargs,
        advanced functionality can be added in subclass"""
        self._oldkwargs = self.kwargs.copy()
        self.cleanOldValues(**kwargs)
        self._bool = True
        txt = kwargs.get("notas", "")
        if txt:
            self.notas = txt
            self.notasPlain = txt

        for key in self.kwargsInteger:
            if key not in self.kwargs_forbidden:
                self.kwargs[key] = int(self.kwargs[key])

    def cleanOldValues(self, **kwargs):
        """Update kwargs with new input kwargs, defined in child class,
        here can be implemented kwarg incompatibiity input and more"""
        self.kwargs.update(kwargs)

    def clear(self):
        """Clear entity and stay as new instance"""
        self.kwargs = self.__class__.kwargs
        self.__dict__.clear()
        self._bool = False

    def __bool__(self):
        return self._bool

    def show(self):
        """General function to show entity properties as key: value text"""
        for key in sorted(self.__dict__):
            print(key, ": ", self.__dict__[key])

    def setNotas(self, html, txt):
        self.notas = html
        self.notasPlain = txt
        if html:
            self._bool = True

    @property
    def numInputs(self):
        """Input count in kwargs"""
        count = 0
        for key, value in self.kwargs.items():
            if key not in self.kwargs_forbidden and value:
                count += 1
        return count

    # Read-Write to file
    def writeToJSON(self, data):
        """Save kwargs properties of entity to file"""
        kwarg = {}
        for key, value in self.kwargs.items():
            if key not in self.kwargs_forbidden and value:
                if isinstance(value, list):
                    self.writeListtoJSON(kwarg, key, value)
                else:
                    kwarg[key] = value
        data["kwarg"] = kwarg

        # Write state
        data["status"] = self.status
        data["msg"] = self.msg
        data["bool"] = self._bool
        data["external_dependences"] = self._dependence
        data["notas"] = self.notas
        data["notasPlain"] = self.notasPlain
        if self.status:
            state = {}
            self.writeStatetoJSON(state)
            data["state"] = state
        else:
            data["state"] = {}

    def writeListtoJSON(self, data, key, value):
        """Write list to file, Customize in entities with complex list"""
        data[key] = value

    def readFromJSON(self, data):
        """Read entity from file"""
        for key, value in data["kwarg"].items():
            if isinstance(self.kwargs[key], list):
                value = self.readListFromJSON(data["kwarg"], key)
            self.kwargs[key] = value

        # Read state
        self.status = data["status"]
        self.msg = data["msg"]
        self._bool = data["bool"]
        try:
            self._dependence = data["external_dependences"]
        except KeyError:
            self._dependence = ""
        self.notas = data["notas"]
        self.notasPlain = data["notasPlain"]
        if self.status:
            self.readStatefromJSON(data["state"])


#        if run:
#            self.__call__()
#            print(self)

    def readListFromJSON(self, data, key):
        """Read list from file, customize in entities with complex list"""
        return data[key]

    def writeStatetoJSON(self, data):
        pass

    def readStatefromJSON(self, data):
        pass

    # Properties
    @classmethod
    def propertiesNames(cls):
        """List with properties availables for show in poput"""
        return []

    def properties(self):
        lista = []
        for name, attr, unit in self.propertiesNames():
            prop = self._prop(attr)
            if isinstance(prop, list):
                if unit == str:
                    lista.append((prop, name))
                elif unit is None:
                    lista.append((["" for f in prop], name))
                else:
                    lista.append(([f.str for f in prop], name))
            elif unit == str:
                lista.append((prop, name, 0))
            elif unit == int or unit == float:
                lista.append((str(prop), name, 1))
            else:
                lista.append((prop.str, name, 1))
        return lista

    def _prop(self, attr):
        if attr == "className":
            prop = self.__class__.__name__
        elif attr == "notasPlain":
            prop = self.notasPlain
        elif type(attr) == tuple:
            prop = self.__getattribute__(attr[0])[self.kwargs[attr[1]]]
        elif attr in self.__dict__:
            prop = self.__getattribute__(attr)
        elif attr in self.kwargs:
            prop = self.kwargs[attr]
        return prop

    def propertiesListTitle(self, index):
        """Define titles for list properties in popup"""

    def propertiesTitle(self):
        return [prop[0] for prop in self.propertiesNames()]

    def propertiesAttribute(self):
        return [prop[1] for prop in self.propertiesNames()]

    def propertiesUnit(self):
        return [prop[2] for prop in self.propertiesNames()]

    def popup(self, preferences, exception=[]):
        """
        Return a list with properties of entity to show in a puput
        preferences: ConfigParser instance with selected properties to show"""
        txt = []
        propiedades = self.properties()
        for i in map(int, preferences.get(
                "TooltipEntity", self.__class__.__name__).split(",")):
            if isinstance(propiedades[i][0], list):
                txt.append((propiedades[i][1], "", 0))
                title = self.propertiesListTitle(i)
                for name, value in zip(title, propiedades[i][0]):
                    txt.append(("%s\t%s" % (name, value), "", 1))
            else:
                txt.append(propiedades[i])
        return txt

    def propertiesToText(self, index=None, linesep=True, suffix="",
                         kwCheck=False, kwKey="", kwSuffix="", kwValue=""):
        """
        Return a string representation of properties for report
        index: Index of properties in propertiesList
            None: Return all properties
            array: Return the selected properties
            int: Return only the indexed property
        linesep: Boolean to add a linesep at end line
        suffix: Optional suffix text
        """
        mask = "%s-%is%ss" % ("%", self.TEXT_FORMATING_LENG, "%")
        if index is None:
            index = range(len(self.propertiesNames()))
        if isinstance(index, int):
            index = [index]

        txt = ""
        for i in index:
            title, prop, unit = self.propertiesNames()[i]
            if not kwKey:
                kwKey = prop
            value = self._prop(prop)
            if unit != str and unit != float and unit != int:
                value = value.str
            elif unit == str:
                value = " " + value
            elif unit == float or unit == int:
                value = " %s" % value.__repr__()

            txt += mask % (title, value)
            if kwCheck:
                if self.kwargs[kwKey] == kwValue:
                    txt += kwSuffix
            if suffix:
                txt += suffix
            if linesep:
                txt += os.linesep
        return txt
