#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, random, sys

#TODO: Borrar cuando no haga falta debug
#os.environ["pychemqt"]="/home/jjgomera/pychemqt/"
#os.environ["freesteam"]="True"
#os.environ["oasa"]="True"
#os.environ["Elemental"]="True"
#os.environ["CoolProp"]="True"
#os.environ["refprop"]="True"
#os.environ["ezodf"]="True"
#os.environ["openpyxl"]="True"
#os.environ["xlwt"]="True"

from ConfigParser import ConfigParser

from PyQt4 import QtCore, QtGui

from sql import databank


conf_dir = os.path.expanduser('~') + os.sep+".pychemqt"+os.sep

Preferences=ConfigParser()
Preferences.read(conf_dir+"pychemqtrc")


def getComponents(solidos=False, config=None, name=True):
    if not config:
        config=getMainWindowConfig()
    if solidos:
        indices=config.get("Components", "Solids")
        if not isinstance(indices, list):
            indices=eval(indices)
    else:
        indices=config.get("Components", "Components")
        if not isinstance(indices, list):
            indices=eval(indices)
    if name:
        nombres=[]
        M=[]
        for componente in indices:
            databank.execute("select nombre, peso_molecular from compuestos where id == %s" %str(componente))
            texto=databank.fetchone()
            nombres.append(texto[0])
            M.append(texto[1])
        return indices, nombres, M
    else:
        return indices

def getMainWindowConfig():
    #FIXME: Al cargar un nuevo proyecto desde cero se carga con la configuracion del proyecto anterior.
    widget=QtGui.QApplication.activeWindow()
    config=None
    if not widget:
        config=ConfigParser()
        config.read(conf_dir+"pychemqtrc_temporal")
    if isinstance(widget, QtGui.QMainWindow) and widget.__class__.__name__=="UI_pychemqt":
        config=widget.currentConfig
    else:
        lista=QtGui.QApplication.topLevelWidgets()
        for widget in lista:
            if isinstance(widget, QtGui.QMainWindow) and widget.__class__.__name__=="UI_pychemqt":
                config=widget.currentConfig
                break
    if not config:
        config=ConfigParser()
        config.read(conf_dir+"pychemqtrc_temporal")
    return config

#indices, nombres, M=getComponents()
#solidos, nombreSolidos, MSolidos=getComponents(solidos=True)


def representacion(float, format=0, total=0, decimales=4, exp=False, tol=4, signo=False, thousand=False):
    """Función que expresa un valor de tipo float en forma de string
    float: numero a representar
    format: tipo de modo
        0   -   fixed point
        1   -   Significant figures
        2   -   Engineering format
    total: numero total de digitos
    decimales: numero de decimales
    exp: boolean que indica si se usa notacion exponencial para numeros grandes y pequeños
    tol: potencia por encima de la cual se usa notacion exponencial
    signo: mostrar signo positivo
    thousand: usa separador para miles
    """
    if type(float) is str:
        return float

    if signo:
        start="{:+"
    else:
        start="{: "

    if thousand:
        coma=",."
    else:
        coma="."

    if exp:
        if -10**tol > float or -10**-tol < float < 10**-tol or float > 10**tol:
            format=2

    if format==1:
        string=start+"{}{:d}g".format(coma, decimales)+"}"
    elif format==2:
        string=start+"{:d}{}{:d}e".format(total, coma, decimales)+"}"
    else:
        string=start+"{:d}{}{:d}f".format(total, coma, decimales)+"}"

    return string.format(float)


def colors(int):
    """Función que devuelve una lista de colores con el número de elementos indicados"""
    r = lambda: random.randint(0,255)
    colores=[]
    for i in range(int):
        colores.append(('#%02X%02X%02X' % (r(),r(),r())))
    return colores


class Entity(object):
    """Clase general que define la funcionalidad comun básica:
    -clear object
    -definition of boolean characteristic of object
    -input used for definition
    -save/load from file
    -properties for report, tooltip

    Clases que heredan de él:
    -Corriente, Mezcla, Solids
    -equipment

    """
    _bool=False
    kwargs_forbidden=["entrada"]
    notas=""
    notasPlain=""

    def clear(self):
        self.kwargs=self.__class__.kwargs
        self.__dict__.clear()
        self._bool=False

    def __nonzero__(self):
        return self._bool

    def show(self):
        for key in sorted(self.__dict__):
            print key, ": ", self.__dict__[key]

    def setNotas(self, html, txt):
        self.notas=html
        self.notasPlain=txt
        if html:
            self._bool=True

    @property
    def numInputs(self):
        """numero de input disponibles en el kwargs"""
        count=0
        for key, value in self.kwargs.iteritems():
            if key not in self.kwargs_forbidden and value:
                count+=1
        return count

#Read-Write to file
    def writeToStream(self, stream):
        """Guarda las propiedades usadas para su definición en kwargs"""
        stream.writeInt32(self.numInputs)
        for key, value in self.kwargs.iteritems():
            if key not in self.kwargs_forbidden and value:
                stream.writeString(key)
                if isinstance(value, float):
                    stream.writeFloat(value)
                elif isinstance(value, int):
                    stream.writeInt32(value)
                elif isinstance(value, str):
                    stream.writeString(value)
                elif isinstance(value, list):
                    self.writeListtoStream(stream, key, value)

    def writeListtoStream(self, stream, key, value):
        """Personalizar en el caso de equipos con listas complejas"""
        stream.writeInt32(len(value))
        for val in value:
            stream.writeFloat(val)

    def readFromStream(self, stream):
        for i in range(stream.readInt32()):
            key=stream.readString()
            if isinstance(self.kwargs[key], float):
                valor=stream.readFloat()
            elif isinstance(self.kwargs[key], int):
                valor=stream.readInt32()
            elif isinstance(self.kwargs[key], str):
                valor=stream.readString()
            elif isinstance(self.kwargs[key], list):
                valor=self.readListFromStream(stream, key)
            self.kwargs[key]=valor
        self.__call__()

    def readListFromStream(self, stream, key):
        """Personalizar en el caso de equipos con listas complejas"""
        valor=[]
        for i in range(stream.readInt32()):
            valor.append(stream.readFloat())
        return valor

#Properties
    @classmethod
    def propertiesNames(cls):
        """Lista de los nombres de las propiedades disponibles para cada entity"""
        return []

    def properties(self):
        lista=[]
        for name, attr, unit in self.propertiesNames():
            prop=self._prop(attr)
            if isinstance(prop, list):
                if unit==str:
                    lista.append((prop, name))
                elif unit==None:
                    lista.append((["" for f in prop], name))
                else:
                    lista.append(([f.str for f in prop], name))
            elif unit == str:
                lista.append((prop, name, 0))
            elif unit == int:
                lista.append((str(prop), name, 1))
            else:
                lista.append((prop.str, name, 1))
        return lista

    def _prop(self, attr):
        if attr=="className":
            prop=self.__class__.__name__
        elif attr=="notasPlain":
            prop=self.notasPlain
        elif type(attr)==tuple:
            prop=self.__getattribute__(attr[0])[self.kwargs[attr[1]]]
        elif attr in self.__dict__:
            prop=self.__getattribute__(attr)
        elif attr in self.kwargs and self.kwargs[attr]:
            prop=self.kwargs[attr]
        return prop

    def propertiesListTitle(self, index):
        """Define los titulos para los popup de listas"""
        pass

    def propertiesTitle(self):
        return [prop[0] for prop in self.propertiesNames()]
    def propertiesAttribute(self):
        return [prop[1] for prop in self.propertiesNames()]
    def propertiesUnit(self):
        return [prop[2] for prop in self.propertiesNames()]

    def popup(self, preferences, exception=[]):
        """A partir de las preferencias dadas, devuelve el texto a ser mostrado en un popup"""
        txt=[]
        propiedades=self.properties()
        for i in eval(preferences.get("TooltipEntity", self.__class__.__name__)):
            if isinstance(propiedades[i][0], list):
                txt.append((propiedades[i][1], "", 0))
                title=self.propertiesListTitle(i)
                for name, value in zip(title, propiedades[i][0]):
                    txt.append(("%s\t%s" %(name, value), "", 1))
            else:
                txt.append(propiedades[i])
        return txt



class fluid(dict):
    """Clase que personaliza una lista con propiedades nulas para las fases nulas en las corrientes"""
    h=0
    s=0
    cp=0
    cv=0
    cp_cv=0
    cp0=0


if __name__ == "__main__":
    import math
    print representacion(math.pi, decimales=6, tol=1)
#    print repr(Configuracion("Density", "DenGas").text())
    print representacion("3232326262")
