#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Module for project definition (pdf of equipment, configuration and many more)
###############################################################################

import os
from configparser import ConfigParser

# from pygraph.classes.graph import graph
# from pygraph.algorithms.cycles import find_cycle
# try:
   # from pygraph.readwrite.dot import write
# except:
   # from pygraph.readwrite.markup import write

from lib.config import conf_dir
from lib.corriente import Corriente
from equipment import equipments
from equipment.flux import Mixer


class Project(object):
    MAGIC_NUMBER = 0x3051E
    FILE_VERSION = 10

    def __init__(self, items={}, streams={}, config=None):
        """
        items: diccionario con los equipos
        stream: diccionario con las corrientes
        """
        self.items = items
        self.out = {}
        if not config:
            config = ConfigParser()
            config.read(conf_dir+"pychemqtrc")
        self.config = config
        self.streams = streams
        # self.graph = self.calGraph()

        self.downToStream = {}
#        import gv
#        for item in items:
#           print gv.tailof(item)

    def __bool__(self):
        return True

#    @property
#    def equipCount(self):
#        return len(self.)
    @property
    def streamCount(self):
        return len(self.streams)

    # def calGraph(self):
        # grafo = graph()
        # for item in self.items:
            # grafo.add_node(item, attrs=[("splines",  "")])
        # for key, stream in self.streams.items():
            # grafo.add_edge((stream[0], stream[1]), attrs=[("splines",  "")])
        # return grafo

    def getObject(self, id):
        if id[0] in ["e", "i", "o"]:
            return self.items[id]
        else:
            return self.getStream(int(id[1:]))

    def setPFD(self, streams):
        self.streams = streams

    def setItems(self, items):
        self.items = items

    def setStreams(self, streams):
        self.streams = streams

    def setConfig(self, config):
        self.config = config

    def addItem(self, id, obj):
        if id not in self.items:
            self.items[id] = obj
            # self.graph.add_node(id)

    def setItem(self, id, obj):
        self.items["e%i" % id] = obj
        self.run("e%i" % id)

    def getItem(self, id):
        return self.items["e%i" % id]

    def setInput(self, id, obj):
        self.items["i%i" % id] = obj
        self.run("i%i" % id)

    def getInput(self, id):
        return self.items["i%i" % id]

    def setOutput(self, id, obj):
        self.items["o%i" % id] = obj

    def getOutput(self, id):
        return self.items["o%i" % id]

    def addStream(self, id, up, down, obj=Corriente(), ind_up=0, ind_down=0):
        if up[0] == "i":
            obj = self.items[up]

        stream = (up, down, ind_up, ind_down, obj)
        if id not in list(self.streams.keys()):
            self.streams[id] = stream
            # self.graph.add_edge((up, down))

        if down[0] == "e":
            eq = self.items[down]
            if isinstance(eq, Mixer):
                kwargs = {"entrada": obj, "id_entrada": ind_down}
            else:
                kwargs = {eq.kwargsInput[ind_down]: obj}
            eq(**kwargs)

    def setStream(self, id, obj):
        stream = self.streams[id]
        self.streams[id] = stream[0:4]+(obj, )
        self.run("s%i" % id)

    def getStream(self, id):
        return self.streams[id][-1]

    def getDownToStream(self, id):
        up, down, ind_up, ind_down, obj = self.streams[id]
        if down[0] == "e":
            return self.getItem(int(down[1]))
        else:
            return obj

    def getDownToEquip(self, str):
        lista = []
        for key, value in self.streams.items():
            if value[0] == str:
                lista.append((key, value))
        return lista

    def run(self, name):
        """Ejecuta el projecto de forma recursiva hasta que encuentra un equipo no resuelto"""
        tipo = name[0]
        ind = int(name[1])

        if tipo == "i":
            obj = self.getInput(ind)
            if obj.status:
                key, (up, down, ind_up, ind_down, oldobj) = self.getDownToEquip(name)[0]
                self.setStream(key, obj)
                self.run("s%i" % key)

        elif tipo == "e":
            obj = self.getItem(ind)
            if obj.status:
                for key, (up, down, ind_up, ind_down, oldobj) in self.getDownToEquip(name):
                    self.setStream(key, obj.salida[ind_up])
                    self.run("s%i" % key)

        elif tipo == "s":
            up, down, ind_up, ind_down, stream = self.streams[ind]
            if stream.status:
                if down[0] == "e":
                    equip = self.items[down]
                    if isinstance(equip, Mixer):
                        kwargs = {"entrada": stream, "id_entrada": ind_down}
                    else:
                        kwargs = {equip.kwargsInput[ind_down]: stream}
                    equip(**kwargs)
                    self.run(down)
                elif down[0] == "o":
                    self.setOutput(int(down[1]), stream)

    def writeToJSON(self, data):
        """Write the project to a dictionary to save to file in json format"""
        # Write configuration
        config = {}
        for section in self.config.sections():
            section_dict = {}
            for option in self.config.options(section):
                section_dict[option] = self.config.get(section, option)
            config[section] = section_dict
        data["config"] = config

        # write equipments
        equipment = {}
        for key, equip in self.items.items():
            if key[0] == "e":
                eq = {}
                eq["id"] = equipments.index(equip.__class__)
                equip.writeToJSON(eq)
                equipment[key] = eq
        data["equipment"] = equipment

        # write streams
        streams = {}
        for id, item in self.streams.items():
            stream = {}
            stream["up"] = item[0]
            stream["down"] = item[1]
            stream["ind_up"] = item[2]
            stream["ind_down"] = item[3]
            item[4].writeToJSON(stream)
            streams[id] = stream
        data["stream"] = streams

    def readFromJSON(self, data, huella=True):
        """Read project from stream
        huella: boolean to save project file to pychemqt_temporal"""
        # read configuration
        config = ConfigParser()
        for section, options in data["config"].items():
            config.add_section(section)
            for option, value in options.items():
                config.set(section, option, value)

        self.setConfig(config)
        if not huella:
            os.rename(conf_dir+"pychemqtrc_temporal", conf_dir+"pychemqtrc_temporal_bak")
        config.write(open(conf_dir+"pychemqtrc_temporal", "w"))

        # read equipments
        items = {}
        for id, equip in data["equipment"].items():
            if id[0] == "e":
                index = equip["id"]
                eq = equipments[index]()
                eq.readFromJSON(equip)
            else:
                eq = None
            items[id] = eq
        self.setItems(items)

        # read streams
        streams = {}
        for id, stream in data["stream"].items():
            id = int(id)
            up = stream["up"]
            down = stream["down"]
            ind_up = stream["ind_up"]
            ind_down = stream["ind_down"]
            obj = Corriente()
            obj.readFromJSON(stream)
            streams[id] = (up, down, ind_up, ind_down, obj)
            if huella:
                if down[0] == "e":
                    equip = self.items[down]
                    if isinstance(equip, Mixer):
                        kwargs = {"entrada": obj, "id_entrada": ind_down}
                        equip.cleanOldValues(**kwargs)
                    else:
                        kwargs = {equip.kwargsInput[ind_down]: obj}
                        equip.kwargs.update(kwargs)
                if up[0] == "e":
                    equip = self.items[up]
                    # Equipment with variable output streams must be corrected
                    while len(equip.salida) <= ind_up:
                        equip.salida.append(None)
                    equip.salida[ind_up] = obj
        self.setStreams(streams)

        if not huella:
            os.rename(conf_dir+"pychemqtrc_temporal_bak", conf_dir+"pychemqtrc_temporal")

    # def printer(self):
        # # Draw as PNG
        # dot = write(self.graph)
        # gvv = gv.readstring(dot)
        # gv.layout(gvv, 'dot')
        # gv.render(gvv, 'png', 'project.png')

    # def cycle(self):
        # cicle = find_cycle(self.graph)
        # if cicle:
            # return cicle
        # else:
            # return None

    # def hasCycle(self):
        # """Detect cycle in project"""
        # cicle = self.cycle()
        # return bool(cicle)


if __name__ == '__main__':

    from .corriente import Corriente

#    project=Project()
#    project.addItem("i1", Corriente())
#    project.addItem("i2", Corriente())
#    Cambiador=heatExchanger.Hairpin()
#    project.addItem("e1", Cambiador)
#    project.addStream(1, "i1", "e1", ind_down=0)
#    project.addStream(2, "i2", "e1", ind_down=1)
#    project.addItem("o1", Corriente())
#    project.addStream(3, "e1", "o1", ind_up=0)
#    project.addItem("o2", Corriente())
#    project.addStream(4, "e1", "o2", ind_up=1)
#    caliente=Corriente(T=140+273.15, P=361540., caudalMasico=1.36, ids=[62], fraccionMolar=[1.])
#    project.setInput(1, caliente)
#    fria=Corriente(T=20+273.15, P=101325., caudalMasico=5000/3600., ids=[62], fraccionMolar=[1.])
#    project.setInput(2, fria)
#    Cambiador(modo=1,
#                      DiTube=0.0525, DeTube=0.0603, DeeTube=0.0779, kTube=54, rTube=0.0459994e-3,
#                      annulliFouling= 0.000352, tubeFouling=0.000176, LTube=2.5)
#    project.setItem(1, Cambiador)
#    print project.getOutput(1),  project.getStream(2)
#    eq=project.getItem(1)
#    print eq.kwargs
#    print eq.status, eq.msg
#    print "Project has cycle: ", project.hasCycle()

#    bomba(entrada=entrada, rendimiento=0.75, deltaP=2)
#
#    items={"i1": entrada, "e1": bomba}
#    streams=[("i1", "e1", 0, 0), ("e1", "o1", 0, 0)]
#
#    project=Project(items=items, streams=streams)
#    project.printer()
#    print dir(gv)
#    print gv.nextin(project.graph)

    project = Project()
    project.addItem("i1", Corriente())
    project.addItem("i2", Corriente())
    mezclador = Mixer()
    project.addItem("e1", mezclador)
    project.addStream(1, "i1", "e1", ind_down=0)
    project.addStream(2, "i2", "e1", ind_down=1)
    project.addItem("o1", Corriente())
    project.addStream(3, "e1", "o1")
