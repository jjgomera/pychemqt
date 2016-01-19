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



###Modulo de equipos

#Inicializa importando las interfaces gráficas de los equipos

#flow
from . import UI_divider
from . import UI_valve
from . import UI_mixer
from . import UI_compressor
from . import UI_turbine
from . import UI_pump
from . import UI_pipe

#Operaciones
from . import UI_flash
from . import UI_columnFUG
from . import UI_heatExchanger
from . import UI_hairpin
from . import UI_shellTube
from . import UI_fireHeater

#solids
from . import UI_ciclon
from . import UI_gravityChamber
from . import UI_baghouse
from . import UI_electricPrecipitator
from . import UI_dryer

#Tools
from . import UI_spreadsheet

#No funcionales
from . import UI_centrifuge
from . import UI_crystallizer
from . import UI_filter
from . import UI_grinder
from . import UI_screen
from . import UI_solidWasher
from . import UI_vacuum
from . import UI_scrubber
from . import UI_tank
from . import UI_tower
from . import UI_reactor


UI_equipments=[UI_divider, UI_valve, UI_mixer, UI_pump, UI_compressor, UI_turbine, UI_pipe, UI_flash, UI_columnFUG, UI_heatExchanger, UI_shellTube, UI_hairpin, UI_fireHeater, UI_ciclon, UI_gravityChamber, UI_baghouse, UI_electricPrecipitator, UI_dryer, UI_scrubber, UI_spreadsheet, UI_reactor]
#, UI_tower, UI_reactor, UI_centrifuge, UI_grinder, UI_solidWasher, UI_vacuum, ]

equipments=[ui.UI_equipment.Equipment.__class__ for ui in UI_equipments]

# To get a list of equipment available to add to lib/firstrun.py file:
# equipos=[equipment.__name__ for equipment in equipments]
# print equipos
