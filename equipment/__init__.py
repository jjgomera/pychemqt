#!/usr/bin/python3
# -*- coding: utf-8 -*-

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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''



###Modulo de equipos

#Inicializa importando las interfaces gráficas de los equipos

# flow
from equipment import UI_divider
from equipment import UI_valve
from equipment import UI_mixer
from equipment import UI_compressor
from equipment import UI_turbine
from equipment import UI_pump
from equipment import UI_pipe

# Operaciones
from equipment import UI_flash
from equipment import UI_columnFUG
from equipment import UI_heatExchanger
from equipment import UI_hairpin
from equipment import UI_shellTube
from equipment import UI_fireHeater

# solids
from equipment import UI_ciclon
from equipment import UI_gravityChamber
from equipment import UI_baghouse
from equipment import UI_electricPrecipitator
from equipment import UI_dryer
from equipment import UI_scrubber

# Tools
from equipment import UI_spreadsheet

# No funcionales
from equipment import UI_centrifuge
from equipment import UI_crystallizer
from equipment import UI_filter
from equipment import UI_grinder
from equipment import UI_screen
from equipment import UI_solidWasher
from equipment import UI_vacuum
from equipment import UI_tank
from equipment import UI_tower
from equipment import UI_reactor


UI_equipments = [UI_divider, UI_valve, UI_mixer, UI_pump, UI_compressor,
                 UI_turbine, UI_pipe, UI_flash, UI_columnFUG, UI_heatExchanger,
                 UI_shellTube, UI_hairpin, UI_fireHeater, UI_ciclon,
                 UI_gravityChamber, UI_baghouse, UI_electricPrecipitator,
                 UI_dryer, UI_scrubber, UI_spreadsheet, UI_reactor]
# UI_tower, UI_reactor, UI_centrifuge, UI_grinder, UI_solidWasher, UI_vacuum, ]

equipments = [ui.UI_equipment.Equipment.__class__ for ui in UI_equipments]

__all__ = [equip.__name__ for equip in equipments]
__all__.sort()


# Import equipment class
from equipment.flux import Divider, Mixer, Valve
from equipment.pump import Pump
from equipment.compressor import Compressor, Turbine
from equipment.pipe import Pipe
from equipment.distillation import Flash, ColumnFUG
from equipment.heatExchanger import (Heat_Exchanger, Shell_Tube, Hairpin,
                                     Fired_Heater)
from equipment.gas_solid import (Ciclon, GravityChamber, Baghouse,
                                 ElectricPrecipitator)
from equipment.gas_solid_liquid import Dryer, Scrubber
from equipment.spreadsheet import Spreadsheet
from equipment.reactor import Reactor

# To get a list of equipment available to add to lib/firstrun.py file:
# equipos=[equipment.__name__ for equipment in equipments]
# print equipos
