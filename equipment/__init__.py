#!/usr/bin/python
# -*- coding: utf-8 -*-

###Modulo de equipos

#Inicializa importando las interfaces gráficas de los equipos

#flow
import UI_divider
import UI_valve
import UI_mixer
import UI_compressor
import UI_turbine
import UI_pump
import UI_pipe

#Operaciones
import UI_flash
import UI_columnFUG
import UI_heatExchanger
import UI_hairpin
import UI_shellTube
import UI_fireHeater

#solids
import UI_ciclon
import UI_gravityChamber
import UI_baghouse
import UI_electricPrecipitator
import UI_dryer

#Tools
import UI_spreadsheet

#No funcionales
import UI_centrifuge
import UI_crystallizer
import UI_filter
import UI_grinder
import UI_screen
import UI_solidWasher
import UI_vacuum
import UI_scrubber
import UI_tank
import UI_tower
import UI_reactor


UI_equipments=[UI_divider, UI_valve, UI_mixer, UI_pump, UI_compressor, UI_turbine, UI_pipe, UI_flash, UI_columnFUG, UI_heatExchanger, UI_shellTube, UI_hairpin, UI_fireHeater, UI_ciclon, UI_gravityChamber, UI_baghouse, UI_electricPrecipitator, UI_dryer, UI_scrubber, UI_spreadsheet]
#, UI_tower, UI_reactor, UI_centrifuge, UI_grinder, UI_solidWasher, UI_vacuum, ]

equipments=[ui.UI_equipment.Equipment.__class__ for ui in UI_equipments]


