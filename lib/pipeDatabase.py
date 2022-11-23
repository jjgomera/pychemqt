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
# Pipe database
# Module with pipedatabase support
# Description translation
###############################################################################

import os
import sqlite3

from tools.qt import QtWidgets


# Standard pipe database

# Material de la tubería:
 # 0 - Nombre
 # 1 - Tipo
 # 2 - Rugosidad, mm
 # 3 - Diametro nominal, mm
 # 4 - Diametro nominal, pulgadas
 # 5 - Espesor de la pared, mm
 # 6 - Diametro externo, mm
 # 7 - Peso, kg/m


path = os.path.join(os.environ["pychemqt"], 'dat', 'pipeDatabase.db')
databank = sqlite3.connect(path).cursor()
databank.execute("select * from Materials")
CATALOG = databank.fetchall()
CATALOG_TRANSLATE = {
    "Cast Iron (Asphalt Dipped)": QtWidgets.QApplication.translate("pychemqt", 'Cast Iron (Asphalt Dipped)'),
    "Cast Iron": QtWidgets.QApplication.translate("pychemqt", "Cast Iron"),
    "Copper": QtWidgets.QApplication.translate("pychemqt", "Copper"),
    "Drains, Waste, Vents": QtWidgets.QApplication.translate("pychemqt", "Drains, Waste, Vents"),
    "Refrigeration Service": QtWidgets.QApplication.translate("pychemqt", "Refrigeration Service"),
    "Copper Tube": QtWidgets.QApplication.translate("pychemqt", "Copper Tube"),
    "PVC (Iron pipe size)": QtWidgets.QApplication.translate("pychemqt", "PVC (Iron pipe size"),
    "PVC (Sewer pipe)": QtWidgets.QApplication.translate("pychemqt", "PVC (Sewer pipe"),
    "Steel (ANSI)": QtWidgets.QApplication.translate("pychemqt", "Steel (ANSI)"),
    "Steel Galvanised (ANSI)": QtWidgets.QApplication.translate("pychemqt", "Steel Galvanised (ANSI)"),
    "Stainless Steel (ANSI)": QtWidgets.QApplication.translate("pychemqt", "Stainless Steel (ANSI)")}


# Pipe fitting K values

# Accesorios de la Tuberia:
#  0 - Diámetro interno, mm
#  1 - Indice imagen
#  2 - Tipo
#  3 - Diámetro interno, pulgadas
#  4 - Descripción
#  5 - K

path = os.environ["pychemqt"]+"dat/pipeDatabase.db"
databank.execute("select * from Fitting")
FITTING = databank.fetchall()
FITTING_DESC = {
    "SB": QtWidgets.QApplication.translate("pychemqt", 'Standard Bend'),
    "LB": QtWidgets.QApplication.translate("pychemqt", 'Long Bend'),
    "PB": QtWidgets.QApplication.translate("pychemqt", 'Pipe Bend'),
    'MB90': QtWidgets.QApplication.translate("pychemqt", "Mitre bend") + " 90\xb0",
    "E45": QtWidgets.QApplication.translate("pychemqt", 'Elbow') + " 90\xb0",
    'MB45': QtWidgets.QApplication.translate("pychemqt", "Mitre bend")+" 45\xb0",
    'RB': QtWidgets.QApplication.translate("pychemqt", "Return bend"),
    'TT': QtWidgets.QApplication.translate("pychemqt", "Through Tee"),
    'BT': QtWidgets.QApplication.translate("pychemqt", "Branch Tee"),
    'St': QtWidgets.QApplication.translate("pychemqt", "Strainer"),
    'Open': QtWidgets.QApplication.translate("pychemqt", "Open pipe Exit"),
    'ExitCon': QtWidgets.QApplication.translate("pychemqt", "Pipe Exit to Container"),
    'EntSharp': QtWidgets.QApplication.translate("pychemqt", "Pipe Entry Sharp Edged"),
    'EntProj': QtWidgets.QApplication.translate("pychemqt", "Pipe Entry Projecting"),
    'Plug': QtWidgets.QApplication.translate("pychemqt", "Plug Valve Straightway"),
    'Globe': QtWidgets.QApplication.translate("pychemqt", "Globe Valve"),
    'Gate': QtWidgets.QApplication.translate("pychemqt", "Gate Valve"),
    'Bfly': QtWidgets.QApplication.translate("pychemqt", "Butterfly Valve"),
    'BallRB': QtWidgets.QApplication.translate("pychemqt", "Ball Valve Reduced Bore"),
    'BallFB': QtWidgets.QApplication.translate("pychemqt", "Ball Valve Full Bore"),
    'Foot': QtWidgets.QApplication.translate("pychemqt", "Foot Valve with Strainer"),
    'Hinged': QtWidgets.QApplication.translate("pychemqt", "Hinged Foot Valve with Strainer"),
    'ChWaf': QtWidgets.QApplication.translate("pychemqt", "Wafer Check Valve"),
    'ChSw': QtWidgets.QApplication.translate("pychemqt", "Check Swing Valve"),
    'TiltCh': QtWidgets.QApplication.translate("pychemqt", "Tilting Disk Check"),
    'LiftCh': QtWidgets.QApplication.translate("pychemqt", "Lift Check Valve"),
    'Angle': QtWidgets.QApplication.translate("pychemqt", "Globe Valve Angled"),
    'AngCh': QtWidgets.QApplication.translate("pychemqt", "Lift Check Angled")}


# Thermal Conductivity of material, W/mK a 0ºC
thermalConductivity = {
    "Titanium": 22.5,
    "Copper": 401.,
    "Ziconium": 22.6,
}

# Density of material, kg/m3 a 0ºC
density = {
    "Titanium": 4510.,
    "Copper": 8960.,
    "Ziconium": 6520.,
}


finnedTube_database = [
#0 - HPT
#1 - Catalogo
#2 - Plain section OD (mm)
#3 - Plain section average wall (mm)
#4 - Wall under fin average (mm)
#5 - Wall under fin minimum (mm)
#6 - Nominal root diameter (mm)
#7 - Fin section ID (mm)
#8 - Outsite Area Ao (m2/m)
#9 - Inside Area Ai (m2/m)
#10 - Area Ratio Ao/Ai
#11 - I.D. Cross Section Area (cm2)
#12 - mean fin thickness (mm)
#13 - Fin height (mm)
    ["HPT", "304028", 15.875, 1.245, .711, .635, 14.249, 12.827, .125, .040, 3.114, 1.292, .279, .813],
    ["HPT", "304035", 15.875, 1.473, .889, .787, 14.249, 12.471,  .125, .039, 3.186, 1.222 , .279, .813],
    ["HPT", "304042", 15.875, 1.651, 1.067, .940, 14.249, 12.116,  .125, .038, 3.288, 1.153 , .279, .813],
    ["HPT", "304049", 15.875, 1.829, 1.245, 1.118, 14.249, 11.760,  .125, .037, 3.397, 1.086 , .279, .813],
    ["HPT", "304065", 15.875, 2.108, 1.651, 1.473, 14.249, 10.947,  .125, .034, 3.637, .941 , .279, .813],
    ["HPT", "305028", 19.050, 1.245, .711, .635,  17.424, 16.002, .152,  .050,  3.030, 2.011 , .279, .813],
    ["HPT", "305035", 19.050, 1.473, .889, .787,  17.424, 15.646, .152,  .049,  3.106, 1.923 , .279, .813],
    ["HPT", "305042", 19.050, 1.651, 1.067, .940,  17.424, 15.291, .152,  .048,  3.165, 1.836 , .279, .813],
    ["HPT", "305049", 19.050, 1.829, 1.245, 1.118,  17.424, 14.935, .152,  .047,  3.247, 1.752 , .279, .813],
    ["HPT", "305065", 19.050, 2.108, 1.651, 1.473,  17.424, 14.122, .152,  .045,  3.425, 1.566 , .279, .813],
    ["HPT", "305083", 19.050, 2.769, 2.108, 1.880,  17.424, 13.208, .152,  .041,  3.676, 1.370 , .279, .813],
    ["HPT", "306028", 22.225, 1.245, .711, .635,  20.599, 19.177, .179,  .060,  2.965, 2.888 , .279, .813],
    ["HPT", "306035", 22.225, 1.473, .889, .787,  20.599, 18.821, .179,  .059,  3.026, 2.782 , .279, .813],
    ["HPT", "306042", 22.225, 1.651, 1.067, .940,  20.599, 18.466, .179,  .058,  3.089, 2.678 , .279, .813],
    ["HPT", "306049", 22.225, 1.829, 1.245, 1.118,  20.599, 18.110, .179,  .057,  3.139, 2.576 , .279, .813],
    ["HPT", "306065", 22.225, 2.108, 1.651, 1.473,  20.599, 17.297, .179,  .054,  3.298, 2.350 , .279, .813],
    ["HPT", "306083", 22.225, 2.769, 2.108, 1.880,  20.599, 16.383, .179,  .052,  3.473, 2.108 , .279, .813],
    ["HPT", "307035", 25.400, 1.473, .889, .787,  23.774, 21.996, .205,  .069,  2.956, 3.800 , .279, .813],
    ["HPT", "307042", 25.400, 1.651, 1.067, .940,  23.774, 21.641, .205,  .068,  3.009, 3.678 , .279, .813],
    ["HPT", "307049", 25.400, 1.829, 1.245, 1.118,  23.774, 21.285, .205,  .067,  3.064, 3.558 , .279, .813],
    ["HPT", "307065", 25.400, 2.108, 1.651, 1.473,  23.774, 20.472, .205,  .064,  3.180, 3.292 , .279, .813],
    ["HPT", "307083", 25.400, 2.769, 2.108, 1.880,  23.774, 19.558, .205,  .062,  3.322, 3.004 , .279, .813],
    ["HPT", "284028", 15.875, 1.245, .711, .635, 14.097, 12.675, .126, .040, 3.153,  1.262 , .305, .889],
    ["HPT", "284035", 15.875, 1.473, .889, .787, 14.097, 12.319, .126, .039, 3.252,  1.192 , .305, .889],
    ["HPT", "284042", 15.875, 1.651, 1.067, .940, 14.097, 11.963, .126, .037, 3.358,  1.124 , .305, .889],
    ["HPT", "284049", 15.875, 1.829, 1.245, 1.118,  14.097, 11.608, .126, .037, 3.442,  1.058 , .305, .889],
    ["HPT", "284065", 15.875, 2.108, 1.651, 1.473,  14.097, 10.795, .126, .034, 3.721,  .915 , .305, .889],
    ["HPT", "285028", 19.050, 1.245, .711, .635, 17.272, 15.850, .153, .050, 3.074,  1.973 , .305, .889],
    ["HPT", "285035", 19.050, 1.473, .889, .787, 17.272, 15.494, .153, .049, 3.131,  1.885 , .305, .889],
    ["HPT", "285042", 19.050, 1.651, 1.067, .940, 17.272, 15.138, .153, .048, 3.212,  1.800 , .305, .889],
    ["HPT", "285049", 19.050, 1.829, 1.245, 1.118,  17.272, 14.783, .153, .046, 3.296,  1.716 , .305, .889],
    ["HPT", "285065", 19.050, 2.108, 1.651, 1.473,  17.272, 13.970, .153, .044, 3.479,  1.533 , .305, .889],
    ["HPT", "285083", 19.050, 2.769, 2.108, 1.880,  17.272, 13.056, .153, .041, 3.711,  1.339 , .305, .889],
    ["HPT", "286028", 22.225, 1.245, .711, .635, 20.447, 19.025, .179, .060, 3.000,  2.843 , .305, .889],
    ["HPT", "286035", 22.225, 1.473, .889, .787, 20.447, 18.669, .179, .059, 3.063,  2.737 , .305, .889],
    ["HPT", "286042", 22.225, 1.651, 1.067, .940, 20.447, 18.313, .179, .058, 3.111,  2.634 , .305, .889],
    ["HPT", "286049", 22.225, 1.829, 1.245, 1.118,  20.447, 17.958, .179, .056, 3.178,  2.533 , .305, .889],
    ["HPT", "286065", 22.225, 2.108, 1.651, 1.473,  20.447, 17.145, .179, .054, 3.322,  2.309 , .305, .889],
    ["HPT", "286083", 22.225, 2.769, 2.108, 1.880,  20.447, 16.231, .179, .051, 3.521,  2.069 , .305, .889],
    ["HPT", "287035", 25.400, 1.473, .889, .787, 23.622, 21.844, .206, .069, 3.004,  3.748 , .305, .889],
    ["HPT", "287042", 25.400, 1.651, 1.067, .940, 23.622, 21.488, .206, .067, 3.059,  3.627 , .305, .889],
    ["HPT", "287049", 25.400, 1.829, 1.245, 1.118,  23.622, 21.133, .206, .066, 3.101,  3.508 , .305, .889],
    ["HPT", "287065", 25.400, 2.108, 1.651, 1.473,  23.622, 20.320, .206, .064, 3.234,  3.243 , .305, .889],
    ["HPT", "287083", 25.400, 2.769, 2.108, 1.880,  23.622, 19.406, .206, .061, 3.380,  2.958 , .305, .889],
    ["HPT", "264028", 15.875,  1.245, .711, .635, 13.386,  11.963,  .150, .037,  4.000, 1.124 , .330, 1.245],
    ["HPT", "264035", 15.875,  1.473, .889, .787, 13.386,  11.608,  .150, .037,  4.100, 1.058 , .330, 1.245],
    ["HPT", "264042", 15.875,  1.651, 1.067,  .940, 13.386,  11.252,  .150, .035,  4.241, .994 , .330, 1.245],
    ["HPT", "264049", 15.875,  1.829, 1.245,  1.118, 13.386,  10.897,  .150, .034,  4.393, .933 , .330, 1.245],
    ["HPT", "264065", 15.875,  2.108, 1.651,  1.473, 13.386,  10.084,  .150, .032,  4.731, .799 , .330, 1.245],
    ["HPT", "265028", 19.050,  1.245, .711, .635, 16.561,  15.138,  .182, .048,  3.821, 1.800 , .330, 1.245],
    ["HPT", "265035", 19.050,  1.473, .889, .787, 16.561,  14.783,  .182, .046,  3.921, 1.716 , .330, 1.245],
    ["HPT", "265042", 19.050,  1.651, 1.067,  .940, 16.561,  14.427,  .182, .045,  4.000, 1.635 , .330, 1.245],
    ["HPT", "265049", 19.050,  1.829, 1.245,  1.118, 16.561,  14.072,  .182, .044,  4.110, 1.555 , .330, 1.245],
    ["HPT", "265065", 19.050,  2.108, 1.651,  1.473, 16.561,  13.259,  .182, .042,  4.350, 1.381 , .330, 1.245],
    ["HPT", "265083", 19.050,  2.769, 2.108,  1.880, 16.561,  12.344,  .182, .039,  4.693, 1.197 , .330, 1.245],
    ["HPT", "266028", 22.225,  1.245, .711, .635, 19.736,  18.313,  .215, .058,  3.725, 2.634 , .330, 1.245],
    ["HPT", "266035", 22.225,  1.473, .889, .787, 19.736,  17.958,  .215, .056,  3.805, 2.533 , .330, 1.245],
    ["HPT", "266042", 22.225,  1.651, 1.067,  .940, 19.736,  17.602,  .215, .055,  3.890, 2.433 , .330, 1.245],
    ["HPT", "266049", 22.225,  1.829, 1.245,  1.118, 19.736,  17.247,  .215, .054,  3.955, 2.336 , .330, 1.245],
    ["HPT", "266065", 22.225,  2.108, 1.651,  1.473, 19.736,  16.434,  .215, .052,  4.166, 2.121 , .330, 1.245],
    ["HPT", "266083", 22.225,  2.769, 2.108,  1.880, 19.736,  15.519,  .215, .049,  4.400, 1.892 , .330, 1.245],
    ["HPT", "267035", 25.400,  1.473, .889, .787, 22.911,  21.133,  .247, .066,  3.720, 3.508 , .330, 1.245],
    ["HPT", "267042", 25.400,  1.651, 1.067,  .940, 22.911,  20.777,  .247, .065,  3.790, 3.391 , .330, 1.245],
    ["HPT", "267049", 25.400,  1.829, 1.245,  1.118, 22.911,  20.422,  .247, .064,  3.862, 3.275 , .330, 1.245],
    ["HPT", "267065", 25.400,  2.108, 1.651,  1.473, 22.911,  19.609,  .247, .062,  4.015, 3.020 , .330, 1.245],
    ["HPT", "267083", 25.400,  2.769, 2.108,  1.880, 22.911,  18.694,  .247, .059,  4.202, 2.745 , .330, 1.245],
    ["HPT", "364028", 15.875, 1.245, .711, .635, 14.554,  13.132, .125, .041, 3.044,  1.354 , .305, .660],
    ["HPT", "364035", 15.875, 1.473, .889, .787, 14.554,  12.776, .125, .040, 3.114,  1.282 , .305, .660],
    ["HPT", "364042", 15.875, 1.651, 1.067, .940, 14.554,  12.421, .125, .039, 3.211,  1.212 , .305, .660],
    ["HPT", "364049", 15.875, 1.829, 1.245, 1.118,  14.554,  12.065, .125, .038, 3.315,  1.143 , .305, .660],
    ["HPT", "364065", 15.875, 2.108, 1.651, 1.473,  14.554,  11.252, .125, .035, 3.543,  .994 , .305, .660],
    ["HPT", "365028", 19.050, 1.245, .711, .635, 17.729,  16.307, .152, .051, 2.976,  2.088 , .305, .660],
    ["HPT", "365035", 19.050, 1.473, .889, .787, 17.729,  15.951, .152, .050, 3.049,  1.998 , .305, .660],
    ["HPT", "365042", 19.050, 1.651, 1.067, .940, 17.729,  15.596, .152, .049, 3.106,  1.910 , .305, .660],
    ["HPT", "365049", 19.050, 1.829, 1.245, 1.118,  17.729,  15.240, .152, .048, 3.185,  1.824 , .305, .660],
    ["HPT", "365065", 19.050, 2.108, 1.651, 1.473,  17.729,  14.427, .152, .045, 3.356,  1.635 , .305, .660],
    ["HPT", "365083", 19.050, 2.769, 2.108, 1.880,  17.729,  13.513, .152, .042, 3.597,  1.434 , .305, .660],
    ["HPT", "366028", 22.225, 1.245, .711, .635, 20.904,  19.482, .179, .061, 2.920,  2.981 , .305, .660],
    ["HPT", "366035", 22.225, 1.473, .889, .787, 20.904,  19.126, .179, .060, 2.980,  2.873 , .305, .660],
    ["HPT", "366042", 22.225, 1.651, 1.067, .940, 20.904,  18.771, .179, .059, 3.041,  2.767 , .305, .660],
    ["HPT", "366049", 22.225, 1.829, 1.245, 1.118,  20.904,  18.415, .179, .058, 3.089,  2.663 , .305, .660],
    ["HPT", "366065", 22.225, 2.108, 1.651, 1.473,  20.904,  17.602, .179, .055, 3.243,  2.433 , .305, .660],
    ["HPT", "366083", 22.225, 2.769, 2.108, 1.880,  20.904,  16.688, .179, .052, 3.413,  2.187 , .305, .660],
    ["HPT", "367035", 25.400, 1.473, .889, .787, 24.079,  22.301, .205, .070, 2.917,  3.906 , .305, .660],
    ["HPT", "367042", 25.400, 1.651, 1.067, .940, 24.079,  21.946, .205, .069, 2.969,  3.783 , .305, .660],
    ["HPT", "367049", 25.400, 1.829, 1.245, 1.118,  24.079,  21.590, .205, .068, 3.009,  3.661 , .305, .660],
    ["HPT", "367065", 25.400, 2.108, 1.651, 1.473,  24.079,  20.777, .205, .065, 3.136,  3.391 , .305, .660],
    ["HPT", "367083", 25.400, 2.769, 2.108, 1.880,  24.079,  19.863, .205, .062, 3.273,  3.099 , .305, .660],
    ["HPT", "435023", 19.050, 0.889, .584, .508, 17.932, 16.764, .152, .053, 2.890, 2.207 , .254, .559],
    ["HPT", "365025-204610", 19.050, 1.245, .635, .559, 17.729, 16.459, .152, .062, 2.463, 2.116 , .305, .66],
    ["HPT", "365028-204610", 19.050, 1.245, .711, .635, 17.729, 16.307, .152, .061, 2.488, 2.077 , .305, .66],
    ["HPT", "285028-204610", 19.050, 1.245, .711, .635, 17.272, 15.850, .153, .059, 2.569, 1.961 , .305, .889],
    ["HPT", "285035-204610", 19.050, 1.473, .889, .787, 17.272, 15.494, .153, .058, 2.623, 1.871 , .305, .889],
    ["HPT", "265028-204610", 19.050, 1.245, .711, .635, 16.561, 15.138, .182, .057, 3.187, 1.787 , .330, 1.245],
    ["HPT", "265035-204610", 19.050, 1.473, .889, .787, 16.561, 14.783, .182, .055, 3.275, 1.703 , .330, 1.245],

#EA TRUFIN
#Nombre
#Catalogo
#OD (inch)
#OD (mm)
#wall (inch)
#wall (mm)
#OD plain (inch)
#OD plain (mm)
#wall plain (inch)
#wall plain (mm)
#fin per inch
#Finished Fin OD (inch)
#Finished Fin OD (mm)
#Min. Wall Under Fins (inch)
#Min. Wall Under Fins (mm)
#Root Diameter (inch)
#Root Diameter (mm)
    ["E/A Trufin", "68-113028", "5/8", 15.88, 0.028, 0.711, 0.583, 14.81, 0.043, 1.09, 11, 0.590, 14.99, 0.025, 0.635, 0.492, 12.50, 0.276, 0.411, 0.436, 11.07, 0.043, 1.092, 0.114, 0.035, 0.133, 0.041, 0.133, 0.041, 0.299, 0.091],
    ["E/A Trufin", "68-115030", "3/4", 19.05, 0.030, 0.762, 0.743, 18.87, 0.043, 1.09, 11, 0.753, 19.13, 0.027, 0.686, 0.657, 16.69, 0.358, 0.533, 0.597, 15.16, 0.043, 1.092, 0.156, 0.048, 0.180, 0.055, 0.180, 0.055, 0.411, 0.125]
]


if __name__ == "__main__":

    import sqlite3
    conn = sqlite3.connect('pipeDatabase2.db')
    curs = conn.cursor()
    curs.execute("""
    CREATE TABLE Fitting (
    id  INTEGER PRIMARY KEY,
    Di            FLOAT,
    type          TEXT,
    Di_nominal_in TEXT,
    K             FLOAT
    )"""
)

    for indice, material in enumerate(Accesorios_Tuberia):
        query = "INSERT INTO Fitting VALUES "
        del material[4]
        del material[1]
        material.insert(0, indice+1)
        curs.execute(query+str(tuple(material)))

    conn.commit()
    conn.close()
