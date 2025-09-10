#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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


Module with the openbabel functionality, optional module for
graphics representation of chemical compound from its smile code

  * :func:`imageFromSmile`: Generate image file with expanded formulae
  * :class:`ConfBabel`: Openbabel image generation configuration options

'''

import os
import tempfile
from configparser import ConfigParser

from tools.qt import QtCore, QtGui, QtWidgets

try:
    from openbabel.pybel import readstring
except ImportError:
    pass

from UI.widgets import ColorSelector


def imageFromSmile(smile, fname, conf):
    """Generate image file with the compound chemical structure

    Parameters
    ----------
    smile : str
        Smile code for compound
    fname : str
        Filename to save the image
    conf : configparser.ConfigParser
        Configuration parameters
    """

    # Get dict parameters for exported image configuration
    # Parameters available from openbabel documentation
    # https://openbabel.org/docs/current/FileFormats/PNG_2D_depiction.html

    # - p <num>: px Image size, default 300
    # - w <pixels>: image width (or from image size)
    # - h <pixels>: image height (or from image size)
    # - c <num>: number of columns in table
    # - r <num>: number of rows in table
    # - N <num>:max number objects to be output
    # + u: no element-specific atom coloring
    # - U: can define alternate color to internally-specified
    # + b: background color
    # + C: Do not draw terminal C (and H) explicitly
    # + a: draw all carbon atoms
    # - d: do not display molecule name
    # - m: do not add margins to the image
    # + s: use asymmetric double bonds
    # + t: use thicker lines
    # - A: display aliases, if present
    # - O <format ID>: Format of embedded text. For example, molfile or smi.
    #                  If there is no parameter, input format is used.
    # - y <additional chunk ID>: Write to a chunk with specified ID

    opt = {}
    opt["m"] = None

    opt["B"] = conf.get("Openbabel", 'BondColor')
    alpha = conf.getboolean("Openbabel", "BackColorTransparent")
    if alpha:
        opt["b"] = "none"
    else:
        opt["b"] = conf.get("Openbabel", 'BackgroundColor')

    if not conf.getboolean("Openbabel", 'AtomsColor'):
        opt["u"] = None

    if conf.getboolean("Openbabel", 'AtomsAll'):
        opt["a"] = None
    elif conf.getboolean("Openbabel", 'AtomsNone'):
        opt["C"] = None

    if conf.getboolean("Openbabel", 'TighBond'):
        opt["t"] = None

    if conf.getboolean("Openbabel", 'AsymetricDouble'):
        opt["s"] = None

    opt["x"] = None

    mol = readstring("smi", smile)
    mol.write("_png2", filename=fname, overwrite=True, opt=opt)


class ConfBabel(QtWidgets.QDialog):
    """Openbabel image generation configuration options"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(self.tr("Bond color:")), 1, 1)
        self.BondColor = ColorSelector()
        self.BondColor.valueChanged.connect(self.updateImage)
        layout.addWidget(self.BondColor, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Background color:")), 2, 1)
        self.BackColor = ColorSelector()
        self.BackColor.valueChanged.connect(self.updateImage)
        layout.addWidget(self.BackColor, 2, 2)
        self.BackColorTransparent = QtWidgets.QCheckBox(
            self.tr("Transparent background color"))
        self.BackColorTransparent.stateChanged.connect(self.updateImage)
        layout.addWidget(self.BackColorTransparent, 3, 1, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1, 1, 2)

        group = QtWidgets.QGroupBox(self.tr("Atom details"))
        layout.addWidget(group, 5, 1, 1, 2)
        lyt = QtWidgets.QVBoxLayout(group)
        self.radioAll = QtWidgets.QRadioButton(self.tr("Show all atoms"))
        self.radioAll.clicked.connect(self.updateImage)
        lyt.addWidget(self.radioAll)
        self.radioEnd = QtWidgets.QRadioButton(
            self.tr("Show only terminal atoms"))
        self.radioEnd.clicked.connect(self.updateImage)
        lyt.addWidget(self.radioEnd)
        self.radioNone = QtWidgets.QRadioButton(self.tr("Do not show atoms"))
        self.radioNone.clicked.connect(self.updateImage)
        lyt.addWidget(self.radioNone)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 1, 1, 2)
        self.checkColor = QtWidgets.QCheckBox(self.tr("Heteroatom in color"))
        self.checkColor.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkColor, 7, 1, 1, 2)
        self.checkTighBond = QtWidgets.QCheckBox(self.tr("Thicker bond lines"))
        self.checkTighBond.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkTighBond, 8, 1, 1, 2)
        self.checkAsym = QtWidgets.QCheckBox(self.tr("Asymetric double bond"))
        self.checkAsym.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkAsym, 9, 1, 1, 2)

        if os.environ["openbabel"] == "False":
            self.example = QtWidgets.QLabel(
                self.tr("Openbabel library don´t found"))
            self.example.setStyleSheet("color: red")
        else:
            self.example = QtWidgets.QLabel()
        layout.addWidget(self.example, 13, 1, 1, 2,
                         QtCore.Qt.AlignmentFlag.AlignCenter)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 14, 1, 1, 4)

        if config and config.has_section("Openbabel"):
            self.BondColor.setColor(config.get("Openbabel", 'BondColor'))

            self.BackColorTransparent.setChecked(
                config.getboolean("Openbabel", 'BackColorTransparent'))
            self.BackColor.setColor(
                config.get("Openbabel", 'BackgroundColor'))
            self.checkColor.setChecked(
                config.getboolean("Openbabel", 'AtomsColor'))
            self.radioAll.setChecked(
                config.getboolean("Openbabel", 'AtomsAll'))
            self.radioEnd.setChecked(
                config.getboolean("Openbabel", 'AtomsEnd'))
            self.radioNone.setChecked(
                config.getboolean("Openbabel", 'AtomsNone'))
            self.checkTighBond.setChecked(
                config.getboolean("Openbabel", 'TighBond'))
            self.checkAsym.setChecked(
                config.getboolean("Openbabel", 'AsymetricDouble'))

            self.updateImage(config)

    def updateImage(self, config=None):
        """Update image shown when changing some of the parameters"""
        if os.environ["openbabel"] == "False":
            return

        if not isinstance(config, ConfigParser):
            config = ConfigParser()
            config = self.value(config)

        with tempfile.NamedTemporaryFile("w", suffix=".png") as imageFile:
            smile = "C=CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
            imageFromSmile(smile, imageFile.name, config)

            self.example.setPixmap(QtGui.QPixmap(imageFile.name))

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("Openbabel"):
            config.add_section("Openbabel")
        config.set("Openbabel", "BondColor", self.BondColor.color.name())
        config.set("Openbabel", "BackgroundColor",
                   self.BackColor.color.name(QtGui.QColor.NameFormat.HexRgb))
        config.set("Openbabel", "BackColorTransparent",
                   str(self.BackColorTransparent.isChecked()))
        config.set("Openbabel", "AtomsColor", str(self.checkColor.isChecked()))
        config.set("Openbabel", "AtomsAll", str(self.radioAll.isChecked()))
        config.set("Openbabel", "AtomsEnd", str(self.radioEnd.isChecked()))
        config.set("Openbabel", "AtomsNone", str(self.radioNone.isChecked()))
        config.set("Openbabel", "TighBond",
                   str(self.checkTighBond.isChecked()))
        config.set("Openbabel", "AsymetricDouble",
                   str(self.checkAsym.isChecked()))
        return config
