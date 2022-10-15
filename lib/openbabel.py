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


###############################################################################
# Module with the openbabel functionality, optional for chemical compound
# graphics representation from its smile code
#
# - imageFromSmile: Generate image file with the compound chemical structure
###############################################################################


from configparser import ConfigParser
import tempfile

from qt import QtCore, QtWidgets, QtSvg, QtSvgWidgets, QtGui

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
    # https://openbabel.org/docs/current/FileFormats/SVG_2D_depiction.html

    # + u: no element-specific atom coloring
    # - U: can define alternate color to internally-specified
    # + b: background color
    # + C: Do not draw terminal C (and H) explicitly
    # + a: draw all carbon atoms
    # - d: do not display molecule name
    # + s: use asymmetric double bonds
    # + t: use thicker lines
    # - e: embed molecule as CML
    # - p <num>: px Scale to bond length(single mol only)
    # - P <num>:px Single mol in defined size image
    # - c <num>: number of columns in table
    # - cols <num>: number of columns in table(not displayed in GUI)
    # - r <num>: number of rows in table
    # - rows <num>: number of rows in table(not displayed in GUI)
    # - N <num>:max number objects to be output
    # - l: draw grid lines
    # + i: add index to each atom
    # - j: do not embed javascript
    # - x: omit XML declaration (not displayed in GUI)
    # - A: display aliases, if present

    opt = {}

    opt["B"] = conf.get("Openbabel", 'BondColor')
    alpha = conf.getint("Openbabel", "BackColorAlpha")
    if alpha:
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

    if conf.getboolean("Openbabel", 'AtomIndex'):
        opt["i"] = None

    mol = readstring("smi", smile)
    mol.write("svg", filename=fname, overwrite=True, opt=opt)

    # Edit the generated file to edit the rect viewBox to fill all the image
    with open(fname, 'r') as file:
        data = file.read()

    # Add too transparence to the background color
    opacity = ' fill-opacity="%0.2f"' % (alpha/255)
    data = data.replace(
        '<rect x="0" y="0" width="100" height="100"',
        '<rect x="0" y="0" width="200" height="200"' + opacity)

    with open(fname, 'w') as file:
        file.write(data)


class ConfBabel(QtWidgets.QDialog):
    """Openbabel image generation configuration options"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Bond color:")), 1, 1)
        self.BondColor = ColorSelector()
        self.BondColor.valueChanged.connect(self.updateImage)
        layout.addWidget(self.BondColor, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Background color:")), 2, 1)
        self.BackColor = ColorSelector(isAlpha=True)
        self.BackColor.valueChanged.connect(self.updateImage)
        layout.addWidget(self.BackColor, 2, 2)
        self.checkColor = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Heteroatom in color"))
        self.checkColor.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkColor, 3, 1, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1, 1, 2)

        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Atom details"))
        layout.addWidget(group, 5, 1, 1, 2)
        lyt = QtWidgets.QVBoxLayout(group)
        self.radioAll = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Show all atoms"))
        self.radioAll.clicked.connect(self.updateImage)
        lyt.addWidget(self.radioAll)
        self.radioEnd = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Show only terminal atoms"))
        self.radioEnd.clicked.connect(self.updateImage)
        lyt.addWidget(self.radioEnd)
        self.radioNone = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Do not show atoms"))
        self.radioNone.clicked.connect(self.updateImage)
        lyt.addWidget(self.radioNone)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 1, 1, 2)
        self.checkTighBond = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Thicker bond lines"))
        self.checkTighBond.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkTighBond, 7, 1, 1, 2)
        self.checkAsym = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate(
                "pychemqt", "Asymetric double bond"))
        self.checkAsym.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkAsym, 8, 1, 1, 2)
        self.checkIndex = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Show atoms index"))
        self.checkIndex.stateChanged.connect(self.updateImage)
        layout.addWidget(self.checkIndex, 9, 1, 1, 2)

        self.example = QtSvgWidgets.QSvgWidget()
        layout.addWidget(self.example, 13, 1, 1, 2,
                         QtCore.Qt.AlignmentFlag.AlignCenter)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 14, 1, 1, 4)

        if config and config.has_section("Openbabel"):
            self.BondColor.setColor(config.get("Openbabel", 'BondColor'))

            alpha = config.getint("Openbabel", 'BackColorAlpha')
            self.BackColor.setColor(
                config.get("Openbabel", 'BackgroundColor'), alpha)
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
            self.checkIndex.setChecked(
                config.getboolean("Openbabel", 'AtomIndex'))

            self.updateImage(config)

    def updateImage(self, config=None):
        if not isinstance(config, ConfigParser):
            config = ConfigParser()
            config = self.value(config)

        imageFile = tempfile.NamedTemporaryFile("w", suffix=".svg")
        smile = "C=CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        imageFromSmile(smile, imageFile.name, config)

        renderer = QtSvg.QSvgRenderer(imageFile.name)
        self.example.load(imageFile.name)
        self.example.renderer().setViewBox(QtCore.QRect(0, 0, 200, 200))
        self.example.setFixedSize(renderer.defaultSize())

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("Openbabel"):
            config.add_section("Openbabel")
        config.set("Openbabel", "BondColor", self.BondColor.color.name())
        config.set("Openbabel", "BackgroundColor",
                   self.BackColor.color.name(QtGui.QColor.NameFormat.HexRgb))
        config.set("Openbabel", "BackColorAlpha",
                   str(self.BackColor.color.alpha()))
        config.set("Openbabel", "AtomsColor", str(self.checkColor.isChecked()))
        config.set("Openbabel", "AtomsAll", str(self.radioAll.isChecked()))
        config.set("Openbabel", "AtomsEnd", str(self.radioEnd.isChecked()))
        config.set("Openbabel", "AtomsNone", str(self.radioNone.isChecked()))
        config.set("Openbabel", "TighBond",
                   str(self.checkTighBond.isChecked()))
        config.set("Openbabel", "AsymetricDouble",
                   str(self.checkAsym.isChecked()))
        config.set("Openbabel", "AtomIndex", str(self.checkIndex.isChecked()))
        return config
