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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from functools import partial
from math import pi, log, cos, atan

from tools.qt import QtCore, QtWidgets, translate

from equipment.widget.gui import ToolGui, CallableEntity
from lib.unidades import Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades


__doi__ = {
    1:
        {"autor": "El-Genk, M.S., Timothy, M.S.",
         "title": "A Review and Correlations for Convection Heat Transfer and "
                  "Pressure Losses in Toroidal and Helically Coiled Tubes",
         "ref": "Heat Transfer Eng. 38(5) (2017) 447-474",
         "doi": "10.1080/01457632.2016.1194693"},
    2:
        {"autor": "Xin, R.C., Ebadian, M.A.  ",
         "title": "The Effects of Prandtl Numbers on Local and Average "
                  "Convective Heat Transfer Characteristics in Helical Pipes",
         "ref": "J. Heat Transfer 119(3) (1997) 467-73",
         "doi": "10.1115/1.2824120."}
    # 3:
        # {"autor": "",
         # "title": "",
         # "ref": "",
         # "doi": ""}
}




# Heat Transfer coefficient correlations
@refDoc(__doi__, [2, 1])
def Nu_turbulent_XinEbadian(Re, Pr, Di, Dc):
    r"""Calculates Nusselt number for internal flow of a helical coil using the
    correlation of Xin-Ebadian (1997)

    .. math::
        Nu = 0.00619Re^{0.92} Pr^{0.4}\left[1 + 3.455\left(\frac{D_i}{D_c}
        \right)\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Laminar flow
    if Re < 5000:
        De = Re*(Di/Dc)**0.5

        # Eq 5
        Nu = (2.153+0.318*De**0.643) * Pr**0.177
    else:
        # Eq 6
        Nu = 0.00619 * Re**0.92 * Pr**0.4 * (1+3.455*Di/Dc)


class Helical(CallableEntity):
    """Helical coil tube used as anhancing heat transfer equipment.

    Parameters
    ----------
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]
    """

    status = 0
    msg = ""
    kw = {
        "methodReCritic": 0,
        "methodFrictionLaminar": 0,
        "methodFrictionTurbulent": 0,
        "methodHeatLaminar": 0,
        "methodHeatTurbulent": 0,

        "H": 0,
        "Di": 0,
        "Dc": 0}

    valueChanged = QtCore.pyqtSignal(object)
    inputChanged = QtCore.pyqtSignal(object)

    @property
    def isCalculable(self):
        """Check if all input are defined"""
        if not self.kw["Di"]:
            self.msg = translate("equipment", "undefined internal diameter")
            self.status = 0
            return False
        if not self.kw["Dc"]:
            self.msg = translate("equipment", "undefined diameter of helix")
            self.status = 0
            return False

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        """Definition of twisted tape inserts for annuli sections"""
        self.Di = self.kw["Di"]
        self.Dc = self.kw["Dc"]

        self.valueChanged.emit(self)

    @property
    def ReCritical(self):
        """Calculate critical Reynolds number to define transition of regimen
        flow from laminar to turbulent"""
        return 2000

    def Nu(self, Re, Pr):
        """Calculate nusselt number"""
        Rec = self.ReCritical

        if Re < Rec:
            # Laminar flow
            pass
        else:
            Nu = Nu_turbulent_XinEbadian(Re, Pr, self.Di, self.Dc)

        return Nu

    def f(self, Re):
        """Calculate friction factor"""
        pass


class UI_Helical(ToolGui):
    """Helical coil dialog"""

    title = translate("equipment", "Use helical coil")

    def loadUI(self):
        """Add widget"""
        self.Entity = Helical()

        lyt = self.wdg.layout()

        # label = QtWidgets.QLabel(self.tr("Tape pitch, H"))
        # label.setToolTip(self.tr("Tape pitch for twist of π radians (180º)"))
        # lyt.addWidget(label, 2, 1)
        # self.H = Entrada_con_unidades(Length)
        # self.H.valueChanged.connect(partial(self.changeParams, "H"))
        # lyt.addWidget(self.H, 2, 2)
        # label = QtWidgets.QLabel("Di")
        # label.setToolTip(self.tr("Internal diameter of annuli section"))
        # lyt.addWidget(label, 3, 1)
        # self.Di = Entrada_con_unidades(Length, "PipeDiameter")
        # self.Di.valueChanged.connect(partial(self.changeParams, "Di"))
        # lyt.addWidget(self.Di, 3, 2)
        # label = QtWidgets.QLabel("Do")
        # label.setToolTip(self.tr("External diameter of annuli section"))
        # lyt.addWidget(label, 4, 1)
        # self.Do = Entrada_con_unidades(Length, "PipeDiameter")
        # self.Do.valueChanged.connect(partial(self.changeParams, "Do"))
        # lyt.addWidget(self.Do, 4, 2)

        # self.angled = QtWidgets.QCheckBox(self.tr("Angled twisted-tape"))
        # self.angled.toggled.connect(self.setEnableOrientation)
        # lyt.addWidget(self.angled, 5, 1, 1, 2)

        # lytH = QtWidgets.QHBoxLayout()
        # self.labelOrientation = QtWidgets.QLabel(self.tr(
            # "Direction of flow relative to the tape curvature"))
        # lytH.addWidget(self.labelOrientation)
        # self.orientation = QtWidgets.QComboBox()
        # for method in TwistedTapeAnnuli.TEXT_ORIENTACION:
            # self.orientation.addItem(method)
        # self.orientation.currentIndexChanged.connect(
            # partial(self.changeParams, "orientation"))
        # lytH.addWidget(self.orientation)
        # lyt.addLayout(lytH, 7, 1, 1, 2)
        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)

    def setEnableOrientation(self, boolean):
        """Change Enable/Disable state for orientation of twisted tape"""
        self.labelOrientation.setEnabled(boolean)
        self.orientation.setEnabled(boolean)
        self.changeParams("angled", boolean)


class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Twisted-tape insert"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_Helical()
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec())
