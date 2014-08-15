#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module with wizard project configuration
# when it creates a new pychemqt project it execute this wizard to configure:
#   -Component list
#   -Thermodynamic properties methods
#   -Transport properties methods
#   -Units prefered
#
# The wizard can be launch whanever you want using its menu entry
###############################################################################

from ConfigParser import ConfigParser
import os

from PyQt4 import QtGui

from lib.config import conf_dir
from lib.unidades import Temperature, Pressure
from tools import (UI_confComponents, UI_confTransport, UI_confThermo,
                   UI_confUnits, UI_confResolution)
from UI.widgets import Entrada_con_unidades


def auto(tmin=None, tmax=None, pmin=None, pmax=None, components=[]):
    # funcion que calcula los métodos termodinámicos más adecuados para el proyecto
    # TODO
    pass


class AutoDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(AutoDialog, self).__init__(parent)
        layout = QtGui.QGridLayout(self)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T<sub>min</sub>")),1,1)
        self.Tmin = Entrada_con_unidades(Temperature)
        layout.addWidget(self.Tmin,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T<sub>max</sub>")),2,1)
        self.Tmax = Entrada_con_unidades(Temperature)
        layout.addWidget(self.Tmax,2,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "P<sub>min</sub>")),3,1)
        self.Pmin = Entrada_con_unidades(Pressure)
        layout.addWidget(self.Pmin,3,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "P<sub>max</sub>")),4,1)
        self.Pmax = Entrada_con_unidades(Pressure)
        layout.addWidget(self.Pmax,4,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,5,1,1,2)


class Wizard(QtGui.QWizard):
    def __init__(self, config=None, parent=None):
        super(Wizard, self).__init__(parent)
        self.config=config
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Configuration wizard..."))
        self.setOptions(QtGui.QWizard.ExtendedWatermarkPixmap|QtGui.QWizard.IndependentPages|QtGui.QWizard.HaveCustomButton1)

        botonAuto=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Auto"))
        botonAuto.setToolTip(QtGui.QApplication.translate("pychemqt", "Choose good values from project components and conditions"))
        self.setButton(QtGui.QWizard.CustomButton1, botonAuto)
#        layout=QtCore.QtList()
#        layout << QtGui.QWizard.CustomButton1 << QtGui.QWizard.Stretch << QtGui.QWizard.BackButton << QtGui.QWizard.NextButton << QtGui.QWizard.FinishButton
#        self.setButtonLayout(layout)
        self.customButtonClicked.connect(self.auto)

        page1_welcome=QtGui.QWizardPage()
        page1_welcome.setTitle(QtGui.QApplication.translate("pychemqt", "Welcome"))
        page1_welcome.setSubTitle(QtGui.QApplication.translate("pychemqt", "That's the configuration wizard of a new project from pychemqt"))
        page1_welcome.setPixmap(QtGui.QWizard.LogoPixmap, QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt_98.png"))
        page1_welcome.setPixmap(QtGui.QWizard.WatermarkPixmap, QtGui.QPixmap(os.environ["pychemqt"]+"/images/logo_2.jpg"))
        lyt = QtGui.QVBoxLayout(page1_welcome)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", """<html><body>
This wizard let's you configure all parameters necessary in a pychemqt's project<br>
All options will be changed later using the options in menu Edit, and this wizard<br>
can be run at any time later.<br>
These are the options you must expecific next:<br>

<ul>
 <li>Component list</li>
 <li>Thermodinamic properties</li>
 <li>Transport properties</li>
 <li>Engineering units preferred</li>
</ul>
<body></html>""")))
        self.addPage(page1_welcome)

        page2_components=QtGui.QWizardPage()
        page2_components.setTitle(QtGui.QApplication.translate("pychemqt", "Define components"))
        page2_components.setSubTitle(QtGui.QApplication.translate("pychemqt", "Add componentes from database"))
        page2_components.setPixmap(QtGui.QWizard.LogoPixmap, QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt_98.png"))
        lyt = QtGui.QVBoxLayout(page2_components)
        self.componentes=UI_confComponents.UI_confComponents_widget(config)
        self.componentes.componentChanged.connect(self.button(QtGui.QWizard.NextButton).setEnabled)

        lyt.addWidget(self.componentes)
        self.addPage(page2_components)

        page3_thermo=QtGui.QWizardPage()
        page3_thermo.setTitle(QtGui.QApplication.translate("pychemqt", "Define thermodynamics procedures"))
        page3_thermo.setSubTitle(QtGui.QApplication.translate("pychemqt", "The thermodynamics properties are the basic of pychemqt, a bad selection would be disastrous for the results"))
        page3_thermo.setPixmap(QtGui.QWizard.LogoPixmap, QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt_98.png"))
        lyt = QtGui.QVBoxLayout(page3_thermo)
        self.thermo=UI_confThermo.UI_confThermo_widget(config)
        lyt.addWidget(self.thermo)
        self.addPage(page3_thermo)

        page4_transport=QtGui.QWizardPage()
        page4_transport.setTitle(QtGui.QApplication.translate("pychemqt", "Define transport procedures"))
        page4_transport.setSubTitle(QtGui.QApplication.translate("pychemqt", "The transport properties are important too for good simulation results"))
        page4_transport.setPixmap(QtGui.QWizard.LogoPixmap, QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt_98.png"))
        lyt = QtGui.QVBoxLayout(page4_transport)
        self.transport=UI_confTransport.UI_confTransport_widget(config)
        lyt.addWidget(self.transport)
        self.addPage(page4_transport)

        page5_units=QtGui.QWizardPage()
        page5_units.setTitle(QtGui.QApplication.translate("pychemqt", "Define preferred units"))
        page5_units.setSubTitle(QtGui.QApplication.translate("pychemqt", "The preferred units are not necessary for the simulation, but a good election let you only focus in simulation"))
        page5_units.setPixmap(QtGui.QWizard.LogoPixmap, QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt_98.png"))
        lyt = QtGui.QVBoxLayout(page5_units)
        self.units=UI_confUnits.UI_confUnits_widget(config)
        lyt.addWidget(self.units)
        self.addPage(page5_units)
        self.currentIdChanged.connect(self.checkComponents)


    def checkComponents(self, id):
        if id==1:
            self.button(QtGui.QWizard.NextButton).setEnabled(len(self.componentes.indices)!=0)
        self.button(QtGui.QWizard.CustomButton1).setVisible(id==2)

    def auto(self):
        dialogo=AutoDialog()
        if dialogo.exec_():
            tmin=dialogo.Tmin.value
            tmax=dialogo.Tmax.value
            pmin=dialogo.Pmin.value
            tmax=dialogo.Pmax.value
            auto(tmin, tmax, pmin, pmax, self.componentes.indices)


    @property
    def value(self):
        config=self.componentes.value(self.config)
        config=self.thermo.value(config)
        config=self.transport.value(config)
        config=self.units.value(config)
        if not config.has_section("PFD"):
            config.add_section("PFD")
            Preferences=ConfigParser()
            Preferences.read(conf_dir+"pychemqtrc")
            config.set("PFD", "x", Preferences.get("PFD", "x"))
            config.set("PFD", "y", Preferences.get("PFD", "y"))
        return config

    @classmethod
    def default(cls):
        config=ConfigParser()
        config=UI_confComponents.UI_confComponents_widget.default(config)
        config=UI_confThermo.UI_confThermo_widget.default(config)
        config=UI_confTransport.UI_confTransport_widget.default(config)
        config=UI_confUnits.UI_confUnits_widget.default(config)
        config=UI_confResolution.UI_confResolution_widget.default(config)
        return config

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    ui=Wizard()
    ui.show()
    sys.exit(app.exec_())
