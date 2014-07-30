#!/usr/bin/python
# -*- coding: utf-8 -*-

from ConfigParser import ConfigParser

from PyQt4 import QtCore, QtGui

from UI.widgets import Entrada_con_unidades
from lib.config import conf_dir

class UI_confResolution_widget(QtGui.QWidget):
    def __init__(self, config=None, parent=None):
        self.standards=[(600, 400), (640, 480), (720, 400), (800, 600), (832, 624), (1024, 768), (1152, 864), (1280, 1024), (1700, 1250), (1900, 1425), (2400, 1800), (4000, 3000)]
        super(UI_confResolution_widget, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Use default resolution:")),0,0)
        self.standard = QtGui.QComboBox()
        self.standard.addItem("")
        for resolucion in self.standards:
            self.standard.addItem("%ix%i" %resolucion)
        self.standard.currentIndexChanged.connect(self.changeResolution)
        layout.addWidget(self.standard,0,1)
        
        self.checkCustom=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use Custom resolution"))
        layout.addWidget(self.checkCustom,1,0,1,2)
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Width:"))
        label.setIndent(50)
        layout.addWidget(label,2,0)
        self.x=Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        layout.addWidget(self.x,2,1)
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Height:"))
        label.setIndent(50)
        layout.addWidget(label,3,0)
        self.y=Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        layout.addWidget(self.y,3,1)
        
        self.checkCustom.toggled.connect(self.x.setEnabled)
        self.checkCustom.toggled.connect(self.y.setEnabled)
        
        if config and config.has_section("PFD"):
            x=config.getint("PFD","x")
            y=config.getint("PFD","y")
            self.x.setValue(x)
            self.y.setValue(y)
            if (x, y) in self.standards:
                self.standard.setCurrentIndex(self.standards.index((x, y))+1)
                self.checkCustom.setChecked(False)
                self.x.setEnabled(False)
                self.y.setEnabled(False)
            else:
                self.standard.setCurrentIndex(0)
                self.checkCustom.setChecked(True)

    def changeResolution(self):
        x, y=self.standard.currentText().split("x")
        self.x.setValue(int(x))
        self.y.setValue(int(y))
        
    def value(self, config):
        if not config.has_section("PFD"):
            config.add_section("PFD")
        config.set("PFD", "x", str(self.x.value))
        config.set("PFD", "y", str(self.y.value))
        return config
    
    @classmethod
    def default(cls, config):
        config.add_section("PFD")
        Preferences=ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        config.set("PFD", "x", Preferences.get("PFD", "x"))
        config.set("PFD", "y", Preferences.get("PFD", "y"))
        return config


class Dialog(QtGui.QDialog):
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Define PFD resolution"))
        layout = QtGui.QVBoxLayout(self)
        self.datos=UI_confResolution_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        config=self.datos.value(config)
        return config


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec_())
