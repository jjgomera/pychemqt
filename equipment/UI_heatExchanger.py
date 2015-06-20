#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Generic heat exchanger equipment dialog
###############################################################################

from functools import partial
import os

from matplotlib import image
from PyQt4 import QtGui
from scipy import arange, arctan, logspace, pi

from lib.unidades import (Temperature, DeltaT, DeltaP, Power, Area,
                          HeatTransfCoef)
from lib.plot import mpl
from UI.widgets import Entrada_con_unidades
from equipment.heatExchanger import Heat_Exchanger, Heat_ExchangerDesign
from equipment.parents import UI_equip


class UI_equipment(UI_equip):
    """Generic heat exchanger equipment edition dialog"""
    Equipment = Heat_Exchanger()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(Heat_Exchanger, entrada=False,
                                           salida=False, parent=parent)

        # Calculate tab
        lyt = QtGui.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Output temperature")), 1, 1)
        self.Tout = Entrada_con_unidades(Temperature)
        self.Tout.valueChanged.connect(partial(self.changeParams, "Tout"))
        lyt.addWidget(self.Tout, 1, 2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Temperature increase")), 2, 1)
        self.DeltaT = Entrada_con_unidades(DeltaT)
        self.DeltaT.valueChanged.connect(partial(self.changeParams, "DeltaT"))
        lyt.addWidget(self.DeltaT, 2, 2)
        lyt.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Heat Duty")), 3, 1)
        self.Heat = Entrada_con_unidades(Power)
        self.Heat.valueChanged.connect(partial(self.changeParams, "Heat"))
        lyt.addWidget(self.Heat, 3, 2)
        group = QtGui.QGroupBox(
            QtGui.QApplication.translate("pychemqt", "Heat Transfer"))
        lyt.addWidget(group, 4, 1, 1, 2)
        lyt = QtGui.QGridLayout(group)
        lyt.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Area")), 1, 1)
        self.A = Entrada_con_unidades(Area)
        self.A.valueChanged.connect(partial(self.changeParams, "A"))
        lyt.addWidget(self.A, 1, 2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Heat Transfer Coefficient")), 2, 1)
        self.U = Entrada_con_unidades(HeatTransfCoef)
        self.U.valueChanged.connect(partial(self.changeParams, "U"))
        lyt.addWidget(self.U, 2, 2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "External Temperature")), 3, 1)
        self.Text = Entrada_con_unidades(Temperature)
        self.Text.valueChanged.connect(partial(self.changeParams, "Text"))
        lyt.addWidget(self.Text, 3, 2)

        lyt.addItem(QtGui.QSpacerItem(
            10, 10, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            5, 0, 1, 3)
        lyt.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Pressure loss")), 6, 1)
        self.DeltaP = Entrada_con_unidades(DeltaP, value=0)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        lyt.addWidget(self.DeltaP, 6, 2)
        lyt.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            7, 0, 1, 3)

        group = QtGui.QGroupBox(
            QtGui.QApplication.translate("pychemqt", "Results"))
        lyt.addWidget(group, 8, 1, 1, 5)
        layout = QtGui.QGridLayout(group)
        layout.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Heat Duty")), 0, 1)
        self.HeatCalc = Entrada_con_unidades(Power, retornar=False)
        self.HeatCalc.setReadOnly(True)
        layout.addWidget(self.HeatCalc, 0, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Output Temperature")), 1, 1)
        self.ToutCalc = Entrada_con_unidades(Temperature, retornar=False)
        self.ToutCalc.setReadOnly(True)
        layout.addWidget(self.ToutCalc, 1, 2)

        lyt.addItem(QtGui.QSpacerItem(
            0, 0, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            9, 0, 1, 3)

        if equipment:
            self.setEquipment(equipment)

    def changeParams(self, parametro, valor):
        if parametro == "Tout":
            self.Heat.clear()
            self.DeltaT.clear()
        elif parametro == "DeltaT":
            self.Heat.clear()
            self.Tout.clear()
        elif parametro == "Heat":
            self.DeltaT.clear()
            self.Tout.clear()
        self.calculo(**{parametro: valor})


class Chart(QtGui.QDialog):
    """Dialog to implement general heat exchanger plot"""
    PosImage=None

    def resizeEvent(self, event):
        """Implement resizeEvent to precalculate the new position of image"""
        self.refixImage()
        QtGui.QDialog.resizeEvent(self, event)

    def refixImage(self, event=None):
        """Recalculate image position as resize/move window"""
        xmin, xmax=self.diagrama.axes2D.get_xlim()
        ymin, ymax=self.diagrama.axes2D.get_ylim()
        a=self.PosLogo[0]*xmax+(1-self.PosLogo[0])*xmin
        b=self.PosLogo[1]*ymax+(1-self.PosLogo[1])*ymin
        x, y=self.diagrama.axes2D.transData.transform_point((a, b))
        hlogo, wlogo=self.logo.get_size()
        himage, wimage=self.image.get_size()
        if self.PosLogo[0]==1:
            x-=wlogo
            xi=x-wimage-20
        else:
            xi=x+wlogo+20
        if self.PosLogo[1]==1:
            yi=y-himage-10
            y-=hlogo
        else:
            yi=y+10

        self.logo.ox=x
        self.logo.oy=y

        if self.PosImage:
            a=self.PosImage[0]*xmax+(1-self.PosImage[0])*xmin
            b=self.PosImage[1]*ymax+(1-self.PosImage[1])*ymin
            xi, yi=self.diagrama.axes2D.transData.transform_point((a, b))
            if self.PosImage[0]==1:
                xi-=wimage+10
            else:
                xi+=10
            if self.PosImage[1]==1:
                yi-=himage+10
            else:
                yi+=10

        self.image.ox=xi
        self.image.oy=yi
        self.diagrama.draw()

    def save(self):
        """Show the dialog to select name to save the image as file"""
        fname=unicode(QtGui.QFileDialog.getSaveFileName(
            self, QtGui.QApplication.translate("pychemqt", "Save chart to file"),
            "./", "Portable Network Graphics (*.png)"))
        self.diagrama.fig.savefig(fname, facecolor='#eeeeee')


class Efectividad(Chart):
    """Heat Exchanger effectiveness plot"""
    title = QtGui.QApplication.translate("pychemqt",
                                         "Heat Exchanger effectiveness")
    flujo = [
        (QtGui.QApplication.translate("pychemqt", "Counterflow"), "CF"),
        (QtGui.QApplication.translate("pychemqt", "Parallelflow"), "PF"),
        (QtGui.QApplication.translate(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (QtGui.QApplication.translate(
            "pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
        (QtGui.QApplication.translate(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (QtGui.QApplication.translate(
            "pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    mezclado = ("Cmin", "Cmax")
    PosLogo = (0, 1)
    PosImage = (1, 0)

    def __init__(self, parent=None):
        super(Efectividad, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Flow Arrangement")), 1, 1)
        self.flow=QtGui.QComboBox()
        for text, flow in self.flujo:
            self.flow.addItem(text)
        self.flow.currentIndexChanged.connect(self.plot)
        layout.addWidget(self.flow, 1, 2)
        self.mixed=QtGui.QComboBox()
        for text in self.mezclado:
            self.mixed.addItem(text)
        self.mixed.currentIndexChanged.connect(self.changeMixed)
        layout.addWidget(self.mixed, 1, 3)
        layout.addItem(QtGui.QSpacerItem(
            10, 10, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed), 1, 4)

        self.diagrama = mpl(self, width=10, height=10)
        self.image=self.diagrama.fig.figimage([[0]], 0, 0, zorder=1)
        logo=image.imread('images/pychemqt_98.png')
        self.logo=self.diagrama.fig.figimage(logo, 0, 0, zorder=1)

        self.botonGuardar=QtGui.QToolButton()
        self.botonGuardar.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")))
        self.botonGuardar.setToolTip(QtGui.QApplication.translate("pychemqt", "Save chart to file"))
        self.botonGuardar.clicked.connect(self.save)
        layout.addWidget(self.botonGuardar,5,1)
        layout.addWidget(self.diagrama,2,1,1,6)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 2, 1, 5)

        self.plot(0)
        self.showMaximized()
        self.refixImage()

    def changeMixed(self, indice):
        modelo=self.flow.currentIndex()
        self.plot(modelo)

    def plot(self, indice):
        self.diagrama.axes2D.clear()
        self.diagrama.axes2D.set_xlim(0, 6)
        self.diagrama.axes2D.set_ylim(0, 1)
        title = QtGui.QApplication.translate("pychemqt",
                                             "Heat Transfer effectiveness")
        self.diagrama.axes2D.set_title(title, size='12')
        self.diagrama.axes2D.set_xlabel("NTU", size='12')
        self.diagrama.axes2D.set_ylabel(u"ε", size='14')

        flujo=self.flujo[indice][1]
        self.mixed.setVisible(flujo=="CrFSMix")
        kw={}
        if flujo=="CrFSMix":
            kw["mixed"]=str(self.mixed.currentText())

        C=[0, 0.2, 0.4, 0.6, 0.8, 1.]

        NTU=arange(0, 6.1, 0.1)
        for  ci in C:
            e=[0]
            for N in NTU[1:]:
                e.append(Heat_ExchangerDesign.efectividad(N, ci, flujo, **kw))
            self.diagrama.plot(NTU, e, "k")

            fraccionx=(NTU[40]-NTU[30])/6
            fracciony=(e[40]-e[30])
            try:
                angle=arctan(fracciony/fraccionx)*360/2/pi
            except ZeroDivisionError:
                angle=90

            self.diagrama.axes2D.annotate(
                "C*=%0.1f" %ci, (NTU[29], e[30]), rotation=angle, size="medium",
                horizontalalignment="left", verticalalignment="bottom")

        self.diagrama.draw()

        img=image.imread('images/equation/%s.png' %flujo)
        self.image.set_data(img)
        self.refixImage()


class TemperatureEfectividad(Chart):
    title=QtGui.QApplication.translate("pychemqt", "Heat Exchanger temperature effectiveness")
    flujo=[(QtGui.QApplication.translate("pychemqt", "Counterflow"), "CF"),
                (QtGui.QApplication.translate("pychemqt", "Parallelflow"), "PF"),
                (QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
                (QtGui.QApplication.translate("pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
                (QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
                (QtGui.QApplication.translate("pychemqt", "1-2 TEMA E"), "1-2TEMAE"),
                (QtGui.QApplication.translate("pychemqt", "1-2 TEMA E, shell fluid divided"), "1-2TEMAE2"),
                (QtGui.QApplication.translate("pychemqt", "1-3 TEMA E"), "1-3TEMAE"),
                (QtGui.QApplication.translate("pychemqt", "1-4 TEMA E"), "1-4TEMAE"),
                (QtGui.QApplication.translate("pychemqt", "1-1 TEMA G"), "1-1TEMAG"),
                (QtGui.QApplication.translate("pychemqt", "1-2 TEMA G"), "1-2TEMAG"),
                (QtGui.QApplication.translate("pychemqt", "1-1 TEMA H"), "1-1TEMAH"),
                (QtGui.QApplication.translate("pychemqt", "1-2 TEMA H"), "1-2TEMAH"),
                (QtGui.QApplication.translate("pychemqt", "1-1 TEMA J"), "1-1TEMAJ"),
                (QtGui.QApplication.translate("pychemqt", "1-2 TEMA J"), "1-2TEMAJ"),
                (QtGui.QApplication.translate("pychemqt", "1-4 TEMA J"), "1-4TEMAJ"),

                ]

    mezclado=("1", "2")
    PosLogo=(0, 1)

    def __init__(self, parent=None):
        super(TemperatureEfectividad, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flow Arrangement")), 1, 1)
        self.flow=QtGui.QComboBox()
        for text, flow in self.flujo:
            self.flow.addItem(text)
        self.flow.currentIndexChanged.connect(self.plot)
        layout.addWidget(self.flow,1,2)
        self.mixed=QtGui.QComboBox()
        for text in self.mezclado:
            self.mixed.addItem(text)
        self.mixed.currentIndexChanged.connect(self.changeMixed)
        layout.addWidget(self.mixed,1,3)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,4)

        self.diagrama = mpl(self, width=10, height=10)
        self.image=self.diagrama.fig.figimage([[0]], 0, 0, zorder=1)
        logo=image.imread('images/pychemqt_98.png')
        self.logo=self.diagrama.fig.figimage(logo, 0, 0, zorder=1)

        self.botonGuardar=QtGui.QToolButton()
        self.botonGuardar.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")))
        self.botonGuardar.setToolTip(QtGui.QApplication.translate("pychemqt", "Save chart to file"))
        self.botonGuardar.clicked.connect(self.save)
        layout.addWidget(self.botonGuardar,5,1)
        layout.addWidget(self.diagrama,2,1,1,6)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 2, 1, 5)

        self.plot(0)
        self.showMaximized()
        self.refixImage()

    def changeMixed(self, indice):
        modelo=self.flow.currentIndex()
        self.plot(modelo)


    def plot(self, indice):
        self.diagrama.axes2D.clear()
        self.diagrama.axes2D.set_xlim(0.1, 10)
        self.diagrama.axes2D.set_ylim(0, 1)
        self.diagrama.axes2D.set_xscale("log")

        self.diagrama.axes2D.set_title(QtGui.QApplication.translate("pychemqt", "Heat Transfer Temperature Effectiveness"), size='12')
        self.diagrama.axes2D.set_xlabel("NTU", size='12')
        self.diagrama.axes2D.set_ylabel("P", size='14')
        self.diagrama.axes2D.set_xticklabels(["0.1", "1.0", "10"])
        self.diagrama.axes2D.set_xticklabels(["0.2", "0.3", "", "0.5", "", "0.7", "", "", "2.0", "3.0", "", "5.0", "", "7.0", "", ""], minor=True)

        flujo=self.flujo[indice][1]
        self.mixed.setVisible(flujo=="CrFSMix")
        kwargs={}
        if flujo=="CrFSMix":
            kwargs["mixed"]=str(self.mixed.currentText())

        R=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 6., 8., 10., 15.]

        NTU=logspace(-1.5, 1, 100)
        for  ri in R:
            e=[0]+[Heat_ExchangerDesign.TemperatureEffectiveness(N, ri, flujo, **kwargs) for N in NTU[1:]]
            self.diagrama.plot(NTU, e, "k")
            self.diagrama.axes2D.annotate(" R=%0.1f" %ri, (NTU[-1], e[-1]), size="medium", horizontalalignment="left", verticalalignment="center")

#        F=[0.3]
#        for f in F:
#            p=[]
#            NTU=[]
#            for r in R:
#                func=lambda P: Heat_ExchangerDesign.CorrectionFactor(P, r, flujo)-f
#                print func(0.)
#                pi=fsolve(func, 0.2)
#                p.append(pi)
#                NTU.append(Heat_ExchangerDesign.NTU_fPR(p, r, flujo))

#            p=[fsolve(lambda P: Heat_ExchangerDesign.CorrectionFactor(P, r, flujo)-f, 0.5)[0] for r in R]
#            NTU=[Heat_ExchangerDesign.NTU_fPR(pi, ri) for pi, ni in zip(p, R)]
#            self.diagrama.plot(NTU, p, "--")

        self.diagrama.draw()

        if flujo=="CrFSMix" and self.mixed.currentIndex():
            img=image.imread('images/equation/%s2.png' %flujo)
        else:
            img=image.imread('images/equation/%s.png' %flujo)
        self.image.set_data(img)
        self.refixImage()


class F(Chart):
    title=QtGui.QApplication.translate("pychemqt", u"ΔT Correction Factor", None, QtGui.QApplication.UnicodeUTF8)

    flujo=[(QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
                (QtGui.QApplication.translate("pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
                (QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
                (QtGui.QApplication.translate("pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    mezclado=("Cmin", "Cmax")
    PosLogo=(1, 1)
    PosImage=(0, 0)

    def __init__(self, parent=None):
        super(F, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flow Arrangement")), 1, 1)
        self.flow=QtGui.QComboBox()
        for text, flow in self.flujo:
            self.flow.addItem(text)
        self.flow.currentIndexChanged.connect(self.plot)
        layout.addWidget(self.flow,1,2)
        self.mixed=QtGui.QComboBox()
        for text in self.mezclado:
            self.mixed.addItem(text)
        self.mixed.currentIndexChanged.connect(self.changeMixed)
        layout.addWidget(self.mixed,1,3)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,4)

        self.diagrama = mpl(self, width=10, height=10)
        self.image=self.diagrama.fig.figimage([[0]], 0, 0, zorder=1)
        logo=image.imread('images/pychemqt_98.png')
        self.logo=self.diagrama.fig.figimage(logo, 0, 0, zorder=1)

        self.botonGuardar=QtGui.QToolButton()
        self.botonGuardar.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")))
        self.botonGuardar.setToolTip(QtGui.QApplication.translate("pychemqt", "Save chart to file"))
        self.botonGuardar.clicked.connect(self.save)
        layout.addWidget(self.botonGuardar,5,1)
        layout.addWidget(self.diagrama,2,1,1,6)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 2, 1, 5)

        self.plot(0)
        self.showMaximized()
        self.refixImage()

    def changeMixed(self, indice):
        modelo=self.flow.currentIndex()
        self.plot(modelo)

    def plot(self, indice):
        self.diagrama.axes2D.clear()
        self.diagrama.axes2D.set_xlim(0, 1)
        self.diagrama.axes2D.set_ylim(0, 1)
        self.diagrama.axes2D.set_title(QtGui.QApplication.translate("pychemqt", u"$\Delta T_{ml}$ Correction Factor", None, QtGui.QApplication.UnicodeUTF8), size='12')
        self.diagrama.axes2D.set_xlabel("$P=\\frac{T_{1o}-T_{1i}}{T_{2i}-T_{1i}}$", size='12')
        self.diagrama.axes2D.set_ylabel("F", size='14')

        flujo=self.flujo[indice][1]
#        self.mixed.setVisible(flujo=="CrFSMix")
        kwargs={}
        if flujo=="CrFSMix":
            kwargs["mixed"]=str(self.mixed.currentText())

        R=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 6, 8, 10, 15, 20]
#        R=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

        P=arange(0, 1.01, 0.01)
        for ri in R:
            f=[Heat_ExchangerDesign.CorrectionFactor(p, ri, flujo, **kwargs) for p in P]
            self.diagrama.plot(P, f, "k")

#            fraccionx=P[90]-P[80]
#            fracciony=f[90]-f[80]
#            try:
#                angle=arctan(fracciony/fraccionx)*360/2/pi
#            except ZeroDivisionError:
#                angle=90
#            self.diagrama.axes2D.annotate("R=%0.1f" %ri, (P[90], f[90]), rotation=angle, size="medium", horizontalalignment="left", verticalalignment="bottom")

        self.diagrama.draw()

        img=image.imread('images/equation/%s.png' %flujo)
        self.image.set_data(img)
        self.refixImage()

class Fi(Chart):
    title=QtGui.QApplication.translate("pychemqt", u"ψ", None, QtGui.QApplication.UnicodeUTF8)

    flujo=[#(QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
#                (QtGui.QApplication.translate("pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
#                (QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
                (QtGui.QApplication.translate("pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    mezclado=("Cmin", "Cmax")
    PosLogo=(1, 1)

    def __init__(self, parent=None):
        super(Fi, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flow Arrangement")), 1, 1)
        self.flow=QtGui.QComboBox()
        for text, flow in self.flujo:
            self.flow.addItem(text)
        self.flow.currentIndexChanged.connect(self.plot)
        layout.addWidget(self.flow,1,2)
        self.mixed=QtGui.QComboBox()
        for text in self.mezclado:
            self.mixed.addItem(text)
        self.mixed.currentIndexChanged.connect(self.changeMixed)
        layout.addWidget(self.mixed,1,3)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,4)

        self.diagrama = mpl(self, width=10, height=10)
        self.image=self.diagrama.fig.figimage([[0]], 0, 0, zorder=1)
        logo=image.imread('images/pychemqt_98.png')
        self.logo=self.diagrama.fig.figimage(logo, 0, 0, zorder=1)

        self.botonGuardar=QtGui.QToolButton()
        self.botonGuardar.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")))
        self.botonGuardar.setToolTip(QtGui.QApplication.translate("pychemqt", "Save chart to file"))
        self.botonGuardar.clicked.connect(self.save)
        layout.addWidget(self.botonGuardar,5,1)
        layout.addWidget(self.diagrama,2,1,1,6)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 2, 1, 5)

        self.plot(0)
        self.showMaximized()
        self.refixImage()

    def changeMixed(self, indice):
        modelo=self.flow.currentIndex()
        self.plot(modelo)

    def plot(self, indice):
        self.diagrama.axes2D.clear()
        self.diagrama.axes2D.set_xlim(0, 1)
        self.diagrama.axes2D.set_ylim(0, 1)
        self.diagrama.axes2D.set_title(QtGui.QApplication.translate("pychemqt", u"$\Delta T_{ml}$ Correction Factor", None, QtGui.QApplication.UnicodeUTF8), size='12')
        self.diagrama.axes2D.set_xlabel("$P=\\frac{T_{1o}-T_{1i}}{T_{2i}-T_{1i}}$", size='12')
        self.diagrama.axes2D.set_ylabel("F", size='14')

        flujo=self.flujo[indice][1]
#        self.mixed.setVisible(flujo=="CrFSMix")
        kwargs={}
        if flujo=="CrFSMix":
            kwargs["mixed"]=str(self.mixed.currentText())

        R=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 6, 8, 10]

        P=arange(0, 1.01, 0.01)
        for ri in R:
            f=[Heat_ExchangerDesign.Fi(p, ri, flujo, **kwargs) for p in P]
            self.diagrama.plot(P, f, "k")

        NTU=[0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.]
        for ntu in NTU:
            self.diagrama.plot([0, 1], [0, 1./ntu], "k", linestyle=":")
#            fraccionx=P[90]-P[80]
#            fracciony=f[90]-f[80]
#            try:
#                angle=arctan(fracciony/fraccionx)*360/2/pi
#            except ZeroDivisionError:
#                angle=90
#            self.diagrama.axes2D.annotate("R=%0.1f" %ri, (P[90], f[90]), rotation=angle, size="medium", horizontalalignment="left", verticalalignment="bottom")

        self.diagrama.draw()

        img=image.imread('images/equation/%s.png' %flujo)
        self.image.set_data(img)
        self.refixImage()


chart=(Efectividad, TemperatureEfectividad, F, Fi)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtGui.QApplication(sys.argv)
#    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
#    cambiador=Heat_Exchanger(entrada=agua, Tout=90+273.15)
    dialogo = UI_equipment()
#    dialogo = Fi()
    dialogo.show()
    sys.exit(app.exec_())
