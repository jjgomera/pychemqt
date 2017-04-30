#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Module for graphics elements of PFD
###############################################################################


from functools import partial
from datetime import datetime
import tempfile
import os
import json
import subprocess
from copy import deepcopy
from xml.dom import minidom

from PyQt5 import QtCore, QtGui, QtSvg, QtWidgets

from lib import unidades
from lib.project import Project
from lib.thread import WaitforClick
from lib.config import Preferences
from lib.corriente import Corriente
from UI import texteditor, UI_corriente
from UI.plots import Plot_Distribucion
from UI.widgets import createAction, Table_Graphics, PathConfig
from UI.prefPFD import ConfLineDialog
from equipment import *
from equipment.parents import equipment

# Value for mouse wheel zoom
factor = 5


class SelectStreamProject(QtWidgets.QDialog):
    project = None

    def __init__(self, parent=None):
        super(SelectStreamProject, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Select stream from file"))

        layout = QtWidgets.QGridLayout(self)
        label = QtWidgets.QApplication.translate("pychemqt", "Project path")
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "Select pychemqt project file")
        patrones = []
        patrones.append(QtWidgets.QApplication.translate(
            "pychemqt", "pychemqt project file") + " (*.pcq)")
        patron = ";;".join(patrones)
        self.filename = PathConfig(label + ":", msg=msg, patron=patron)
        self.filename.valueChanged.connect(self.changeproject)
        layout.addWidget(self.filename, 1, 1, 1, 3)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Streams")), 2, 1)
        self.stream = QtWidgets.QComboBox()
        layout.addWidget(self.stream, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 3, 3)

        self.status = QtWidgets.QLabel()
        layout.addWidget(self.status, 10, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(False)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 2, 1, 2)

    def changeproject(self, path):
        st = QtWidgets.QApplication.translate("pychemqt", "Loading project...")
        self.status.setText(st)
        QtWidgets.QApplication.processEvents()
        try:
            with open(path, "r") as file:
                self.project = json.load(file)
        except Exception as e:
            print(e)
            self.status.setText(QtGui.QApplication.translate(
                "pychemqt", "Failed to loading project..."))
        self.buttonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(True)
        self.status.setText(QtWidgets.QApplication.translate(
            "pychemqt", "Project loaded succesfully"))
        self.stream.clear()
        for stream in sorted(self.project["stream"].keys()):
            self.stream.addItem(stream)


class TextItemDlg(QtWidgets.QDialog):
    """Dialogo de edición de textos del diagrama de flujo"""

    def __init__(self, text=None, parent=None):
        super(TextItemDlg, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        self.editor = texteditor.TextEditor()
        self.editor.notas.textChanged.connect(self.updateUi)
        layout.addWidget(self.editor, 1, 1, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 2, 1, 1, 1)
        self.editor.notas.setFocus()
        if text:
            self.editor.setText(text)
        self.setWindowTitle(QtWidgets.QApplication.translate("pychemqt", "Edit text"))
        self.updateUi()

    def updateUi(self):
        self.buttonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(bool(self.editor.notas.toPlainText()))


class GeometricItem(object):
    """Clase genérica con la funcionalidad comun de los elementos geométricos"""

    def __init__(self, parent=None):
        super(GeometricItem, self).__init__(parent)
        self.setPen(self._pen())
        self.setFlags(
            QtWidgets.QGraphicsItem.ItemIsSelectable | QtWidgets.QGraphicsItem.ItemIsMovable | QtWidgets.QGraphicsItem.ItemSendsGeometryChanges | QtWidgets.QGraphicsItem.ItemIsFocusable)
        self.setZValue(-1)

    def _pen(self):
        pen = QtGui.QPen(QtGui.QColor(Preferences.get("PFD", 'Color_Stream')))
        pen.setWidthF(Preferences.getfloat("PFD", 'Width'))
        pen.setJoinStyle(
            [QtCore.Qt.MiterJoin, QtCore.Qt.BevelJoin, QtCore.Qt.RoundJoin][Preferences.getint("PFD", 'Union')])
        pen.setMiterLimit(Preferences.getfloat("PFD", 'Miter_limit'))
        pen.setCapStyle(
            [QtCore.Qt.FlatCap, QtCore.Qt.RoundCap, QtCore.Qt.SquareCap][Preferences.getint("PFD", 'Punta')])
        pen.setStyle([QtCore.Qt.SolidLine, QtCore.Qt.DashLine, QtCore.Qt.DotLine, QtCore.Qt.DashDotLine,
                      QtCore.Qt.DashDotDotLine][Preferences.getint("PFD", 'Guion')])
        pen.setDashOffset(Preferences.getfloat("PFD", 'Dash_offset'))
        return pen

    def delete(self):
        self.scene().delete(self)

    def format(self):
        dialog = ConfLineDialog(self.pen())
        if dialog.exec_():
            pen = dialog.pen()
            self.setPen(pen)
            self.itemChange(QtWidgets.QGraphicsItem.ItemPositionChange, 0)

    def contextMenu(self):
        contextMenu = QtWidgets.QMenu("%s Item" % self.type, self.scene().parent())
        contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(self.icon)))
        contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/editDelete.png"),
                              QtWidgets.QApplication.translate("pychemqt", "Delete"), self.delete)
        contextMenu.addSeparator()
        contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Appearance"), self.format)
        return contextMenu

    def itemChange(self, change, variant):
        if self.scene():
            if change == QtWidgets.QGraphicsItem.ItemPositionChange:
                if self.scene().parent().dirty[self.scene().parent().idTab] == False:
                    self.scene().parent().dirty[self.scene().parent().idTab] = True
                    self.scene().parent().activeControl(True)
                    self.scene().parent().tabModified(self.scene().parent().idTab)
        return QtWidgets.QGraphicsItem.itemChange(self, change, variant)

    def keyPressEvent(self, event):
        if event.modifiers() & QtCore.Qt.ShiftModifier:
            if event.key() == QtCore.Qt.Key_Up:
                rect = self.rect()
                rect.setBottom(self.rect().bottom() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key_Down:
                rect = self.rect()
                rect.setBottom(self.rect().bottom() + factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key_Left:
                rect = self.rect()
                rect.setRight(self.rect().right() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key_Right:
                rect = self.rect()
                rect.setRight(self.rect().right() + factor)
                self.setRect(rect)
        else:
            if event.key() == QtCore.Qt.Key_Delete or event.key() == QtCore.Qt.Key_Backspace:
                self.delete()
            elif event.key() == QtCore.Qt.Key_Escape:
                self.setSelected(False)
            elif event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
                self.setCurrentCell(self.currentRow() - 1, self.currentColumn())
            elif event.key() == QtCore.Qt.Key_Up:
                rect = self.rect()
                rect.moveTop(self.rect().y() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key_Down:
                rect = self.rect()
                rect.moveTop(self.rect().y() + factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key_Left:
                rect = self.rect()
                rect.moveLeft(self.rect().x() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key_Right:
                rect = self.rect()
                rect.moveLeft(self.rect().x() + factor)
                self.setRect(rect)


class RectItem(GeometricItem, QtWidgets.QGraphicsRectItem):
    """Clase que define un rectangulo"""
    type = "square"
    icon = os.environ["pychemqt"] + "/images/equipment/square.png"


class EllipseItem(GeometricItem, QtWidgets.QGraphicsEllipseItem):
    """Clase que define una circunferencia"""
    type = "ellipse"
    icon = os.environ["pychemqt"] + "/images/equipment/circle.png"


class TextItem(QtWidgets.QGraphicsTextItem):
    """Clase que define los textos en el diagrama de flujo"""
    type = "txt"

    def __init__(self, text, parent=None, position=QtCore.QPointF(0, 0), transform=QtGui.QTransform(), selectable=True):
        super(TextItem, self).__init__(parent=parent)
        if selectable:
            self.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable | QtWidgets.QGraphicsItem.ItemIsMovable |
                          QtWidgets.QGraphicsItem.ItemSendsGeometryChanges | QtWidgets.QGraphicsItem.ItemIsFocusable)
        else:
            self.setFlags(QtWidgets.QGraphicsItem.ItemIsMovable)
        self.setHtml(text)
        self.setPos(position)
        self.setTransform(transform)
        self.selectable = selectable

    def delete(self):
        self.scene().delete(self)

    def mouseDoubleClickEvent(self, event=None):
        dialog = TextItemDlg(self.toHtml())
        if dialog.exec_():
            self.setHtml(dialog.editor.texto)
            self.itemChange(QtWidgets.QGraphicsItem.ItemPositionChange, 0)

    def contextMenu(self):
        if self.selectable:
            contextMenu = QtWidgets.QMenu(
                QtWidgets.QApplication.translate("pychemqt", "Text Item: %s" % self.toPlainText()),
                self.scene().parent())
            contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"] + "/images/equipment/text.png")))
            contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/editDelete.png"),
                                  QtWidgets.QApplication.translate("pychemqt", "Delete"), self.delete)
            contextMenu.addSeparator()
            contextMenu.addAction("Edit", self.mouseDoubleClickEvent)
            return contextMenu

    def itemChange(self, change, variant):
        if self.scene():
            if change == QtWidgets.QGraphicsItem.ItemPositionChange:
                if self.scene().parent().dirty[self.scene().parent().idTab] == False:
                    self.scene().parent().dirty[self.scene().parent().idTab] = True
                    self.scene().parent().activeControl(True)
                    self.scene().parent().tabModified(self.scene().parent().idTab)
        return QtWidgets.QGraphicsItem.itemChange(self, change, variant)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Delete or event.key() == QtCore.Qt.Key_Backspace:
            self.delete()
        elif event.key() == QtCore.Qt.Key_Escape:
            self.setSelected(False)
        elif event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
            self.mouseDoubleClickEvent()
        elif event.key() == QtCore.Qt.Key_Up:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y() - factor))
        elif event.key() == QtCore.Qt.Key_Down:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y() + factor))
        elif event.key() == QtCore.Qt.Key_Left:
            self.setPos(QtCore.QPointF(self.pos().x() - factor, self.pos().y()))
        elif event.key() == QtCore.Qt.Key_Right:
            self.setPos(QtCore.QPointF(self.pos().x() + factor, self.pos().y()))


class GraphicsEntity(object):
    """Clase que modela la funcionalidad comun a corrientes y equipos en el PFD"""

    def view(self):
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".txt", encoding="utf-8") as temp:
            temp.write(QtWidgets.QApplication.translate("pychemqt", "Project Name") + ": " +
                       self.scene().parent().currentFilename + os.linesep)
            if isinstance(self.entity, Corriente):
                temp.write(QtWidgets.QApplication.translate("pychemqt", "Stream Id"))
            else:
                temp.write(QtWidgets.QApplication.translate("pychemqt", "Equipment Id"))
            temp.write(": %i" % self.id + os.linesep)
            ahora = datetime.today()
            temp.write(QtWidgets.QApplication.translate("pychemqt", "Report generated at") + " %s - %s" % (
                ahora.strftime("%H:%M:%S"), ahora.strftime("%d/%m/%Y")) + os.linesep)
            temp.write(self.entity.txt())
            subprocess.Popen([Preferences.get("Applications", 'TextViewer'), temp.name])

    def exportExcel(self):
        msg = QtWidgets.QApplication.translate("pychemqt", "Select Spreadsheet")
        patrones = []
        if os.environ["ezodf"]:
            patrones.append(QtWidgets.QApplication.translate("pychemqt", "Libreoffice spreadsheet files") + " (*.ods)")
        if os.environ["xlwt"]:
            patrones.append(
                QtWidgets.QApplication.translate("pychemqt", "Microsoft Excel 97/2000/XP/2003 XMLL") + " (*.xls)")
        if os.environ["openpyxl"]:
            patrones.append(QtWidgets.QApplication.translate("pychemqt", "Microsoft Excel 2007/2010 XML") + " (*.xlsx)")
        patron = ";;".join(patrones)
        dir = os.path.dirname(str(self.scene().parent().currentFilename))
        ruta = str(QtWidgets.QFileDialog.getSaveFileName(self.scene().parent(), msg, dir, patron)[0])
        if ruta:
            name, ext = os.path.splitext(ruta)
            if not ext or ext not in (".ods", ".xlsx", ".xls"):
                ruta += "." + str(patrones[0]).split(".")[-1][:-1]

            if ruta[-3:] == "ods":
                import ezodf
                templatefile = os.environ["pychemqt"] + os.sep + "dat" + os.sep + "templates" + os.sep + self.entity.__class__.__name__.lower() + ".ots"
                if os.path.isfile(templatefile):
                    spreadsheet = ezodf.newdoc("ods", ruta, templatefile)
                    sheet = spreadsheet.sheets[0]
                    for attr, type, cell in self.entity.datamap2xls():
                        prop = self.entity._prop(attr)
                        if type == "value":
                            value = prop.config()
                        else:
                            value = prop.text()
                        sheet[cell].set_value(value)
                else:
                    spreadsheet = ezodf.newdoc("ods", ruta)
                    sheets = spreadsheet.sheets
                    sheet = ezodf.Table('pychemqt - s%i' % self.id)
                    sheets += sheet
                    propiedades = self.entity.properties()
                    sheet.reset(size=(len(propiedades) + 1, 10))

                    for i, (name, attr, unit) in enumerate(self.entity.propertiesNames()):
                        value = self.entity._prop(attr)
                        txt = ""
                        if unit in unidades._all and not isinstance(value, list):
                            txt = value.text()
                            value = value.config()
                        elif isinstance(value, list):
                            txt = value[0].text()

                        sheet["B%i" % (i + 2)].set_value(name)
                        sheet["C%i" % (i + 2)].set_value(value)
                        sheet["D%i" % (i + 2)].set_value(txt)

                spreadsheet.save()
            elif ruta[-4:] == "xlsx":
                print(ruta, "is xlsx")
            elif ruta[-3:] == "xls":
                print(ruta, "is xls")


class StreamItem(GeometricItem, QtWidgets.QGraphicsPathItem, GraphicsEntity):
    """Clase que define una corriente gráficamente"""
    up = None
    down = None
    id = 0
    free_id = []
    type = "stream"

    def __init__(self, parent=None):
        super(StreamItem, self).__init__()
        self.parent = parent
        self.setPen(self._pen())
        qp = QtGui.QPainterPath()
        self.setPath(qp)
        self.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable | QtWidgets.QGraphicsItem.ItemIsFocusable)
        if StreamItem.free_id:
            self.id = StreamItem.free_id.pop(0)
        else:
            self.id = StreamItem.id + 1
            StreamItem.id += 1
        self.idLabel = TextItem("S%i" % self.id, self, selectable=False)
        self.idLabel.setZValue(2)
        self.setAcceptHoverEvents(True)

    @property
    def corriente(self):
        return self.scene().project.getStream(self.id)

    @property
    def entity(self):
        return self.corriente

    def setCorriente(self, corriente):
        self.scene().project.setStream(self.id, corriente)
        kwargs = {"entrada": corriente}
        if isinstance(self.scene().project.getDownToStream(self.id), flux.Mixer):
            kwargs["id_entrada"] = self.scene().project.streams[self.id][3] + 1
        equip = self.scene().project.getDownToStream(self.id)
        if isinstance(equip, equipment):
            equip(**kwargs)
        pen = self.pen()
        if corriente.status == 1:
            pen.setColor(QtGui.QColor("blue"))
        else:
            pen.setColor(QtGui.QColor("red"))
        self.setPen(pen)
        self.itemChange(QtWidgets.QGraphicsItem.ItemPositionChange, 0)

    def mouseDoubleClickEvent(self, event=None):
        dialog = UI_corriente.Corriente_Dialog(self.corriente)
        if dialog.exec_():
            self.setCorriente(dialog.corriente)

    def copyFromProject(self):
        dialog = SelectStreamProject()
        if dialog.exec_():
            indice = dialog.stream.currentText()
            data = dialog.project["stream"][indice]
            corriente = Corriente()
            corriente.readFromJSON(data)
            self.setCorriente(corriente)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Delete or event.key() == QtCore.Qt.Key_Backspace:
            self.delete()
        elif event.key() == QtCore.Qt.Key_Escape:
            self.setSelected(False)
        elif event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
            self.mouseDoubleClickEvent()

    def hoverEnterEvent(self, event):
        if not (self.scene().addObj and self.scene().addType == "stream"):
            self.tabla = Table_Graphics(self.corriente, self.id, self.scene().parent().Preferences)
            self.tabla.move(event.screenPos())
            self.tabla.show()

    def hoverLeaveEvent(self, event):
        if not (self.scene().addObj and self.scene().addType == "stream"):
            self.tabla.hide()
            self.tabla.deleteLater()

            #    def hoverMoveEvent(self, event):
            #        self.tabla.move(event.screenPos())
            #        self.tabla.show()

    def contextMenu(self):
        ViewAction = createAction(QtWidgets.QApplication.translate("pychemqt", "View Properties"), slot=self.view,
                                  parent=self.scene())
        SolidDistributionAction = createAction(QtWidgets.QApplication.translate("pychemqt", "Solid Distribution Fit"),
                                               slot=self.solidFit, parent=self.scene())
        if self.corriente:
            if not self.corriente.solido:
                SolidDistributionAction.setEnabled(False)
        else:
            ViewAction.setEnabled(False)
            SolidDistributionAction.setEnabled(False)

        contextMenu = QtWidgets.QMenu("Stream %i" % self.id, self.scene().parent())
        contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"] + "/images/equipment/stream.png")))
        contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Copy from another project"),
                              self.copyFromProject)
        contextMenu.addAction(SolidDistributionAction)
        contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Edit"), self.mouseDoubleClickEvent)
        contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/editDelete.png"),
                              QtWidgets.QApplication.translate("pychemqt", "Delete"), self.delete)
        contextMenu.addSeparator()
        contextMenu.addAction(ViewAction)
        contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Export to spreadsheet"), self.exportExcel)
        contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Show/Hide Id Label"),
                              self.idLabelVisibility)
        contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Appearance"), self.format)
        return contextMenu

    def redraw(self, entrada=None, salida=None):
        if entrada:
            self.entrada = entrada
        if salida:
            self.salida = salida

        max_height = max(self.up.boundingRect().height(), self.down.boundingRect().height())
        max_width = max(self.up.boundingRect().width(), self.down.boundingRect().width())
        y_up = min(self.up.pos().y(), self.down.pos().y())
        y_down = max(self.up.pos().y(), self.down.pos().y())
        if self.up.pos().y() == y_up:
            height_sup = self.up.boundingRect().height()
        else:
            height_sup = self.down.boundingRect().height()
        x_mean = (self.entrada.x() + self.salida.x()) / 2.
        Xdist_entrada = 0
        Ydist_entrada = 0
        Xdist_salida = 0
        Ydist_salida = 0
        if self.Ang_entrada == 0:
            Xdist_entrada = 20
        elif self.Ang_entrada == 180:
            Xdist_entrada = -20
        elif self.Ang_entrada == 90:
            Ydist_entrada = 20
        elif self.Ang_entrada == 360:
            Ydist_entrada = -20

        if self.Ang_salida == 0:
            Xdist_salida = -20
        elif self.Ang_salida == 180:
            Xdist_salida = 20

        qp = QtGui.QPainterPath()
        qp.moveTo(self.entrada)
        if self.Ang_entrada == self.Ang_salida:
            if self.salida.x() > self.entrada.x() + 10:
                qp.lineTo(QtCore.QPointF(x_mean, self.entrada.y()))
                qp.lineTo(QtCore.QPointF(x_mean, self.salida.y()))
            else:
                if abs(self.entrada.y() - self.salida.y()) > max_height:  # cabe la linea entre ambos equipos
                    y_mean = (y_up + y_down + height_sup) / 2.
                    qp.lineTo(QtCore.QPointF(self.entrada.x() + Xdist_entrada, self.entrada.y() + Ydist_entrada))
                    qp.lineTo(QtCore.QPointF(self.entrada.x() + Xdist_entrada, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x() + Xdist_salida, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x() + Xdist_salida, self.salida.y() + Ydist_entrada))
                else:  # sacamos la linea por encima de los equipos
                    y_mean = y_up - 20
                    qp.lineTo(QtCore.QPointF(self.entrada.x() + Xdist_entrada, self.entrada.y() + Ydist_entrada))
                    qp.lineTo(QtCore.QPointF(self.entrada.x() + Xdist_entrada, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x() + Xdist_salida, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x() + Xdist_salida, self.salida.y() + Ydist_salida))
        else:
            x_mean = max(self.entrada.x(), self.salida.x()) + Xdist_salida
            qp.lineTo(QtCore.QPointF(x_mean, self.entrada.y()))
            qp.lineTo(QtCore.QPointF(x_mean, self.salida.y()))

        qp.lineTo(self.salida)
        self.prepareGeometryChange()
        self.setPath(qp)
        if abs(self.entrada.y() - self.salida.y()) <= 30:
            self.idLabel.setPos(x_mean, max(self.entrada.y(), self.salida.y()))
        else:
            self.idLabel.setPos(x_mean, (self.entrada.y() + self.salida.y()) / 2. - 10)

    def postDelete(self):
        StreamItem.free_id.append(self.id)
        self.up.down_used -= 1
        self.down.up_used -= 1
        self.up.down.remove(self)
        self.down.up.remove(self)

    def idLabelVisibility(self):
        self.idLabel.setVisible(not self.idLabel.isVisible())

    def solidFit(self):
        if self.corriente.solido:
            dialog = Plot_Distribucion(self.id, self.corriente.solido)
            self.scene().parent().currentMdi.addSubWindow(dialog)
            dialog.show()


class EquipmentItem(QtSvg.QGraphicsSvgItem, GraphicsEntity):
    """Clase que define los equipos en el diagrama de flujo"""
    up = []
    down = []
    up_used = 0
    down_used = 0
    id = 0
    id_in = 0
    id_out = 0
    type = "equip"

    def __init__(self, name, dialogoId, parent=None):
        self.name = name
        imagen = os.environ["pychemqt"] + "images/equipment/%s.svg" % name
        super(EquipmentItem, self).__init__(imagen, parent=parent)
        self.dialogoId = dialogoId
        self.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable | QtWidgets.QGraphicsItem.ItemIsMovable |
                      QtWidgets.QGraphicsItem.ItemSendsGeometryChanges | QtWidgets.QGraphicsItem.ItemIsFocusable)
        self.imagen = imagen
        self.angle = 0
        self.setAcceptHoverEvents(True)

        if dialogoId != None:
            self.dialogo = UI_equipments[dialogoId].UI_equipment
            EquipmentItem.id += 1
            self.id = EquipmentItem.id
            self.tipo = "e"
            self.idLabel = TextItem("E%i" % self.id, self, selectable=False)
            self.idLabel.setPos(self.boundingRect().width() / 3., -20)

        else:
            self.dialogo = UI_corriente.Corriente_Dialog
            if name == "in":
                EquipmentItem.id_in += 1
                self.id = EquipmentItem.id_in
                self.tipo = "i"
            else:
                EquipmentItem.id_out += 1
                self.id = EquipmentItem.id_out
                self.tipo = "o"

        output = []
        input = []
        doc = minidom.parse(imagen)

        for entrada in doc.getElementsByTagName("ins")[0].childNodes:
            if isinstance(entrada, minidom.Element):
                if entrada.tagName == "in":
                    x = float(entrada.getAttribute("x"))
                    y = float(entrada.getAttribute("y"))
                    d = float(entrada.getAttribute("d"))
                    input.append([x, y, d])

        for salida in doc.getElementsByTagName("outs")[0].childNodes:
            if isinstance(salida, minidom.Element):
                if salida.tagName == "out":
                    x = float(salida.getAttribute("x"))
                    y = float(salida.getAttribute("y"))
                    d = float(salida.getAttribute("d"))
                    output.append([x, y, d])
        doc.unlink()

        self.input = []
        if input:
            for entrada in input:
                obj = QtWidgets.QGraphicsEllipseItem(self)
                obj.setRect(entrada[0] * self.boundingRect().width() - 5,
                            entrada[1] * self.boundingRect().height() - 5,
                            10, 10)
                obj.direction = int(entrada[2])
                obj.setPen(QtGui.QColor(255, 255, 255))
                obj.setBrush(QtGui.QColor(Preferences.get("PFD", 'Color_Entrada')))
                self.input.append(obj)
        self.output = []
        if output:
            for salida in output:
                obj = QtWidgets.QGraphicsEllipseItem(self)
                obj.setRect(salida[0] * self.boundingRect().width() - 5,
                            salida[1] * self.boundingRect().height() - 5,
                            10, 10)
                obj.direction = int(salida[2])
                obj.setPen(QtGui.QColor(255, 255, 255))
                obj.setBrush(QtGui.QColor(Preferences.get("PFD", 'Color_Salida')))
                self.output.append(obj)
        self.showInput(False)

    @property
    def equipment(self):
        return self.scene().project.getItem(self.id)

    @property
    def entity(self):
        return self.equipment

    def mouseDoubleClickEvent(self, event=None):
        if self.dialogoId != None:
            kwarg = {"equipment": self.equipment}
            # Aditional parameters for selected equipment
            if isinstance(self.equipment, flux.Divider):
                # Divider
                if not len(self.down):
                    return
                kwarg["salidas"] = len(self.down)
            elif isinstance(self.equipment, flux.Mixer):
                # mixer
                if not len(self.up):
                    return
                kwarg["entradas"] = len(self.up)
            elif isinstance(self.equipment, spreadsheet.Spreadsheet):
                # spreadsheet
                self.equipment(project=self.scene().project)
                kwarg["project"] = self.scene().project

            dialog = self.dialogo(**kwarg)
            if dialog.exec_():
                self.scene().project.setItem(self.id, dialog.Equipment)
                #                self.up[0].setCorriente(dialog.Equipment.entrada)
                for i, corriente in enumerate(dialog.Equipment.salida):
                    self.down[i].setCorriente(corriente)

                self.itemChange(QtWidgets.QGraphicsItem.ItemPositionChange, 0)
        else:
            if self.output:
                self.down[0].mouseDoubleClickEvent()
            else:
                self.up[0].mouseDoubleClickEvent()

    def mousePressEvent(self, event):
        QtSvg.QGraphicsSvgItem.mousePressEvent(self, event)
        if self.scene().addObj:
            if self.scene().addType == "stream":
                if len(self.scene().Pos) == 0:
                    punto = self.output
                    self.scene().up = self
                    x = self.down_used
                else:
                    punto = self.input
                    self.scene().down = self
                    x = self.up_used
                self.scene().Pos.append(self.mapToScene(punto[x].rect().center()))
                self.scene().points.append(punto[x])
            else:
                self.scene().Pos.append(event.pos())

    def mouseMoveEvent(self, event=None):
        if event:
            QtWidgets.QGraphicsPixmapItem.mouseMoveEvent(self, event)
        for i, corriente in enumerate(self.up):
            corriente.redraw(salida=self.mapToScene(self.input[i].rect().center()))
        for i, corriente in enumerate(self.down):
            corriente.redraw(entrada=self.mapToScene(self.output[i].rect().center()))

    def showInput(self, bool):
        for entrada in self.input:
            entrada.setVisible(bool)
        for salida in self.output:
            salida.setVisible(bool)

    def hoverEnterEvent(self, event):
        if self.scene().addObj and self.scene().addType == "stream":
            self.showInput(True)
        else:
            if self.dialogoId != None:
                self.tabla = Table_Graphics(self.equipment, self.id, self.scene().parent().Preferences)
            else:
                if self.output:
                    self.tabla = Table_Graphics(self.down[0].corriente, self.down[0].id,
                                                self.scene().parent().Preferences)
                else:
                    self.tabla = Table_Graphics(self.up[0].corriente, self.up[0].id,
                                                self.scene().parent().Preferences)
            self.tabla.move(event.screenPos())
            self.tabla.show()

    def hoverLeaveEvent(self, event):
        self.showInput(False)
        if not self.scene().addObj and self.scene().addType == "stream":
            self.tabla.hide()
            self.tabla.deleteLater()

            #    def hoverMoveEvent(self, event):
            #        self.tabla.move(event.screenPos())
            #        self.tabla.show()
            #

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Delete or event.key() == QtCore.Qt.Key_Backspace:
            self.delete()
        elif event.key() == QtCore.Qt.Key_Escape:
            self.setSelected(False)
        elif event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
            self.mouseDoubleClickEvent()
        elif event.key() == QtCore.Qt.Key_Up:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y() - factor))
            self.mouseMoveEvent()
        elif event.key() == QtCore.Qt.Key_Down:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y() + factor))
            self.mouseMoveEvent()
        elif event.key() == QtCore.Qt.Key_Left:
            self.setPos(QtCore.QPointF(self.pos().x() - factor, self.pos().y()))
            self.mouseMoveEvent()
        elif event.key() == QtCore.Qt.Key_Right:
            self.setPos(QtCore.QPointF(self.pos().x() + factor, self.pos().y()))
            self.mouseMoveEvent()

    def itemChange(self, change, variant):
        if self.scene():
            if change == QtWidgets.QGraphicsItem.ItemPositionChange:
                if self.scene().parent().dirty[self.scene().parent().idTab] == False:
                    self.scene().parent().dirty[self.scene().parent().idTab] = True
                    self.scene().parent().activeControl(True)
                    self.scene().parent().tabModified(self.scene().parent().idTab)
        return QtWidgets.QGraphicsItem.itemChange(self, change, variant)

    def contextMenu(self):
        if self.dialogoId != None:
            ViewAction = createAction(QtWidgets.QApplication.translate("pychemqt", "View Properties"),
                                      slot=self.view,
                                      parent=self.scene())
            ViewAction.setEnabled(self.equipment.status)

            contextMenu = QtWidgets.QMenu("Equipment %i" % self.id, self.scene().parent())
            contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"] +
                                                          "/images/equipment/%s.svg" % self.name)))
            contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Edit"), self.mouseDoubleClickEvent)
            contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/editDelete.png"),
                                  QtWidgets.QApplication.translate("pychemqt", "Delete"), self.delete)
            contextMenu.addSeparator()
            contextMenu.addAction(ViewAction)
            contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Export to spreadsheet"),
                                  self.exportExcel)
            contextMenu.addAction(QtWidgets.QApplication.translate("pychemqt", "Show/Hide Id Label"),
                                  self.idLabelVisibility)
        # contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Appearance"), self.format)
        #            contextMenu.addSeparator()
        #            contextMenu.addAction("Run", self.mouseDoubleClickEvent)
        else:
            if self.output:
                contextMenu = self.down[0].contextMenu()
            else:
                contextMenu = self.up[0].contextMenu()

        self.menuTransform = QtWidgets.QMenu(QtWidgets.QApplication.translate("pychemqt", "Transform"))
        self.menuTransform.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/transform_rotate_90.png"),
                                     QtWidgets.QApplication.translate("pychemqt", "Rotate by 90"),
                                     partial(self.rotate, 90))
        self.menuTransform.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/transform_rotate_180.png"),
                                     QtWidgets.QApplication.translate("pychemqt", "Rotate by 180"),
                                     partial(self.rotate, 180))
        self.menuTransform.addAction(QtGui.QIcon(os.environ["pychemqt"] + "/images/button/transform_rotate_270.png"),
                                     QtWidgets.QApplication.translate("pychemqt", "Rotate by 270"),
                                     partial(self.rotate, 270))
        self.menuTransform.addSeparator()
        self.menuTransform.addAction(QtWidgets.QApplication.translate("pychemqt", "Mirror about X"),
                                     partial(self.rotate, 270))
        self.menuTransform.addAction(QtWidgets.QApplication.translate("pychemqt", "Mirror about Y"),
                                     partial(self.rotate, 270))
        contextMenu.addAction(self.menuTransform.menuAction())

        return contextMenu

    def delete(self):
        self.scene().delete(self)

    def format(self):
        pass

    def postDelete(self):
        while self.down:
            stream = self.down.pop()
            self.scene().delete(stream)
        while self.up:
            stream = self.up.pop()
            self.scene().delete(stream)

    def idLabelVisibility(self):
        self.idLabel.setVisible(not self.idLabel.isVisible())

    def rotate(self, angle):
        self.angle = angle
        transform = self.transform()
        transform.rotate(angle)
        self.setTransform(transform)
        self.mouseMoveEvent()
        for i, entrada in enumerate(self.up):
            new_angle = (self.input[i].direction + angle) % 360
            self.input[i].direction = new_angle
            entrada.Ang_salida = new_angle
            entrada.redraw()
        for i, salida in enumerate(self.down):
            new_angle = (self.output[i].direction + angle) % 360
            self.output[i].direction = new_angle
            salida.Ang_entrada = new_angle
            salida.redraw()


class GraphicsView(QtWidgets.QGraphicsView):
    mouseMove = QtCore.pyqtSignal(QtCore.QPointF)

    #    mouseClick = QtCore.pyqtSignal("QtCore.QPointF")
    def __init__(self, PFD=True, parent=None):
        super(GraphicsView, self).__init__(parent)
        self.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setRenderHint(QtGui.QPainter.TextAntialiasing)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("#aaaaaa"), QtCore.Qt.Dense7Pattern))
        self.setMouseTracking(True)
        # self.setAcceptDrops(True)
        self.PFD = PFD

    #    def wheelEvent(self, event):
    #        print event.delta()
    #
    #        factor = 1.41 ** (-event.delta() / 240.0)
    #        self.zoom(factor)

    def mouseMoveEvent(self, event):
        QtWidgets.QGraphicsView.mouseMoveEvent(self, event)
        self.mouseMove.emit(event.globalPos())

    def mousePressEvent(self, event):
        QtWidgets.QGraphicsView.mousePressEvent(self, event)
        if not self.PFD:
            self.scene().views()[0].centerOn(self.mapToScene(event.pos()))

    def closeEvent(self, event):
        if self.PFD:
            event.ignore()

    # def dragEnterEvent(self, event):
        # if event.mimeData().hasFormat("application/x-equipment"):
            # event.accept()
        # else:
            # event.ignore()

    # def dragMoveEvent(self, event):
        # if event.mimeData().hasFormat("application/x-equipment"):
            # event.setDropAction(QtCore.Qt.CopyAction)
            # event.accept()
        # else:
            # event.ignore()

    # def dropEvent(self, event):
        # if event.mimeData().hasFormat("application/x-equipment"):
            # data = event.mimeData().data("application/x-equipment")
            # stream = QtCore.QDataStream(data, QtCore.QIODevice.ReadOnly)
            # icon = QtGui.QIcon()
            # stream >> icon
            # event.setDropAction(QtCore.Qt.CopyAction)
            # print(event.pos())
            # event.accept()
            # self.updateGeometry()
            # self.update()
        # else:
            # event.ignore()

    def zoom(self, value):
        factor = value / 100.0
        self.resetMatrix()
        self.scale(factor, factor)


class GraphicsScene(QtWidgets.QGraphicsScene):
    copiedItem = QtCore.QByteArray()
    pasteOffset = 5
    points = []
    addObj = False
    addType = ""
    project = Project()
    objects = {"txt": [], "square": [], "ellipse": [], "stream": {}, "in": {}, "out": {}, "equip": {}}

    def __init__(self, parent=None):
        super(GraphicsScene, self).__init__(parent)

    def mousePressEvent(self, event):
        QtWidgets.QGraphicsScene.mousePressEvent(self, event)
        if self.addObj and self.addType != "stream":
            self.Pos.append(event.scenePos())

    def addActions(self, menu, pos=None):
        menu.addAction(QtWidgets.QApplication.translate("pychemqt", "Redraw"), self.update)
        menu.addSeparator()
        menu.addAction(createAction(QtWidgets.QApplication.translate("pychemqt", "Select All"),
                                    slot=self.selectAll,
                                    shortcut=QtGui.QKeySequence.SelectAll,
                                    icon=os.environ["pychemqt"] + "/images/button/selectAll",
                                    parent=self))
        menu.addSeparator()
        actionCut = createAction(QtWidgets.QApplication.translate("pychemqt", "Cut"),
                                 slot=self.cut,
                                 shortcut=QtGui.QKeySequence.Cut,
                                 icon=os.environ["pychemqt"] + "/images/button/editCut",
                                 parent=self)
        menu.addAction(actionCut)
        actionCopy = createAction(QtWidgets.QApplication.translate("pychemqt", "Copy"),
                                  slot=self.copy,
                                  shortcut=QtGui.QKeySequence.Copy,
                                  icon=os.environ["pychemqt"] + "/images/button/editCopy",
                                  parent=self)
        menu.addAction(actionCopy)
        actionPaste = createAction(QtWidgets.QApplication.translate("pychemqt", "Paste"),
                                   slot=partial(self.paste, pos),
                                   shortcut=QtGui.QKeySequence.Paste,
                                   icon=os.environ["pychemqt"] + "/images/button/editPaste",
                                   parent=self)
        menu.addAction(actionPaste)
        actionDelete = createAction(QtWidgets.QApplication.translate("pychemqt", "Delete All"),
                                    slot=self.delete,
                                    shortcut=QtGui.QKeySequence.Delete,
                                    icon=os.environ["pychemqt"] + "/images/button/editDelete",
                                    parent=self)
        menu.addAction(actionDelete)
        menu.addSeparator()

        if self.copiedItem.isEmpty():
            actionPaste.setEnabled(False)
        items = self.selectedItems()
        if not items:
            actionCut.setEnabled(False)
            actionCopy.setEnabled(False)
            actionDelete.setEnabled(False)

        for item in items:
            menuEl = item.contextMenu()
            menu.addAction(menuEl.menuAction())
        return menu

    def contextMenuEvent(self, event):
        item = self.itemAt(event.scenePos(), self.views()[0].transform())
        if item:
            item.setSelected(True)
        contextMenu = QtWidgets.QMenu()
        self.addActions(contextMenu, event.scenePos())
        contextMenu.exec_(event.screenPos())

    def selectAll(self):
        for item in list(self.items()):
            item.setSelected(True)

    def copy(self, item=None):
        if not item:
            item = self.selectedItems()[0]
        self.copiedItem.clear()
        self.pasteOffset = 5
        stream = QtCore.QDataStream(self.copiedItem, QtCore.QIODevice.WriteOnly)
        self.writeItemToStream(stream, item)

    def cut(self):
        item = self.selectedItems()[0]
        self.copy(item)
        self.removeItem(item)
        del item

    def paste(self, pos=None):
        stream = QtCore.QDataStream(self.copiedItem, QtCore.QIODevice.ReadOnly)
        item = self.readItemFromStream(stream)
        if pos:
            item.setPos(pos)
        else:
            item.setPos(item.pos() + QtCore.QPointF(self.pasteOffset, self.pasteOffset))
            self.pasteOffset += 5
        self.addItem(item)

    def delete(self, items=None):
        if items:
            items = [items]
        else:
            items = self.selectedItems()
        for item in items:
            tipo = item.type
            if tipo in ["stream", "equip"]:
                item.postDelete()
                del self.objects[tipo][item.id]
            else:
                self.objects[tipo].remove(item)
            self.removeItem(item)
        self.update()
        self.parent().list.updateList(self.objects)

    def waitClick(self, numClick, type, object):
        self.object = object
        self.addType = type
        self.addObj = True
        self.views()[0].viewport().setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                                             "Click in desire text position in screen"))
        self.Pos = []
        self.clickCollector = WaitforClick(numClick, self)
        self.clickCollector.finished.connect(self.click)
        self.clickCollector.start()

    def click(self):
        if self.addType in ["equip", "in", "out", "txt"]:
            self.object.setPos(self.Pos[0])
        elif self.addType == "stream":
            self.object.up = self.up
            self.object.down = self.down
            self.up.down = self.up.down + [self.object]
            self.down.up = self.down.up + [self.object]
            self.up.down_used += 1
            self.down.up_used += 1
            self.object.Ang_entrada = self.points[0].direction
            self.object.Ang_salida = self.points[1].direction
            self.object.redraw(self.Pos[0], self.Pos[1])
        elif self.addType in ["square", "ellipse"]:
            rect = QtCore.QRectF(self.Pos[0], self.Pos[1])
            self.object.setRect(rect)

        self.addItem(self.object)
        if self.addType == "equip":
            self.project.addItem("e%i" % self.object.id, self.object.dialogo.Equipment)
        elif self.addType == "in":
            self.project.addItem("i%i" % self.object.id, self.object.dialogo.corriente)
        elif self.addType == "out":
            self.project.addItem("o%i" % self.object.id, self.object.dialogo.corriente)
        elif self.addType == "stream":
            self.project.addStream(self.object.id, "%s%i" % (self.up.tipo, self.up.id),
                                   "%s%i" % (self.down.tipo, self.down.id), Corriente(),
                                   self.up.down_used - 1, self.down.up_used - 1)

        self.parent().dirty[self.parent().idTab] = True
        self.parent().saveControl()
        self.update()
        self.object.setSelected(True)

        self.parent().statusbar.clearMessage()
        self.addObj = False
        if self.addType in ("txt", "square", "ellipse"):
            self.objects[self.addType].append(self.object)
        else:
            id = self.object.id
            self.objects[self.addType][id] = self.object
        self.parent().list.updateList(self.objects)

        self.views()[0].viewport().setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))

    #        for item in self.items():
    #            if isinstance(item, QtSvg.QGraphicsSvgItem):
    #                item.hoverLeaveEvent(None)
    #            item.setAcceptHoverEvents(False)




    def readItemFromStream(self, stream):
        type = QtCore.QString()
        matrix = QtGui.QTransform()
        stream >> type >> matrix
        if type == "txt":
            text = QtCore.QString()
            stream >> text
            item = TextItem(text)
        elif type == "square":
            rect = QtCore.QRectF()
            pen = QtGui.QPen()
            stream >> rect >> pen
            item = RectItem()
            item.setRect(rect)
            item.setPen(pen)
        elif type == "ellipse":
            rect = QtCore.QRectF()
            pen = QtGui.QPen()
            stream >> rect >> pen
            item = EllipseItem()
            item.setRect(rect)
            item.setPen(pen)
        elif type == "equip":
            name = QtCore.QString()
            stream >> name
            dialogoid = stream.readInt32()
            item = EquipmentItem(name, dialogoid)

        item.setTransform(matrix)
        return item

    def writeItemToStream(self, stream, item):
        stream << QtCore.QString(item.type) << item.transform()
        if isinstance(item, TextItem):
            stream << item.toHtml()
        elif isinstance(item, EllipseItem):
            stream << item.rect() << item.pen()
        elif isinstance(item, RectItem):
            stream << item.rect() << item.pen()
        elif isinstance(item, EquipmentItem):
            stream.writeString(item.name)
            # stream << QtCore.QString(item.name)
            stream.writeInt32(item.dialogoId)

    def readFromJSON(self, data):
        self.objects = deepcopy(GraphicsScene.objects)

        for text in data["PFD"]["txt"].values():
            txt = text["txt"]
            s = TextItem(txt)
            x = text["x"]
            y = text["y"]
            pos = QtCore.QPoint(x, y)
            s.setPos(pos)
            self.objects["txt"].append(s)
            self.addItem(s)

        for obj in data["PFD"]["square"].values():
            s = RectItem()
            x = obj["x"]
            y = obj["y"]
            width = obj["width"]
            height = obj["height"]
            rect = QtCore.QRect(x, y, width, height)
            s.setRect(rect)

            pen = QtGui.QPen(QtGui.QColor(obj["color"]))
            pen.setWidthF(obj["width"])
            pen.setJoinStyle(obj["joinStyle"])
            pen.setMiterLimit(obj["miterLimit"])
            pen.setCapStyle(obj["capStyle"])
            pen.setStyle(obj["style"])
            pen.setDashOffset(obj["dashOffset"])
            s.setPen(pen)

            self.objects["square"].append(s)
            self.addItem(s)

        for obj in data["PFD"]["ellipse"].values():
            s = EllipseItem()
            x = obj["x"]
            y = obj["y"]
            width = obj["width"]
            height = obj["height"]
            rect = QtCore.QRect(x, y, width, height)
            s.setRect(rect)

            pen = QtGui.QPen(QtGui.QColor(obj["color"]))
            pen.setWidthF(obj["width"])
            pen.setJoinStyle(obj["joinStyle"])
            pen.setMiterLimit(obj["miterLimit"])
            pen.setCapStyle(obj["capStyle"])
            pen.setStyle(obj["style"])
            pen.setDashOffset(obj["dashOffset"])
            s.setPen(pen)

            self.objects["ellipse"].append(s)
            self.addItem(s)

        id_stream = []
        up_stream = {}
        down_stream = {}
        for id, obj in data["PFD"]["stream"].items():
            id = int(id)
            s = StreamItem()
            in_x = obj["input_x"]
            in_y = obj["input_y"]
            entrada = QtCore.QPointF(in_x, in_y)
            out_x = obj["output_x"]
            out_y = obj["output_y"]
            salida = QtCore.QPointF(out_x, out_y)

            pen = QtGui.QPen(QtGui.QColor(obj["pen"]["color"]))
            pen.setWidthF(obj["pen"]["width"])
            pen.setJoinStyle(obj["pen"]["joinStyle"])
            pen.setMiterLimit(obj["pen"]["miterLimit"])
            pen.setCapStyle(obj["pen"]["capStyle"])
            pen.setStyle(obj["pen"]["style"])
            pen.setDashOffset(obj["pen"]["dashOffset"])
            s.setPen(pen)

            up_type = obj["up_type"]
            down_type = obj["down_type"]
            up_id = obj["up_id"]
            down_id = obj["down_id"]
            id_stream.append(id)
            up_stream[id] = up_type, up_id
            down_stream[id] = down_type, down_id

            s.id = id
            s.entrada = entrada
            s.salida = salida
            s.Ang_entrada = obj["input_angle"]
            s.Ang_salida = obj["output_angle"]
            self.objects["stream"][id_stream[-1]] = s
            self.addItem(s)

            txt = obj["label"]
            x = obj["label_x"]
            y = obj["label_y"]
            pos = QtCore.QPointF(x, y)
            s.idLabel.setPos(pos)
            s.idLabel.setHtml(txt)
            visible = obj["label_visible"]
            s.idLabel.setVisible(visible)

        angle_in = {}
        for id, obj in data["PFD"]["in"].items():
            id = int(id)
            s = EquipmentItem("in", None)
            x = obj["x"]
            y = obj["y"]
            pos = QtCore.QPointF(x, y)
            s.setPos(pos)

            angle_in[id] = obj["angle"]
            down = [self.objects["stream"][obj["down_id"]]]
            s.down = down
            self.objects["in"][id] = s
            self.addItem(s)

        angle_out = {}
        for id, obj in data["PFD"]["out"].items():
            id = int(id)
            s = EquipmentItem("out", None)
            x = obj["x"]
            y = obj["y"]
            pos = QtCore.QPointF(x, y)
            s.setPos(pos)

            angle_out[id] = obj["angle"]
            up = [self.objects["stream"][obj["up_id"]]]
            s.up = up
            self.objects["out"][id] = s
            self.addItem(s)

        angle_equip = {}
        for id, obj in data["PFD"]["equip"].items():
            id = int(id)

            name = obj["name"]
            dialogoId = obj["dialogo_id"]
            s = EquipmentItem(name, dialogoId)
            s.id = id

            x = obj["x"]
            y = obj["y"]
            pos = QtCore.QPointF(x, y)
            s.setPos(pos)

            angle_equip[id] = obj["angle"]
            up = [self.objects["stream"][i] for i in obj["up_id"]]
            s.up = up
            down = [self.objects["stream"][i] for i in obj["down_id"]]
            s.down = down

            self.objects["equip"][id] = s
            self.addItem(s)

            txt = obj["label"]
            x = obj["label_x"]
            y = obj["label_y"]
            pos = QtCore.QPointF(x, y)
            s.idLabel.setPos(pos)
            s.idLabel.setHtml(txt)
            visible = obj["label_visible"]
            s.idLabel.setVisible(visible)

        for id in id_stream:
            tipo, i = up_stream[id]
            self.objects["stream"][id].up = self.getObject(tipo, i)
            tipo, i = down_stream[id]
            self.objects["stream"][id].down = self.getObject(tipo, i)
            self.objects["stream"][id].redraw()

        for id, angle in angle_in.items():
            self.objects["in"][id].rotate(angle)
        for id, angle in angle_out.items():
            self.objects["out"][id].rotate(angle)
        for id, angle in angle_equip.items():
            self.objects["equip"][id].rotate(angle)

    def writeToJSON(self, data):
        txts = {}
        for i, obj in enumerate(self.objects["txt"]):
            txt = {}
            txt["txt"] = obj.toHtml()
            txt["x"] = obj.pos().x()
            txt["y"] = obj.pos().y()
            txts[i] = txt
        data["txt"] = txts

        squares = {}
        for i, obj in enumerate(self.objects["square"]):
            square = {}
            square["x"] = obj.rect().x()
            square["y"] = obj.rect().y()
            square["width"] = obj.rect().width()
            square["height"] = obj.rect().height()

            pen = {}
            pen["color"] = obj.pen().color().name()
            pen["width"] = obj.pen().widthF()
            pen["joinStyle"] = obj.pen().joinStyle()
            pen["miterLimit"] = obj.pen().miterLimit()
            pen["capStyle"] = obj.pen().capStyle()
            pen["style"] = obj.pen().style()
            pen["dashOffset"] = obj.pen().dashOffset()
            square["pen"] = pen
            squares[i] = square
        data["square"] = squares

        ellipses = {}
        for i, obj in enumerate(self.objects["ellipse"]):
            ellipse = {}
            ellipse["x"] = obj.rect().x()
            ellipse["y"] = obj.rect().y()
            ellipse["width"] = obj.rect().width()
            ellipse["height"] = obj.rect().height()
            pen = {}
            pen["color"] = obj.pen().color().name()
            pen["width"] = obj.pen().widthF()
            pen["joinStyle"] = obj.pen().joinStyle()
            pen["miterLimit"] = obj.pen().miterLimit()
            pen["capStyle"] = obj.pen().capStyle()
            pen["style"] = obj.pen().style()
            pen["dashOffset"] = obj.pen().dashOffset()
            ellipse["pen"] = pen
            ellipses[i] = ellipse
        data["ellipse"] = ellipses

        streams = {}
        for id, obj in self.objects["stream"].items():
            stream = {}
            stream["input_x"] = obj.entrada.x()
            stream["input_y"] = obj.entrada.y()
            stream["output_x"] = obj.salida.x()
            stream["output_y"] = obj.salida.y()
            pen = {}
            pen["color"] = obj.pen().color().name()
            pen["width"] = obj.pen().widthF()
            pen["joinStyle"] = obj.pen().joinStyle()
            pen["miterLimit"] = obj.pen().miterLimit()
            pen["capStyle"] = obj.pen().capStyle()
            pen["style"] = obj.pen().style()
            pen["dashOffset"] = obj.pen().dashOffset()
            stream["pen"] = pen

            stream["up_id"] = obj.up.id
            stream["up_type"] = obj.up.tipo
            stream["down_id"] = obj.down.id
            stream["down_type"] = obj.down.tipo
            stream["input_angle"] = obj.Ang_entrada
            stream["output_angle"] = obj.Ang_salida
            stream["label"] = obj.idLabel.toHtml()
            stream["label_x"] = obj.idLabel.pos().x()
            stream["label_y"] = obj.idLabel.pos().y()
            stream["label_visible"] = int(obj.idLabel.isVisible())
            streams[id] = stream
        data["stream"] = streams

        ins = {}
        for id, obj in self.objects["in"].items():
            in_ = {}
            in_["x"] = obj.x()
            in_["y"] = obj.y()
            in_["angle"] = obj.angle

            if obj.down:
                in_["down_id"] = obj.down[0].id
            else:
                in_["down_id"] = None
            ins[id] = in_
        data["in"] = ins

        outs = {}
        for id, obj in self.objects["out"].items():
            out = {}
            out["x"] = obj.x()
            out["y"] = obj.y()
            out["angle"] = obj.angle

            if obj.up:
                out["up_id"] = obj.up[0].id
            else:
                out["up_id"] = None
            outs[id] = out
        data["out"] = outs

        equipments = {}
        for id, obj in self.objects["equip"].items():
            equip = {}
            equip["name"] = obj.name
            equip["dialogo_id"] = obj.dialogoId
            equip["x"] = obj.pos().x()
            equip["y"] = obj.pos().y()
            equip["angle"] = obj.angle

            ups = []
            for up in obj.up:
                ups.append(up.id)
            equip["up_id"] = ups
            downs = []
            for down in obj.down:
                downs.append(down.id)
            equip["down_id"] = downs

            equip["label"] = obj.idLabel.toHtml()
            equip["label_x"] = obj.idLabel.pos().x()
            equip["label_y"] = obj.idLabel.pos().y()
            equip["label_visible"] = int(obj.idLabel.isVisible())
            equipments[id] = equip
        data["equip"] = equipments

    def getObject(self, tipo, id):
        if tipo == "e":
            lista = self.objects["equip"]
        elif tipo == "i":
            lista = self.objects["in"]
        elif tipo == "o":
            lista = self.objects["out"]
        elif tipo == "s":
            lista = self.objects["stream"]
        else:
            raise Exception
        return lista[id]


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    #    dialogo = ConfLineDialog()
    dialogo = SelectStreamProject()
    dialogo.show()
    sys.exit(app.exec_())
