#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=too-many-lines

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


Module for graphics elements in PFD, defining the main qt graphics class:

    * :class:`GraphicsView`: PFD viewer
    * :class:`GraphicsScene`: PFD scene to manage the graphical elements

and the different graphical items:

    * :class:`GeometricItem`: Common functionality for geometric PFD elements
    * :class:`RectItem`: Class to plot a rectangular item
    * :class:`EllipseItem`: Class to plot a circle or oval item
    * :class:`TextItem`: Class to plot a text item
    * :class:`GraphicsEntity`: Common functionality for Entity in PFD
    * :class:`StreamItem`: Class to plot a mass stream
    * :class:`EquipmentItem`: Class to plot equipment item

Dialog related defined:

    * :class:`SelectStreamProject`: Dialog to select a stream from a project
    * :class:`TextItemDlg`: Dialog to edit texts in PFD
'''


from ast import literal_eval
from copy import deepcopy
from datetime import datetime
from functools import partial
import json
import os
from random import randint
import subprocess
import tempfile
from xml.dom import minidom

try:
    import ezodf
except ImportError:
    pass

from lib import unidades
from lib.config import Preferences, IMAGE_PATH
from lib.corriente import Corriente
from lib.EoS import K, H
from lib.project import Project
from lib.thread import WaitforClick
from UI import texteditor, UI_corriente
from UI.plots import Plot_Distribucion
from UI.prefPFD import ConfLine, ConfLineDialog, Dialog, BrushCombo
from UI.widgets import createAction, Table_Graphics, PathConfig, ClickableLabel
from equipment import flux, spreadsheet, UI_equipments
from equipment.parents import equipment
from tools import UI_confResolution, UI_confThermo
from tools.qt import QtCore, QtGui, QtSvgWidgets, QtWidgets, translate


# Value for keyboard navigation, unnecessary add to configuration
factor = 5


class GraphicsView(QtWidgets.QGraphicsView):
    """Class for PFD representation"""
    mouseMove = QtCore.pyqtSignal(QtCore.QPointF)
    zoomChanged = QtCore.pyqtSignal(int)

    def __init__(self, PFD=True, parent=None):
        """
        Parameters
        ----------
        PFD : boolean
            Set if this GraphicsView is the PFD of project, in that case the
            windows can't be closed
        """
        super().__init__(parent)
        self.PFD = PFD
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.RubberBandDrag)
        self.setRenderHints(QtGui.QPainter.RenderHint.Antialiasing
                            | QtGui.QPainter.RenderHint.TextAntialiasing)

        brushColor = Preferences.get("PFD", "brushColor")
        brushPattern = BrushCombo.BRUSH[Preferences.getint("PFD", "brush")]
        self.setBackgroundBrush(
            QtGui.QBrush(QtGui.QColor(brushColor), brushPattern))

        self.setMouseTracking(True)
        # self.setAcceptDrops(True)

        # Widgets to show in the statusbar of mainwindow
        self.statusWidget = []
        self.statusPosition = ClickableLabel()
        self.statusPosition.setFrameShape(QtWidgets.QFrame.Shape.WinPanel)
        self.statusPosition.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        self.statusWidget.append(self.statusPosition)
        self.statusThermo = ClickableLabel()
        self.statusThermo.setFrameShape(QtWidgets.QFrame.Shape.WinPanel)
        self.statusThermo.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        self.statusWidget.append(self.statusThermo)

        self.statusPosition.clicked.connect(
            partial(parent.dialogConfig, UI_confResolution))
        self.statusThermo.clicked.connect(
            partial(parent.dialogConfig, UI_confThermo))
        self.mouseMove.connect(self.updatePosition)
        self.leaveEvent()

    def mouseMoveEvent(self, event):
        """Reimplement to print the mouse position on statusbar"""
        QtWidgets.QGraphicsView.mouseMoveEvent(self, event)
        self.mouseMove.emit(self.mapToScene(event.pos()))

    def mousePressEvent(self, event):
        """Reimplement to use in a overview window the mouse click position as
        the center position in the real PFD window"""
        QtWidgets.QGraphicsView.mousePressEvent(self, event)
        if not self.PFD:
            self.scene().views()[0].centerOn(self.mapToScene(event.pos()))

    def closeEvent(self, event):
        """Reimplement to avoid close window if it's the PFD window"""
        if self.PFD:
            event.ignore()

    def wheelEvent(self, event):
        """Change zoom of window, only work with PFD window, not in overview"""
        if self.PFD:
            ratio = 1.41 ** (-event.angleDelta().y() / 240.0)
            self.scale(ratio, ratio)
            self.zoomChanged.emit(self.transform().m11()*100)

    def leaveEvent(self, event=None):
        """Reimplement to set value of status position to the total size"""
        x = Preferences.getint("PFD", "x")
        y = Preferences.getint("PFD", "y")
        self.statusPosition.setText("%i, %i" % (x, y))

    def zoom(self, value):
        """Apply zoom value to a pfd view"""
        value /= 100.0
        self.resetTransform()
        self.scale(value, value)

    # def dragEnterEvent(self, event):
        # if event.mimeData().hasFormat("application/x-equipment"):
            # event.accept()
        # else:
            # event.ignore()

    # def dragMoveEvent(self, event):
        # if event.mimeData().hasFormat("application/x-equipment"):
            # event.setDropAction(QtCore.Qt.DropAction.CopyAction)
            # event.accept()
        # else:
            # event.ignore()

    # def dropEvent(self, event):
        # if event.mimeData().hasFormat("application/x-equipment"):
            # data = event.mimeData().data("application/x-equipment")
            # stream = QtCore.QDataStream(data, QtCore.QIODevice.ReadOnly)
            # icon = QtGui.QIcon()
            # stream >> icon
            # event.setDropAction(QtCore.Qt.DropAction.CopyAction)
            # print(event.pos())
            # event.accept()
            # self.updateGeometry()
            # self.update()
        # else:
            # event.ignore()

    def updatePosition(self, event):
        """Update text with mouse position"""
        self.statusPosition.setText("%i, %i" % (event.x(), event.y()))

    def changeStatusThermo(self, config):
        """Show thermodynamic option in statusbar"""
        if config.has_section("Thermo") and \
                config.has_section("Components"):
            components = literal_eval(config.get("Components", "components"))
            if config.getboolean("Thermo", "iapws") and \
                    config.getboolean("Thermo", "freesteam") and \
                    len(components) == 1 and components[0] == 62:
                txt = "Freesteam"
            elif config.getboolean("Thermo", "iapws") and \
                    len(components) == 1 and components[0] == 62:
                txt = "iapws97"
            elif config.getboolean("Thermo", "meos") and \
                    config.getboolean("Thermo", "refprop"):
                txt = "Refprop"
            elif config.getboolean("Thermo", "meos") and \
                    config.getboolean("Thermo", "coolprop"):
                txt = "CoolProp"
            elif config.getboolean("Thermo", "meos"):
                txt = "MEoS"
            else:
                txt = "K: %s  H: %s" % (
                    K[config.getint("Thermo", "K")].__status__,
                    H[config.getint("Thermo", "H")].__status__)

            self.statusThermo.setText(txt)


class GraphicsScene(QtWidgets.QGraphicsScene):
    """Graphics PFD scene too manage the graphical elements"""
    copiedItem = QtCore.QByteArray()
    pasteOffset = 5
    points = []
    addObj = False
    addType = ""
    project = Project()
    objects = {"txt": [], "square": [], "ellipse": [], "stream": {}, "in": {},
               "out": {}, "equip": {}}

    def __init__(self, parent=None):
        super().__init__(parent)

        self.popup = Table_Graphics()
        proxy = QtWidgets.QGraphicsProxyWidget()
        proxy.setWidget(self.popup)
        self.addItem(proxy)
        # self.popup.move(0, -self.popup.height())
        self.popup.hide()


    def mousePressEvent(self, event):
        """Save event pos to locate item"""
        QtWidgets.QGraphicsScene.mousePressEvent(self, event)
        if self.addObj and self.addType != "stream":
            self.Pos.append(event.scenePos())

    def addActions(self, menu, pos=None):
        """Define actions and state of its to add to context menuw"""
        actionCut, actionCopy, actionPaste = self.defineShortcut(pos)
        menu.addAction(self.tr("Redraw"), self.update)
        menu.addAction(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "configure.png")),
            self.tr("Configure"),
            self.configure)
        menu.addSeparator()
        menu.addAction(createAction(
            self.tr("Select All"),
            slot=self.selectAll,
            shortcut=QtGui.QKeySequence.StandardKey.SelectAll,
            icon=os.path.join("button", "selectAll.png"), parent=self))
        menu.addSeparator()
        actionCut = createAction(
            self.tr("Cut"),
            slot=self.cut,
            shortcut=QtGui.QKeySequence.StandardKey.Cut,
            icon=os.path.join("button", "editCut.png"), parent=self)
        menu.addAction(actionCut)
        actionCopy = createAction(
            self.tr("Copy"),
            slot=self.copy,
            shortcut=QtGui.QKeySequence.StandardKey.Copy,
            icon=os.path.join("button", "editCopy.png"), parent=self)
        menu.addAction(actionCopy)
        actionPaste = createAction(
            self.tr("Paste"),
            slot=partial(self.paste, pos),
            shortcut=QtGui.QKeySequence.StandardKey.Paste,
            icon=os.path.join("button", "editPaste.png"), parent=self)
        menu.addAction(actionPaste)
        actionDelete = createAction(
            self.tr("Delete All"),
            slot=self.delete,
            shortcut=QtGui.QKeySequence.StandardKey.Delete,
            icon=os.path.join("button", "editDelete.png"), parent=self)
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
        menu.addSeparator()

        return menu

    def defineShortcut(self, pos=None):
        """Add actions to QGraphicsView to enable keyboard shortcut"""
        # Delete actions of last run
        for action in self.views()[0].actions():
            self.views()[0].removeAction(action)

        # Add actions to QGraphicsView to enable key shortcut
        actionCut = createAction(
            self.tr("Cut"),
            slot=self.cut,
            shortcut=QtGui.QKeySequence.StandardKey.Cut,
            icon=os.path.join("button", "editCut.png"), parent=self)
        actionCopy = createAction(
            self.tr("Copy"),
            slot=self.copy,
            shortcut=QtGui.QKeySequence.StandardKey.Copy,
            icon=os.path.join("button", "editCopy.png"), parent=self)
        actionPaste = createAction(
            self.tr("Paste"),
            slot=partial(self.paste, pos),
            shortcut=QtGui.QKeySequence.StandardKey.Paste,
            icon=os.path.join("button", "editPaste.png"), parent=self)
        self.views()[0].addAction(actionCopy)
        self.views()[0].addAction(actionCut)
        self.views()[0].addAction(actionPaste)
        return actionCut, actionCopy, actionPaste

    def contextMenuEvent(self, event):
        """Create the context menu to show then right click"""
        item = self.itemAt(event.scenePos(), self.views()[0].transform())
        if item:
            item.setSelected(True)
        contextMenu = QtWidgets.QMenu()
        self.addActions(contextMenu, event.scenePos())
        contextMenu.exec(event.screenPos())

    def selectAll(self):
        """Select all itemin scene"""
        for item in list(self.items()):
            item.setSelected(True)

    def copy(self, items=None):
        """Copy selected items to internal clickboard, StreamItem is not
        suppoerted"""
        if not items:
            items = self.selectedItems()
        self.copiedItem.clear()
        self.pasteOffset = 50
        st = QtCore.QDataStream(
            self.copiedItem, QtCore.QIODevice.OpenModeFlag.WriteOnly)

        st.writeInt32(len(items))
        for item in items:
            self.writeItemToStream(st, item)

    def cut(self):
        """Copy selected items to internal clickboard and delete of scene"""
        item = self.selectedItems()[0]
        self.copy(item)
        self.removeItem(item)
        del item

    def paste(self, pos=None):
        """Paste item saved in internal clipboard to the scene"""
        st = QtCore.QDataStream(
            self.copiedItem, QtCore.QIODevice.OpenModeFlag.ReadOnly)
        count = st.readInt32()
        for it in range(count):
            item = self.readItemFromStream(st)

            # Discard stream item
            if isinstance(item, StreamItem):
                continue

            if pos:
                item.setPos(pos)
            else:
                offset = QtCore.QPointF(
                    randint(-item.pos().x(), item.pos().x()),
                    randint(-item.pos().y(), item.pos().y()))
                item.setPos(item.pos() + offset)
                self.pasteOffset += 5
            self.addItem(item)

            # Add new equipment to project
            if item.tipo == "e":
                if item.name == "in":
                    self.project.copyInput(item.id-1)
                elif item.name == "out":
                    self.project.setOutput(item.id, Corriente())
                else:
                    self.project.copyItem(item.id-1)

    def delete(self, items=None):
        if items:
            items = [items]
        else:
            items = self.selectedItems()
        for item in items:
            tipo = item.tipo
            if tipo in ["stream", "equip"]:
                item.postDelete()
                del self.objects[tipo][item.id]
            else:
                self.objects[tipo].remove(item)
            self.removeItem(item)
        self.update()
        self.parent().list.updateList(self.objects)

    def configure(self):
        dlg = Dialog(Preferences)
        if dlg.exec():
            preferences = dlg.value(Preferences)
            self.parent().updatePreferences(preferences)

    def waitClick(self, numClick, tipo, object):
        self.Pos = []
        self.object = object
        self.addType = tipo
        self.addObj = True
        self.views()[0].viewport().setCursor(
            QtGui.QCursor(QtCore.Qt.CursorShape.CrossCursor))
        self.parent().statusBar().showMessage(
            self.tr("Click in desire text position in screen"))
        clickCollector = WaitforClick(numClick, self)
        clickCollector.finished.connect(self.click)
        clickCollector.start()

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

            # Unchecked stream button in toolbox
            self.parent().botonCorriente.setChecked(False)

        elif self.addType in ["square", "ellipse"]:
            rect = QtCore.QRectF(self.Pos[0], self.Pos[1])
            self.object.setRect(rect)

        self.addItem(self.object)
        if self.addType == "equip":
            self.project.addItem(
                "e%i" % self.object.id, self.object.dialogo.Equipment)
        elif self.addType == "in":
            self.project.addItem(
                "i%i" % self.object.id, self.object.dialogo.corriente)
        elif self.addType == "out":
            self.project.addItem(
                "o%i" % self.object.id, self.object.dialogo.corriente)
        elif self.addType == "stream":
            self.project.addStream(
                self.object.id, "%s%i" % (self.up.tipo, self.up.id),
                "%s%i" % (self.down.tipo, self.down.id), Corriente(),
                self.up.down_used - 1, self.down.up_used - 1)

        self.parent().dirty[self.parent().idTab] = True
        self.parent().saveControl()
        self.update()
        self.object.setSelected(True)

        self.parent().statusBar().clearMessage()
        self.addObj = False
        if self.addType in ("txt", "square", "ellipse"):
            self.objects[self.addType].append(self.object)
        else:
            id = self.object.id
            self.objects[self.addType][id] = self.object
        self.parent().list.updateList(self.objects)

        self.views()[0].viewport().setCursor(
            QtGui.QCursor(QtCore.Qt.CursorShape.ArrowCursor))
    #        for item in self.items():
    #            if isinstance(item, QtSvg.QGraphicsSvgItem):
    #                item.hoverLeaveEvent(None)
    #            item.setAcceptHoverEvents(False)

    def writeItemToStream(self, stream, item):
        stream.writeQString(item.tipo)
        stream << item.transform() << item.pos()
        if isinstance(item, TextItem):
            stream.writeQString(item.toHtml())
        elif isinstance(item, (RectItem, EllipseItem)):
            stream << item.rect() << item.pen()
        elif isinstance(item, EquipmentItem):
            stream.writeQString(item.name)
            if item.tipo == "e":
                stream.writeInt32(item.dialogoId)

    def readItemFromStream(self, stream):
        tipo = stream.readQString()
        matrix = QtGui.QTransform()
        stream >> matrix
        pos = QtCore.QPointF()
        stream >> pos
        if tipo == "txt":
            txt = stream.readQString()
            item = TextItem(txt)
        elif tipo == "square":
            rect = QtCore.QRectF()
            pen = QtGui.QPen()
            stream >> rect >> pen
            item = RectItem()
            item.setRect(rect)
            item.setPen(pen)
        elif tipo == "ellipse":
            rect = QtCore.QRectF()
            pen = QtGui.QPen()
            stream >> rect >> pen
            item = EllipseItem()
            item.setRect(rect)
            item.setPen(pen)
        elif tipo in ("e", "i", "o"):
            name = stream.readQString()
            if tipo == "e":
                dialogoid = stream.readInt32()
            else:
                dialogoid = None
            item = EquipmentItem(name, dialogoid)
        elif tipo == "stream":
            item = StreamItem()

        item.setTransform(matrix)
        item.setPos(pos)
        return item

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
            pen.setJoinStyle(ConfLine.join[obj["joinStyle"]])
            pen.setMiterLimit(obj["miterLimit"])
            pen.setCapStyle(ConfLine.cap[obj["capStyle"]])
            pen.setStyle(ConfLine.line[obj["style"]])
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
            pen.setJoinStyle(ConfLine.join[obj["joinStyle"]])
            pen.setMiterLimit(obj["miterLimit"])
            pen.setCapStyle(ConfLine.cap[obj["capStyle"]])
            pen.setStyle(ConfLine.line[obj["style"]])
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

            pen.setJoinStyle(ConfLine.join[obj["pen"]["joinStyle"]])
            pen.setMiterLimit(obj["pen"]["miterLimit"])
            pen.setCapStyle(ConfLine.cap[obj["pen"]["capStyle"]])
            pen.setStyle(ConfLine.line[obj["pen"]["style"]])
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
        """Save json format to write to file"""
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
            square["pen"] = obj.getPen()
            squares[i] = square
        data["square"] = squares

        ellipses = {}
        for i, obj in enumerate(self.objects["ellipse"]):
            ellipse = {}
            ellipse["x"] = obj.rect().x()
            ellipse["y"] = obj.rect().y()
            ellipse["width"] = obj.rect().width()
            ellipse["height"] = obj.rect().height()
            ellipse["pen"] = obj.getPen()
            ellipses[i] = ellipse
        data["ellipse"] = ellipses

        streams = {}
        for idx, obj in self.objects["stream"].items():
            stream = {}
            stream["input_x"] = obj.entrada.x()
            stream["input_y"] = obj.entrada.y()
            stream["output_x"] = obj.salida.x()
            stream["output_y"] = obj.salida.y()
            stream["pen"] = obj.getPen()

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
            streams[idx] = stream
        data["stream"] = streams

        ins = {}
        for idx, obj in self.objects["in"].items():
            in_ = {}
            in_["x"] = obj.x()
            in_["y"] = obj.y()
            in_["angle"] = obj.angle

            if obj.down:
                in_["down_id"] = obj.down[0].id
            else:
                in_["down_id"] = None
            ins[idx] = in_
        data["in"] = ins

        outs = {}
        for idx, obj in self.objects["out"].items():
            out = {}
            out["x"] = obj.x()
            out["y"] = obj.y()
            out["angle"] = obj.angle

            if obj.up:
                out["up_id"] = obj.up[0].id
            else:
                out["up_id"] = None
            outs[idx] = out
        data["out"] = outs

        equipments = {}
        for idx, obj in self.objects["equip"].items():
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
            equipments[idx] = equip
        data["equip"] = equipments

    def getObject(self, tipo, id):
        """Return the object of project"""
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


class GeometricItem():
    """Generic class with common functionality for geometric PFD elements"""
    tipo = ""
    icon = None

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setPen(self._pen())
        self.setFlags(
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsMovable
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsFocusable)
        self.setZValue(-1)

    def _pen(self):
        """Load pen properties from preferences"""
        pen = QtGui.QPen(QtGui.QColor(Preferences.get("PFD", 'Color_Stream')))
        pen.setWidthF(Preferences.getfloat("PFD", 'Width'))
        pen.setJoinStyle(ConfLine.join[Preferences.getint("PFD", 'Union')])
        pen.setMiterLimit(Preferences.getfloat("PFD", 'Miter_limit'))
        pen.setCapStyle(ConfLine.cap[Preferences.getint("PFD", 'Punta')])
        pen.setStyle(ConfLine.line[Preferences.getint("PFD", 'Guion')])
        pen.setDashOffset(Preferences.getfloat("PFD", 'Dash_offset'))
        return pen

    def getPen(self):
        """Get item pen properties in a json serializable format"""
        pen = {}
        pen["color"] = self.pen().color().name()
        pen["width"] = self.pen().widthF()
        pen["joinStyle"] = ConfLine.join.index(self.pen().joinStyle())
        pen["miterLimit"] = self.pen().miterLimit()
        pen["capStyle"] = ConfLine.cap.index(self.pen().capStyle())
        pen["style"] = ConfLine.line.index(self.pen().style())
        pen["dashOffset"] = self.pen().dashOffset()
        return pen

    def delete(self):
        """Delete item from scene"""
        self.scene().delete(self)

    def format(self):
        """Configure formating of line"""
        dialog = ConfLineDialog(self.pen())
        if dialog.exec():
            pen = dialog.pen()
            self.setPen(pen)

    def contextMenu(self):
        contextMenu = QtWidgets.QMenu(
            "%s Item" % self.tipo, self.scene().parent())
        contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(self.icon)))
        contextMenu.addAction(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "editDelete.png")),
            self.tr("Delete"),
            self.delete)
        contextMenu.addSeparator()
        contextMenu.addAction(
            self.tr("Appearance"),
            self.format)
        return contextMenu

    def itemChange(self, key, value):
        if self.scene() and key == \
                QtWidgets.QGraphicsItem.GraphicsItemChange.ItemPositionChange:
            mainwindow = self.scene().parent()
            if mainwindow.dirty[mainwindow.idTab] == False:
                mainwindow.dirty[mainwindow.idTab] = True
                mainwindow.activeControl(True)
                mainwindow.tabModified(mainwindow.idTab)
        return QtWidgets.QGraphicsItem.itemChange(self, key, value)

    def keyPressEvent(self, event):
        if event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier:
            if event.key() == QtCore.Qt.Key.Key_Up:
                rect = self.rect()
                rect.setBottom(self.rect().bottom() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key.Key_Down:
                rect = self.rect()
                rect.setBottom(self.rect().bottom() + factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key.Key_Left:
                rect = self.rect()
                rect.setRight(self.rect().right() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key.Key_Right:
                rect = self.rect()
                rect.setRight(self.rect().right() + factor)
                self.setRect(rect)
        else:
            if event.key() == QtCore.Qt.Key.Key_Delete or \
                    event.key() == QtCore.Qt.Key.Key_Backspace:
                self.delete()
            elif event.key() == QtCore.Qt.Key.Key_Escape:
                self.setSelected(False)
            elif event.key() == QtCore.Qt.Key.Key_Return or \
                    event.key() == QtCore.Qt.Key.Key_Enter:
                self.setCurrentCell(self.currentRow()-1, self.currentColumn())
            elif event.key() == QtCore.Qt.Key.Key_Up:
                rect = self.rect()
                rect.moveTop(self.rect().y() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key.Key_Down:
                rect = self.rect()
                rect.moveTop(self.rect().y() + factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key.Key_Left:
                rect = self.rect()
                rect.moveLeft(self.rect().x() - factor)
                self.setRect(rect)
            elif event.key() == QtCore.Qt.Key.Key_Right:
                rect = self.rect()
                rect.moveLeft(self.rect().x() + factor)
                self.setRect(rect)


class RectItem(GeometricItem, QtWidgets.QGraphicsRectItem):
    """Class to plot a rectangular item"""
    tipo = "square"
    icon = os.path.join(IMAGE_PATH, "equipment", "square.png")


class EllipseItem(GeometricItem, QtWidgets.QGraphicsEllipseItem):
    """Class to plot a circle or oval item"""
    tipo = "ellipse"
    icon = os.path.join(IMAGE_PATH, "equipment", "cirle.png")


class TextItem(QtWidgets.QGraphicsTextItem):
    """Class to plot a text item"""
    tipo = "txt"

    def __init__(self, text, parent=None, position=QtCore.QPointF(0, 0),
                 transform=QtGui.QTransform(), selectable=True):
        super().__init__(parent=parent)
        if selectable:
            self.setFlags(
                QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
                | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsMovable
                | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges
                | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsFocusable)
        else:
            self.setFlags(
                QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setHtml(text)
        self.setPos(position)
        self.setTransform(transform)
        self.selectable = selectable

    def delete(self):
        """Delete item from scene"""
        self.scene().delete(self)

    def mouseDoubleClickEvent(self, event=None):
        """Show dialog to change text"""
        dialog = TextItemDlg(self.toHtml())
        if dialog.exec():
            self.setHtml(dialog.editor.texto)

    def contextMenu(self):
        """Define context menu of item"""
        if self.selectable:
            contextMenu = QtWidgets.QMenu(
                self.tr("Text Item: %s" % self.toPlainText()),
                self.scene().parent())
            contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(
                os.path.join(IMAGE_PATH, "equipment", "text.png"))))
            contextMenu.addAction(
                QtGui.QIcon(os.path.join(
                    IMAGE_PATH, "button", "editDelete.png")),
                self.tr("Delete"),
                self.delete)
            contextMenu.addSeparator()
            contextMenu.addAction("Edit", self.mouseDoubleClickEvent)
            return contextMenu

    def itemChange(self, key, value):
        if self.scene() and key == \
                QtWidgets.QGraphicsItem.GraphicsItemChange.ItemPositionChange:
            mainwindow = self.scene().parent()
            if mainwindow.dirty[mainwindow.idTab] is False:
                mainwindow.dirty[mainwindow.idTab] = True
                mainwindow.activeControl(True)
                mainwindow.tabModified(mainwindow.idTab)
        return QtWidgets.QGraphicsItem.itemChange(self, key, value)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key.Key_Delete or \
                event.key() == QtCore.Qt.Key.Key_Backspace:
            self.delete()
        elif event.key() == QtCore.Qt.Key.Key_Escape:
            self.setSelected(False)
        elif event.key() == QtCore.Qt.Key.Key_Return or \
                event.key() == QtCore.Qt.Key.Key_Enter:
            self.mouseDoubleClickEvent()
        elif event.key() == QtCore.Qt.Key.Key_Up:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()-factor))
        elif event.key() == QtCore.Qt.Key.Key_Down:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()+factor))
        elif event.key() == QtCore.Qt.Key.Key_Left:
            self.setPos(QtCore.QPointF(self.pos().x()-factor, self.pos().y()))
        elif event.key() == QtCore.Qt.Key.Key_Right:
            self.setPos(QtCore.QPointF(self.pos().x()+factor, self.pos().y()))


class GraphicsEntity():
    """Class with common functionality for Entity in PFD"""

    def view(self):
        """Generate text report with properties calculated of entity"""
        with tempfile.NamedTemporaryFile(
                "w", delete=False, suffix=".txt", encoding="utf-8") as temp:
            title = self.tr("Project Name") + ": "
            title += self.scene().parent().currentFilename
            temp.write(title + os.linesep)

            if isinstance(self.entity, Corriente):
                temp.write(self.tr("Stream Id"))
            else:
                temp.write(self.tr("Equipment Id"))
            temp.write(": %i" % self.id + os.linesep)

            now = datetime.today()
            temp.write(self.tr("Report generated at"))
            temp.write(now.strftime("%H:%M:%S - %d/%m/%Y") + os.linesep)
            temp.write(self.entity.txt())
            subprocess.Popen(
                [Preferences.get("Applications", 'TextViewer'), temp.name])

    def exportExcel(self):
        """Export data to spreadsheet file"""
        msg = self.tr("Select Spreadsheet")

        ext = []
        if os.environ["ezodf"]:
            ext.append(
                self.tr("Libreoffice spreadsheet files") + " (*.ods)")
        if os.environ["xlwt"]:
            ext.append(
                self.tr("Microsoft Excel 97/2000/XP/2003 XMLL")
                + " (*.xls)")
        if os.environ["openpyxl"]:
            ext.append(
                self.tr("Microsoft Excel 2007/2010 XML") + " (*.xlsx)")

        patron = ";;".join(ext)
        folder = os.path.dirname(str(self.scene().parent().currentFilename))
        ruta = QtWidgets.QFileDialog.getSaveFileName(
            self.scene().parent(), msg, folder, patron)[0]

        if ruta:
            name, ext = os.path.splitext(ruta)
            if not ext or ext not in (".ods", ".xlsx", ".xls"):
                ruta += "." + str(ext[0]).split(".")[-1][:-1]

            if ruta[-3:] == "ods":
                templatefile = os.path.join(
                    os.environ["pychemqt"], "dat", "templates",
                    self.entity.__class__.__name__.lower()) + ".ots"

                if os.path.isfile(templatefile):
                    spreadsheet = ezodf.newdoc("ods", ruta, templatefile)
                    sheet = spreadsheet.sheets[0]
                    for attr, tipo, cell in self.entity.datamap2xls():
                        prop = self.entity._prop(attr)
                        if tipo == "value":
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

                    for i, (name, attr, unit) in enumerate(
                            self.entity.propertiesNames()):
                        value = self.entity._prop(attr)
                        txt = ""
                        if isinstance(value, list):
                            print(value)
                            # txt = value[0].text()
                        elif unit in unidades._all:
                            txt = value.text()
                            value = value.config()

                        sheet["B%i" % (i + 2)].set_value(name)
                        sheet["C%i" % (i + 2)].set_value(value)
                        sheet["D%i" % (i + 2)].set_value(txt)

                spreadsheet.save()

            elif ruta[-4:] == "xlsx":
                # TODO:
                print(ruta, "is xlsx")

            elif ruta[-3:] == "xls":
                print(ruta, "is xls")


class StreamItem(GeometricItem, QtWidgets.QGraphicsPathItem, GraphicsEntity):
    """Class to plot a mass stream"""
    up = None
    down = None
    id = 0
    free_id = []
    tipo = "stream"

    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent
        self.setPen(self._pen())
        qp = QtGui.QPainterPath()
        self.setPath(qp)
        self.setFlags(
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsFocusable)
        if StreamItem.free_id:
            self.id = StreamItem.free_id.pop(0)
        else:
            self.id = StreamItem.id + 1
            StreamItem.id += 1
        self.idLabel = TextItem("S%i" % self.id, self, selectable=False)
        self.idLabel.setZValue(2)
        self.setAcceptHoverEvents(True)
        self.tr = partial(translate, "StreamItem")

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

    def mouseDoubleClickEvent(self, event=None):
        dialog = UI_corriente.Corriente_Dialog(self.corriente)
        if dialog.exec():
            self.setCorriente(dialog.corriente)

    def copyFromProject(self):
        dialog = SelectStreamProject()
        if dialog.exec():
            indice = dialog.stream.currentText()
            data = dialog.project["stream"][indice]
            corriente = Corriente()
            corriente.readFromJSON(data)
            self.setCorriente(corriente)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key.Key_Delete or \
                event.key() == QtCore.Qt.Key.Key_Backspace:
            self.delete()
        elif event.key() == QtCore.Qt.Key.Key_Escape:
            self.setSelected(False)
        elif event.key() == QtCore.Qt.Key.Key_Return or \
                event.key() == QtCore.Qt.Key.Key_Enter:
            self.mouseDoubleClickEvent()

    def hoverEnterEvent(self, event):
        if not (self.scene().addObj and self.scene().addType == "stream"):
            self.scene().popup.populate(self.corriente, self.id, Preferences)
            point = self.mapToParent(event.scenePos())
            self.scene().popup.move(int(point.x()), int(point.y()))
            self.scene().popup.show()

    def hoverLeaveEvent(self, event):
        self.scene().popup.hide()

    def hoverMoveEvent(self, event):
        point = self.mapToParent(event.scenePos())
        self.scene().popup.move(int(point.x()), int(point.y()))
        self.scene().popup.show()

    def contextMenu(self):
        ViewAction = createAction(
            self.tr("View Properties"),
            slot=self.view, parent=self.scene())
        SolidDistributionAction = createAction(
            self.tr("Solid Distribution Fit"),
            slot=self.solidFit, parent=self.scene())
        if self.corriente:
            if not self.corriente.solido:
                SolidDistributionAction.setEnabled(False)
        else:
            ViewAction.setEnabled(False)
            SolidDistributionAction.setEnabled(False)

        contextMenu = QtWidgets.QMenu(
            "Stream %i" % self.id, self.scene().parent())
        contextMenu.setIcon(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "equipment", "stream.png")))
        contextMenu.addAction(
            self.tr("Copy from another project"),
            self.copyFromProject)
        contextMenu.addAction(SolidDistributionAction)
        contextMenu.addAction(
            self.tr("Edit"),
            self.mouseDoubleClickEvent)
        contextMenu.addAction(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "editDelete.png")),
            self.tr("Delete"),
            self.delete)
        contextMenu.addSeparator()
        contextMenu.addAction(ViewAction)
        contextMenu.addAction(
            self.tr("Export to spreadsheet"), self.exportExcel)
        contextMenu.addAction(
            self.tr("Show/Hide Id Label"),
            self.idLabelVisibility)
        contextMenu.addAction(
            self.tr("Appearance"),
            self.format)
        return contextMenu

    def redraw(self, entrada=None, salida=None):
        """Recalcule the stream path in screen"""
        if entrada:
            self.entrada = entrada
        if salida:
            self.salida = salida

        max_height = max(self.up.boundingRect().height(),
                         self.down.boundingRect().height())
        max_width = max(self.up.boundingRect().width(),
                        self.down.boundingRect().width())
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
                if abs(self.entrada.y() - self.salida.y()) > max_height:
                    # cabe la linea entre ambos equipos
                    y_mean = (y_up + y_down + height_sup) / 2.
                    qp.lineTo(QtCore.QPointF(
                        self.entrada.x() + Xdist_entrada,
                        self.entrada.y() + Ydist_entrada))
                    qp.lineTo(QtCore.QPointF(
                        self.entrada.x() + Xdist_entrada, y_mean))
                    qp.lineTo(QtCore.QPointF(
                        self.salida.x() + Xdist_salida, y_mean))
                    qp.lineTo(QtCore.QPointF(
                        self.salida.x() + Xdist_salida,
                        self.salida.y() + Ydist_entrada))

                else:
                    # sacamos la linea por encima de los equipos
                    y_mean = y_up - 20
                    qp.lineTo(QtCore.QPointF(
                        self.entrada.x() + Xdist_entrada,
                        self.entrada.y() + Ydist_entrada))
                    qp.lineTo(QtCore.QPointF(
                        self.entrada.x() + Xdist_entrada, y_mean))
                    qp.lineTo(QtCore.QPointF(
                        self.salida.x() + Xdist_salida, y_mean))
                    qp.lineTo(QtCore.QPointF(
                        self.salida.x() + Xdist_salida,
                        self.salida.y() + Ydist_salida))
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
            self.idLabel.setPos(x_mean, (self.entrada.y()+self.salida.y())/2-10)

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


class EquipmentItem(QtSvgWidgets.QGraphicsSvgItem, GraphicsEntity):
    """Class to plot equipment item"""
    up = []
    down = []
    up_used = 0
    down_used = 0
    id = 0
    id_in = 0
    id_out = 0
    tipo = "equip"

    def __init__(self, name, dialogoId, parent=None):
        self.name = name
        imagen = os.environ["pychemqt"] + "images/equipment/%s.svg" % name
        super().__init__(imagen, parent=parent)
        self.dialogoId = dialogoId
        self.setFlags(
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsMovable
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges
            | QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIsFocusable)
        self.imagen = imagen
        self.angle = 0
        self.setAcceptHoverEvents(True)

        if dialogoId is not None:
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
                obj.setBrush(QtGui.QColor(
                    Preferences.get("PFD", 'Color_Entrada')))
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
                obj.setBrush(QtGui.QColor(
                    Preferences.get("PFD", 'Color_Salida')))
                self.output.append(obj)
        self.showInput(False)

    @property
    def equipment(self):
        return self.scene().project.getItem(self.id)

    @property
    def entity(self):
        return self.equipment

    def mouseDoubleClickEvent(self, event=None):
        if self.dialogoId is not None:
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
            if dialog.exec():
                self.scene().project.setItem(self.id, dialog.Equipment)
                # self.up[0].setCorriente(dialog.Equipment.entrada)
                for i, corriente in enumerate(dialog.Equipment.salida):
                    self.down[i].setCorriente(corriente)

        else:
            if self.output:
                self.down[0].mouseDoubleClickEvent()
            else:
                self.up[0].mouseDoubleClickEvent()

    def mousePressEvent(self, event):
        QtSvgWidgets.QGraphicsSvgItem.mousePressEvent(self, event)
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
                self.scene().Pos.append(
                    self.mapToScene(punto[x].rect().center()))
                self.scene().points.append(punto[x])
            else:
                self.scene().Pos.append(event.pos())

    def mouseMoveEvent(self, event=None):
        if event:
            QtWidgets.QGraphicsPixmapItem.mouseMoveEvent(self, event)
        for i, corriente in enumerate(self.up):
            corriente.redraw(
                salida=self.mapToScene(self.input[i].rect().center()))
        for i, corriente in enumerate(self.down):
            corriente.redraw(
                entrada=self.mapToScene(self.output[i].rect().center()))

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
                self.scene().popup.populate(self.equipment, self.id, Preferences)
            else:
                if self.output:
                    self.scene().popup.populate(
                        self.down[0].corriente, self.down[0].id, Preferences)
                else:
                    self.scene().popup.populate(
                        self.up[0].corriente, self.up[0].id, Preferences)
            point = event.scenePos()
            self.scene().popup.move(int(point.x()), int(point.y()))
            self.scene().popup.show()

    def hoverLeaveEvent(self, event):
        self.showInput(False)
        self.scene().popup.hide()

    def hoverMoveEvent(self, event):
        point = event.scenePos()
        self.scene().popup.move(int(point.x()), int(point.y()))

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key.Key_Delete or \
                event.key() == QtCore.Qt.Key.Key_Backspace:
            self.delete()
        elif event.key() == QtCore.Qt.Key.Key_Escape:
            self.setSelected(False)
        elif event.key() == QtCore.Qt.Key.Key_Return or \
                event.key() == QtCore.Qt.Key.Key_Enter:
            self.mouseDoubleClickEvent()
        elif event.key() == QtCore.Qt.Key.Key_Up:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()-factor))
            self.mouseMoveEvent()
        elif event.key() == QtCore.Qt.Key.Key_Down:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()+factor))
            self.mouseMoveEvent()
        elif event.key() == QtCore.Qt.Key.Key_Left:
            self.setPos(QtCore.QPointF(self.pos().x()-factor, self.pos().y()))
            self.mouseMoveEvent()
        elif event.key() == QtCore.Qt.Key.Key_Right:
            self.setPos(QtCore.QPointF(self.pos().x()+factor, self.pos().y()))
            self.mouseMoveEvent()

    def itemChange(self, key, value):
        if self.scene() and key == \
                QtWidgets.QGraphicsItem.GraphicsItemChange.ItemPositionChange:
            mainwindow = self.scene().parent()
            if mainwindow.dirty[mainwindow.idTab] == False:
                mainwindow.dirty[mainwindow.idTab] = True
                mainwindow.activeControl(True)
                mainwindow.tabModified(mainwindow.idTab)
        return QtWidgets.QGraphicsItem.itemChange(self, key, value)

    def contextMenu(self):
        if self.dialogoId != None:
            ViewAction = createAction(
                self.tr("View Properties"),
                slot=self.view, parent=self.scene())
            ViewAction.setEnabled(self.equipment.status)

            contextMenu = QtWidgets.QMenu(
                "Equipment %i" % self.id, self.scene().parent())
            contextMenu.setIcon(QtGui.QIcon(
                os.path.join(IMAGE_PATH, "equipment", "%s.svg" % self.name)))
            contextMenu.addAction(self.tr("Edit"), self.mouseDoubleClickEvent)
            contextMenu.addAction(
                QtGui.QIcon(os.path.join(
                    IMAGE_PATH, "button", "editDelete.png")),
                self.tr("Delete"),
                self.delete)
            contextMenu.addSeparator()
            contextMenu.addAction(ViewAction)
            contextMenu.addAction(
                self.tr("Export to spreadsheet"), self.exportExcel)
            contextMenu.addAction(
                self.tr("Show/Hide Id Label"), self.idLabelVisibility)

        # contextMenu.addAction(
        # self.tr("Appearance"), self.format)
        #            contextMenu.addSeparator()
        #            contextMenu.addAction("Run", self.mouseDoubleClickEvent)

        else:
            if self.output and self.down:
                contextMenu = self.down[0].contextMenu()
            elif not self.output and self.up:
                contextMenu = self.up[0].contextMenu()
            else:
                contextMenu = QtWidgets.QMenu()

        self.menuTransform = QtWidgets.QMenu(self.tr("Transform"))
        self.menuTransform.addAction(
            QtGui.QIcon(os.path.join(
                IMAGE_PATH, "button", "transform_rotate_90.png")),
            self.tr("Rotate by 90"),
            partial(self.rotate, 90))
        self.menuTransform.addAction(
            QtGui.QIcon(os.path.join(
                IMAGE_PATH, "button", "transform_rotate_180.png")),
            self.tr("Rotate by 180"),
            partial(self.rotate, 180))
        self.menuTransform.addAction(
            QtGui.QIcon(os.path.join(
                IMAGE_PATH, "button", "transform_rotate_270.png")),
            self.tr("Rotate by 270"),
            partial(self.rotate, 270))
        self.menuTransform.addSeparator()
        self.menuTransform.addAction(
            self.tr("Mirror about X"),
            partial(self.rotate, 270))
        self.menuTransform.addAction(
            self.tr("Mirror about Y"),
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


class SelectStreamProject(QtWidgets.QDialog):
    """Dialog to select a stream from any external project"""
    project = None

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Select stream from file"))

        layout = QtWidgets.QVBoxLayout(self)
        label = self.tr("Project path")
        msg = self.tr("Select pychemqt project file")
        patrones = []
        patrones.append(self.tr("pychemqt project file") + " (*.pcq)")
        patron = ";;".join(patrones)
        self.filename = PathConfig(label + ":", msg=msg, patron=patron)
        self.filename.valueChanged.connect(self.changeproject)
        layout.addWidget(self.filename)

        lyt1 = QtWidgets.QHBoxLayout()
        lyt1.addWidget(QtWidgets.QLabel(self.tr("Streams")))
        self.stream = QtWidgets.QComboBox()
        lyt1.addWidget(self.stream)
        lyt1.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        layout.addLayout(lyt1)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding))

        lyt2 = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel()
        lyt2.addWidget(self.status)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttonBox.button(
            QtWidgets.QDialogButtonBox.StandardButton.Ok).setEnabled(False)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        lyt2.addWidget(self.buttonBox)
        layout.addLayout(lyt2)

    def changeproject(self, path):
        """Upgrade dialog with the streams of loaded file"""
        self.stream.clear()
        st = self.tr("Loading project...")
        self.status.setText(st)
        QtWidgets.QApplication.processEvents()
        try:
            with open(path, "r") as file:
                self.project = json.load(file)
        except Exception as e:
            self.status.setText(self.tr("Failed to loading project..."))
            raise e
        self.buttonBox.button(
            QtWidgets.QDialogButtonBox.StandardButton.Ok).setEnabled(True)
        self.status.setText(self.tr("Project loaded succesfully"))
        for stream in sorted(self.project["stream"].keys()):
            self.stream.addItem(stream)


class TextItemDlg(QtWidgets.QDialog):
    """Dialog to edit texts in PFD"""

    def __init__(self, text=None, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        self.editor = texteditor.TextEditor()
        self.editor.notas.textChanged.connect(self.updateUi)
        layout.addWidget(self.editor, 1, 1, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 2, 1, 1, 1)
        self.editor.notas.setFocus()
        if text:
            self.editor.setText(text)
        self.setWindowTitle(self.tr("Edit text"))
        self.updateUi()

    def updateUi(self):
        """Set enable/disable OK button"""
        self.buttonBox.button(
            QtWidgets.QDialogButtonBox.StandardButton.Ok).setEnabled(
                bool(self.editor.notas.toPlainText()))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    #    dialogo = ConfLineDialog()
    dialogo = SelectStreamProject()
    dialogo.show()
    sys.exit(app.exec())
