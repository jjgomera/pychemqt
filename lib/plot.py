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


###############################################################################
# Plot library with all matplotlib functionality
###############################################################################


from functools import partial
import os

from matplotlib import style
from matplotlib.backends import backend_qtagg
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties

from lib.config import Preferences
from tools.qt import QtCore, QtGui, QtWidgets, translate
from UI.widgets import (ColorSelector, Entrada_con_unidades, InputFont,
                        LineStyleCombo, MarkerCombo, createAction)


RCParams = {
    # Dict with options, the keys are the keyword matplotlib
    # Each option must define:
    #   - tooltip for widget with a tiny explanation
    #   - Widget class to use
    #   - Optional parameters
    #       - QComboBox: options to populate widget
    #       - QSpinBox: two values for define range (default 0-1)
    #       - QDoubleSpinBox: two values for define range (default 0-1),
    #           one more to define step (default 0.01 for range 0-1 or 1 in
    #           other cases) and one more for define decimals (default 2)

    "axes.facecolor": (
        translate("plot", "Axes background color"), ColorSelector),
    "axes.edgecolor": (translate("plot", "Axes edge color"), ColorSelector),
    "axes.linewidth": (
        translate("plot", "Edge line width"), QtWidgets.QDoubleSpinBox,
        0, 5, 0.1, 1),
    "axes.grid": (
        translate("plot", "Display grid or not"), QtWidgets.QCheckBox),
    "axes.grid.axis": (
        translate("plot", "Which axis the grid should apply to"),
        QtWidgets.QComboBox, "both", "x", "y"),
    "axes.grid.which": (
        translate("plot", "Grid lines at {major, minor, both} ticks"),
        QtWidgets.QComboBox, "both", "major", "minor"),
    "axes.titlelocation": (
        translate("plot", "Alignment of the title"), QtWidgets.QComboBox,
        "left", "right", "center"),
    "axes.titlesize": (
        translate("plot", "Font size of the axes title"),
        QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium',
        'large', 'x-large', 'xx-large'),
    "axes.titleweight": (
        translate("plot", "Font weight of title"),
        QtWidgets.QComboBox, "normal", "bold"),
    "axes.titlecolor": (
        translate("plot", "Color of axes title"), ColorSelector),
    "axes.titley": (
        translate("plot", "Position title (axes relative units)"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "axes.titlepad": (
        translate("plot", "Pad between axes and title in points"),
        QtWidgets.QSpinBox, 0, 20),
    "axes.labelsize": (
        translate("plot", "Font size of the x and y labels"),
        QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium',
        'large', 'x-large', 'xx-large'),
    "axes.labelpad": (
        translate("plot", "Pad between label and axis"),
        QtWidgets.QSpinBox, 0, 20),
    "axes.labelweight": (
        translate("plot", "Weight of the x and y labels"),
        QtWidgets.QComboBox, "normal", "bold"),
    "axes.labelcolor": (
        translate("plot", "Axes label color"), ColorSelector),
    "axes.axisbelow": (
        translate("plot", "Draw axis gridlines and ticks"),
        QtWidgets.QComboBox, "True", "line", "False"),
    "axes.formatter.limits": (translate(
        "plot", "Use scientific notation if log10 of the axis range "
        "is larger than this value"),
        QtWidgets.QSpinBox, 0, 10),
    "axes.formatter.use_locale": (translate(
        "plot", "Format tick labels according to the user's locale"),
        QtWidgets.QCheckBox),
    "axes.formatter.use_mathtext": (
        translate("plot", "Use mathtext for scientific notation"),
        QtWidgets.QCheckBox),
    "axes.formatter.min_exponent": (
        translate("plot", "Minimum exponent to format in scientific notation"),
        QtWidgets.QSpinBox, 0, 10),
    "axes.formatter.useoffset": (translate(
        "plot", "The tick label formatter will default to labeling "
        "ticks relative to an offset when the data range is small compared to "
        "the minimum absolute value of the data."), QtWidgets.QCheckBox),
    "axes.formatter.offset_threshold": (translate(
        "plot", "When useoffset is True, the offset will be used when it can "
        "remove at least this number of significant digits from tick labels"),
        QtWidgets.QSpinBox, 0, 10),
    "axes.spines.left": (
        translate("plot", "Display axis spines"), QtWidgets.QCheckBox),
    "axes.spines.bottom": (
        translate("plot", "Display axis spines"), QtWidgets.QCheckBox),
    "axes.spines.top": (
        translate("plot", "Display axis spines"), QtWidgets.QCheckBox),
    "axes.spines.right": (
        translate("plot", "Display axis spines"), QtWidgets.QCheckBox),
    "axes.unicode_minus": (translate(
        "plot", "Use Unicode for the minus symbol rather than hyphen"),
        QtWidgets.QCheckBox),
    "axes.xmargin": (
        translate("plot", "X margin"), QtWidgets.QDoubleSpinBox, 0, 1),
    "axes.ymargin": (
        translate("plot", "Y margin"), QtWidgets.QDoubleSpinBox, 0, 1),
    "axes.zmargin": (
        translate("plot", "Z margin"), QtWidgets.QDoubleSpinBox, 0, 1),

    "axes3d.grid": (
        translate("plot", "Display grid on 3D axes"), QtWidgets.QCheckBox),

    "figure.titlesize": (
        translate("plot", "Size of the figure title"),
        QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium',
        'large', 'x-large', 'xx-large'),
    "figure.titleweight": (
        translate("plot", "Weight of the figure title"),
        QtWidgets.QComboBox, "normal", "bold"),
    "figure.labelsize": (
        translate("plot", "Size of the figure label"),
        QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium',
        'large', 'x-large', 'xx-large'),
    "figure.labelweight": (
        translate("plot", "Weight of figure label"),
        QtWidgets.QComboBox, "normal", "bold"),
    # ("figure.figsize": (,     6.4, 4.8  # figure size in inches
    "figure.dpi": (
        translate("plot", "Figure dots per inch"),
        QtWidgets.QSpinBox, 10, 200),
    "figure.facecolor": (
        translate("plot", "Figure face color"), ColorSelector),
    "figure.edgecolor": (
        translate("plot", "Figure edge color"), ColorSelector),
    "figure.frameon": (
        translate("plot", "Enable figure frame"), QtWidgets.QCheckBox),
    "figure.subplot.left": (
        translate("plot", "The left side of the subplots of the figure"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.subplot.right":
        (translate("plot", "The right side of the subplots of the figure"),
         QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.subplot.bottom": (
        translate("plot", "The bottom of the subplots of the figure"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.subplot.top": (
        translate("plot", "The top of the subplots of the figure"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.subplot.wspace": (
        translate("plot", "Width reserved for space between subplots"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.subplot.hspace": (
        translate("plot", "Height reserved for space between subplots"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.autolayout": (
        translate("plot", "Automatically adjust subplot"),
        QtWidgets.QCheckBox),
    "figure.constrained_layout.h_pad": (
        translate("plot", "Padding around axes objects. Float representing"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.constrained_layout.w_pad": (
        translate("plot", "Padding around axes objects. Float representing"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.constrained_layout.hspace": (
        translate("plot", "Space between subplot groups. Float representing"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "figure.constrained_layout.wspace": (
        translate("plot", "Space between subplot groups. Float representing"),
        QtWidgets.QDoubleSpinBox, 0, 1),

    "font.family": (
        "", QtWidgets.QComboBox,
        "serif", "sans-serif", "cursive", "fantasy", "monospace"),
    "font.style": ("", QtWidgets.QComboBox, "normal", "italic", "oblique"),
    "font.variant": ("", QtWidgets.QComboBox, "normal", "small-caps"),
    "font.weight": ("", QtWidgets.QComboBox, "normal", "bold"),
    "font.stretch": (
        "", QtWidgets.QComboBox, "ultra-condensed", "extra-condensed",
        "condensed", "semi-condensed", "normal", "semi-expanded", "expanded",
        "extra-expanded", "ultra-expanded", "wider", "narrower"),
    "font.size": ("", QtWidgets.QDoubleSpinBox, 5, 20, 0.5, 1),

    "grid.color": (translate("plot", "Grid color"), ColorSelector),
    "grid.linestyle": (translate("plot", "Grid line style"), LineStyleCombo),
    "grid.linewidth": (
        translate("plot", "Grid line width in points"),
        QtWidgets.QDoubleSpinBox, 0, 5, 0.1, 1),
    "grid.alpha": (
        translate("plot", "Grid lines transparency"),
        QtWidgets.QDoubleSpinBox, 0, 1),

    "hatch.linewidth": (
        translate("plot", "Line width in points"),
        QtWidgets.QDoubleSpinBox, 0, 5, 0.1, 1),
    "hatch.color": (translate("plot", "Hatch color"), ColorSelector),

    "legend.loc": (
        translate("plot", "Location of legend in axes"), QtWidgets.QComboBox,
        'best', 'center', 'center left', 'center right', 'lower center',
        'lower left', 'lower right', 'right', 'upper center', 'upper left',
        'upper right'),
    "legend.frameon": (
        translate("plot", "Draw the legend on a background patch"),
        QtWidgets.QCheckBox),
    "legend.framealpha": (
        translate("plot", "Legend patch transparency"),
        QtWidgets.QDoubleSpinBox, 0, 1),
    "legend.facecolor": (
        translate("plot", "Legend patch color"), ColorSelector),
    "legend.edgecolor": (
        translate("plot", "Background patch boundary color"), ColorSelector),
    "legend.fancybox": (translate(
        "plot",
        "Use a rounded box for the legend background, else a rectangle"),
        QtWidgets.QCheckBox),
    "legend.shadow": (
        translate("plot", "Give background a shadow effect"),
        QtWidgets.QCheckBox),
    "legend.numpoints": (
        translate("plot", "Number of marker points in the legend line"),
        QtWidgets.QSpinBox, 1, 10),
    "legend.scatterpoints": (
        translate("plot", "Number of scatter points"),
        QtWidgets.QSpinBox, 1, 10),
    "legend.markerscale": (
        translate("plot", "Relative size of legend markers vs. original"),
        QtWidgets.QDoubleSpinBox, 0, 2, 0.1, 1),
    "legend.fontsize": (
        "", QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium',
        'large', 'x-large', 'xx-large'),
    "legend.labelcolor": ("", ColorSelector),
    "legend.title_fontsize": (
        "", QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium',
        'large', 'x-large', 'xx-large'),
    "legend.borderpad": (translate(
        "plot", "Dimensions as fraction of font size for border whitespace"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "legend.labelspacing": (translate(
        "plot", "Dimensions as fraction of font size for the vertical space "
        "between the legend entries"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "legend.handlelength": (translate(
        "plot", "Dimensions as fraction of font size for the length of the "
        "legend lines"), QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "legend.handleheight": (translate(
        "plot", "Dimensions as fraction of font size for the height of the "
        "legend handle"), QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "legend.handletextpad": (translate(
        "plot", "Dimensions as fraction of font size for the space between "
        "the legend line and legend text"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "legend.borderaxespad": (translate(
        "plot", "Dimensions as fraction of font size for the border between "
        "the axes and legend edge"), QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "legend.columnspacing": (translate(
        "plot", "Dimensions as fraction of font size for column separation"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),

    "lines.linewidth": (
        translate("plot", "Line width in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "lines.linestyle": (translate("plot", "Line style"), LineStyleCombo),
    "lines.color": (translate("plot", "Line color"), ColorSelector),
    "lines.marker": (translate("plot", "Marker style"), MarkerCombo),
    "lines.markerfacecolor": (
        translate("plot", "Marker face color"), ColorSelector),
    "lines.markeredgecolor": (
        translate("plot", "Marker edge color"), ColorSelector),
    "lines.markeredgewidth": (
        translate("plot", "Line width around the marker symbol"),
        QtWidgets.QDoubleSpinBox, 0, 5, 0.1, 1),
    "lines.markersize": (
        translate("plot", "Marker size, in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "lines.dash_joinstyle": (
        "", QtWidgets.QComboBox, "miter", "round", "bevel"),
    "lines.dash_capstyle": (
        "", QtWidgets.QComboBox, "butt", "round", "projecting"),
    "lines.solid_joinstyle": (
        "", QtWidgets.QComboBox, "miter", "round", "bevel"),
    "lines.solid_capstyle": (
        "", QtWidgets.QComboBox, "butt", "round", "projecting"),
    "lines.antialiased": (
        translate("plot", "Render antialiased"), QtWidgets.QCheckBox),
    "lines.scale_dashes": ("", QtWidgets.QCheckBox),

    "patch.linewidth": (
        translate("plot", "Edge width in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "patch.facecolor": (translate("plot", "Patch face color"), ColorSelector),
    "patch.edgecolor": (translate("plot", "Patch edge color"), ColorSelector),
    "patch.force_edgecolor": (
        translate("plot", "Always use edgecolor"), QtWidgets.QCheckBox),
    "patch.antialiased": (
        translate("plot", "Render patch antialiased"), QtWidgets.QCheckBox),

    "xtick.top": (
        translate("plot", "Draw ticks on the top side"), QtWidgets.QCheckBox),
    "xtick.bottom": (
        translate("plot", "Draw ticks on the bottom side"),
        QtWidgets.QCheckBox),
    "xtick.labeltop": (
        translate("plot", "Draw label on the top"), QtWidgets.QCheckBox),
    "xtick.labelbottom": (
        translate("plot", "Draw label on the bottom"), QtWidgets.QCheckBox),
    "xtick.major.size": (
        translate("plot", "Major tick size in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "xtick.minor.size": (
        translate("plot", "Minor tick size in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "xtick.major.width": (
        translate("plot", "Major tick width in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "xtick.minor.width": (
        translate("plot", "Minor tick width in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "xtick.major.pad": (
        translate("plot", "Distance to major tick label in points"),
        QtWidgets.QDoubleSpinBox, 0, 20, 0.1, 1),
    "xtick.minor.pad": (
        translate("plot", "Distance to the minor tick label in points"),
        QtWidgets.QDoubleSpinBox, 0, 20, 0.1, 1),
    "xtick.color": (translate("plot", "Color of the ticks"), ColorSelector),
    "xtick.labelcolor": (translate(
        "plot", "Color of the tick labels or inherit from xtick.color"),
        ColorSelector),
    "xtick.labelsize": (
        translate("plot", "Font size of the tick labels"),
        QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium', 'large',
        'x-large', 'xx-large'),
    "xtick.direction": (
        translate("plot", "Direction"), QtWidgets.QComboBox,
        "in", "out", "inout"),
    "xtick.minor.visible": (
        translate("plot", "Visibility of minor ticks on x-axis"),
        QtWidgets.QCheckBox),
    "xtick.major.top": (
        translate("plot", "Draw x axis top major ticks"), QtWidgets.QCheckBox),
    "xtick.major.bottom": (
        translate("plot", "Draw x axis bottom major ticks"),
        QtWidgets.QCheckBox),
    "xtick.minor.top": (
        translate("plot", "Draw x axis top minor ticks"), QtWidgets.QCheckBox),
    "xtick.minor.bottom": (
        translate("plot", "Draw x axis bottom minor ticks"),
        QtWidgets.QCheckBox),
    "xtick.alignment": (
        translate("plot", "Alignment of ticks"), QtWidgets.QComboBox,
        "left", "center", "right"),

    "ytick.left": (
        translate("plot", "Draw ticks on the left side"), QtWidgets.QCheckBox),
    "ytick.right": (
        translate("plot", "Draw ticks on the right side"),
        QtWidgets.QCheckBox),
    "ytick.labelleft": (
        translate("plot", "Draw label on the left"), QtWidgets.QCheckBox),
    "ytick.labelright": (
        translate("plot", "Draw label on the right"), QtWidgets.QCheckBox),
    "ytick.major.size": (
        translate("plot", "Major tick size in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "ytick.minor.size": (
        translate("plot", "Minor tick size in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "ytick.major.width": (
        translate("plot", "Major tick width in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "ytick.minor.width": (
        translate("plot", "Minor tick width in points"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "ytick.major.pad": (
        translate("plot", "Distance to major tick label in points"),
        QtWidgets.QDoubleSpinBox, 0, 20, 0.1, 1),
    "ytick.minor.pad": (
        translate("plot", "Distance to the minor tick label in points"),
        QtWidgets.QDoubleSpinBox, 0, 20, 0.1, 1),
    "ytick.color": (translate("plot", "Color of the ticks"), ColorSelector),
    "ytick.labelcolor": (translate(
        "plot", "Color of the tick labels or inherit from ytick.color"),
        ColorSelector),
    "ytick.labelsize": (
        translate("plot", "Font size of the tick labels"),
        QtWidgets.QComboBox, 'xx-small', 'x-small', 'small', 'medium', 'large',
        'x-large', 'xx-large'),
    "ytick.direction": (
        translate("plot", "Direction"), QtWidgets.QComboBox,
        "in", "out", "inout"),
    "ytick.minor.visible": (
        translate("plot", "Visibility of minor ticks on y-axis"),
        QtWidgets.QCheckBox),
    "ytick.major.left": (
        translate("plot", "Draw y axis left major ticks"),
        QtWidgets.QCheckBox),
    "ytick.major.right": (
        translate("plot", "Draw y axis right major ticks"),
        QtWidgets.QCheckBox),
    "ytick.minor.left": (
        translate("plot", "Draw y axis left minor ticks"),
        QtWidgets.QCheckBox),
    "ytick.minor.right": (
        translate("plot", "Draw y axis right minor ticks"),
        QtWidgets.QCheckBox),
    "ytick.alignment": (
        translate("plot", "Alignment of ticks"), QtWidgets.QComboBox,
        'bottom', 'baseline', 'center', 'center_baseline', 'top'),

    "xaxis.labellocation": (
        translate("plot", "Alignment of the xaxis label"),
        QtWidgets.QComboBox, "center", "left", "right"),
    "yaxis.labellocation": (
        translate("plot", "Alignment of the yaxis label"),
        QtWidgets.QComboBox, "center", "bottom", "top"),

    "savefig.dpi": (
        translate("plot", "Figure dots per inch"), QtWidgets.QSpinBox, 0, 100),
    "savefig.facecolor": (
        translate("plot", "Figure facecolor when saving"), ColorSelector),
    "savefig.edgecolor": (
        translate("plot", "Figure edgecolor when saving"), ColorSelector),
    "savefig.format": (
        translate("plot", "File format so save"), QtWidgets.QComboBox,
        "png", "ps", "pdf", "svg"),
    "savefig.bbox": ("", QtWidgets.QComboBox, "standard", "tight"),
    "savefig.pad_inches": (translate(
        "plot", "Padding to be used, when bbox is set to 'tight'"),
        QtWidgets.QDoubleSpinBox, 0, 10, 0.1, 1),
    "savefig.transparent": (
        translate("plot", "Figures are saved with a transparent background"),
        QtWidgets.QCheckBox)}


class PlotWidget(backend_qtagg.FigureCanvasQTAgg):
    """QWidget with matplotlib integration"""
    def __init__(self, dim=2, width=15, height=5, dpi=100, parent=None):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)

        self.dim = dim
        self.setParent(parent)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding)

        if dim == 2:
            self.ax = self.fig.add_subplot(111)

        elif dim == 3:
            self.ax = self.fig.add_subplot(projection="3d")
            self.ax.mouse_init(rotate_btn=1, zoom_btn=2)
        # In any other case the axes are not defined and can be configured at
        # used site

    def plot(self, *args, **kwargs):
        """Direct accesst to ax plot procedure"""
        self.ax.plot(*args, **kwargs)

    def savePNG(self):
        """Save chart image to png file"""
        fmt = "Portable Network Graphics (*.png)"
        fname, ext = QtWidgets.QFileDialog.getSaveFileName(
            self, translate("Save chart to file"), "./", fmt)
        if fname and ext == fmt:
            if fname.split(".")[-1] != "png":
                fname += ".png"
            self.fig.savefig(fname, facecolor='#fafafa')

    # def plot_3D(self, labels, xdata, ydata, zdata, config=None):
        # """Método que dibuja la matriz de datos"""
        # self.ax.clear()
        # self.data = {"x": xdata[0], "y": ydata[:, 0], "z": zdata}

        # if config and config.getboolean("MEOS", "surface"):
            # self.ax.plot_surface(xdata, ydata, zdata, rstride=1, cstride=1)
        # else:
            # self.ax.plot_wireframe(xdata, ydata, zdata, rstride=1, cstride=1)

        # self.ax.set_xlabel(labels[0])
        # self.ax.set_ylabel(labels[1])
        # self.ax.set_zlabel(labels[2])
        # self.ax.mouse_init(rotate_btn=1, zoom_btn=2)


class PlotDialog(QtWidgets.QDialog):
    """QDialog including Plotwidget, navigationtoolbar and a button to close"""
    def __init__(self, accept=False, cancel=True, parent=None):
        super().__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)

        self.plot = PlotWidget()
        gridLayout.addWidget(self.plot, 1, 1, 1, 2)
        self.toolbar = backend_qtagg.NavigationToolbar2QT(self.plot, self.plot)
        gridLayout.addWidget(self.toolbar, 2, 1)
        btonBox = QtWidgets.QDialogButtonBox()
        if accept:
            btonBox.addButton(QtWidgets.QDialogButtonBox.StandardButton.Ok)
            btonBox.accepted.connect(self.accept)
        if cancel:
            btonBox.addButton(QtWidgets.QDialogButtonBox.StandardButton.Cancel)
            btonBox.rejected.connect(self.reject)
        btonBox.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum,
                              QtWidgets.QSizePolicy.Policy.Maximum)
        gridLayout.addWidget(btonBox, 2, 2)

    def addText(self, *args, **kwargs):
        """Direct access to ax text procedure"""
        self.plot.ax.text(*args, **kwargs)

    def addData(self, *args, **kwargs):
        """Direct access to ax plot procedure"""
        self.plot.ax.plot(*args, **kwargs)


# Load style defined in preferences
if Preferences.getboolean("Plot", 'customize'):
    rc = {}
    for key, (tip, widget, *args) in RCParams.items():
        if widget == QtWidgets.QDoubleSpinBox:
            rc[key] = Preferences.getfloat("Plot", key)
        elif widget == QtWidgets.QSpinBox:
            if key == "axes.formatter.limits":
                rc[key] = Preferences.get("Plot", key)
            else:
                rc[key] = Preferences.getint("Plot", key)
        elif widget == QtWidgets.QCheckBox:
            rc[key] = Preferences.getboolean("Plot", key)
        else:
            rc[key] = Preferences.get("Plot", key)
    style.use(rc)
elif Preferences.getint("Plot", 'style') == 0:
    style.use("default")
else:
    style.use(style.available[Preferences.getint("Plot", 'style')-1])


if __name__ == "__main__":
    import os
    import sys
    from configparser import ConfigParser
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    app = QtWidgets.QApplication(sys.argv)

    conf = ConfigParser()
    conf.read(conf_dir+"pychemqtrc")
    dialogo = ConfPlot(conf)

    dialogo.show()
    sys.exit(app.exec())
