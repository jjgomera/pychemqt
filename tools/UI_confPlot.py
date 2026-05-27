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

Dialog for configuration of matplotlib plot widget

'''

from matplotlib import rcParams, rcdefaults, style
from matplotlib.colors import to_rgb
from numpy import arange

from lib.plot import PlotWidget, RCParams
from tools.qt import QtWidgets
from UI.widgets import ColorSelector, LineStyleCombo, MarkerCombo


class ConfPlot(QtWidgets.QDialog):
    """Matplotlib configuration"""

    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(self.tr("Matplotlib Style:")), 1, 1)
        self.style = QtWidgets.QComboBox()
        layout.addWidget(self.style, 1, 2)
        self.style.addItem("default")
        for sty in style.available:
            self.style.addItem(sty)
        self.style.currentTextChanged.connect(self.updateStyle)

        self.customize = QtWidgets.QCheckBox(self.tr("Costomize Style:"))
        layout.addWidget(self.customize, 2, 1, 1, 3)
        self.tabRcParams = QtWidgets.QTabWidget()
        self.tabRcParams.setEnabled(False)
        layout.addWidget(self.tabRcParams, 3, 1, 1, 3)
        self.customize.toggled.connect(self.tabRcParams.setEnabled)
        self.customize.toggled.connect(self.updatePlot)

        tab = []
        tab_widgets = []
        for label, (tooltip, wdg, *opt) in sorted(RCParams.items()):
            title = label.split(".")[0]

            # Regroup tiny option in axes tab
            if title in ("xaxis", "yaxis", "axes3d"):
                title = "axes"

            # Add new tab if not exist
            if title not in tab:
                tab.append(title)
                tab_widgets.append(QtWidgets.QWidget())
                QtWidgets.QVBoxLayout(tab_widgets[-1])

            # Get index of current tab to add widget
            idx = tab.index(title)

            # Create widget with a horizontal layout to add label and
            # appropiate input widget
            rc_widget = QtWidgets.QWidget()
            lyt_wdg = QtWidgets.QHBoxLayout(rc_widget)
            lyt_wdg.addWidget(QtWidgets.QLabel(label, rc_widget))

            # Create the appropiate input widget as defined in rc dict
            if wdg == QtWidgets.QDoubleSpinBox:
                wdg = QtWidgets.QDoubleSpinBox()
                wdg.setRange(*opt[:2])
                if len(opt) >= 3:
                    wdg.setSingleStep(opt[2])
                if len(opt) == 4:
                    wdg.setDecimals(opt[3])
                elif opt[1] == 1 and opt[0] == 0:
                    wdg.setSingleStep(0.01)
                wdg.valueChanged.connect(self.updatePlot)
            elif wdg == QtWidgets.QSpinBox:
                wdg = QtWidgets.QSpinBox()
                wdg.setRange(*opt)
                wdg.valueChanged.connect(self.updatePlot)
            elif wdg == QtWidgets.QComboBox:
                wdg = QtWidgets.QComboBox()
                for value in opt:
                    wdg.addItem(value)
                wdg.currentIndexChanged.connect(self.updatePlot)
            elif wdg == LineStyleCombo:
                wdg = LineStyleCombo()
                wdg.currentIndexChanged.connect(self.updatePlot)
            elif wdg == MarkerCombo:
                wdg = MarkerCombo()
                wdg.currentIndexChanged.connect(self.updatePlot)
            elif wdg == ColorSelector:
                wdg = ColorSelector()
                wdg.valueChanged.connect(self.updatePlot)
            elif wdg == QtWidgets.QCheckBox:
                wdg = QtWidgets.QCheckBox()
                wdg.toggled.connect(self.updatePlot)

            wdg.setToolTip(tooltip)
            lyt_wdg.addWidget(wdg)
            lyt_wdg.addItem(QtWidgets.QSpacerItem(
                10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Fixed))

            tab_widgets[idx].layout().addWidget(rc_widget)

        for title, wdg in zip(tab, tab_widgets):
            # Add a spacer item at end of each tabwidget
            wdg.layout().addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding))

            # Now when the widget is pupulate add the scrollArea
            scroll = QtWidgets.QScrollArea()
            scroll.setFrameStyle(QtWidgets.QFrame.Shape.NoFrame)
            scroll.setWidget(wdg)
            self.tabRcParams.addTab(scroll, title)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 3)

        if config.has_section("Plot"):
            self.style.setCurrentIndex(config.getint("Plot", 'style'))
            self.customize.setChecked(config.getboolean("Plot", 'customize'))

            for idx in range(self.tabRcParams.count()):
                tab = self.tabRcParams.widget(idx).widget()
                for parent_wdg in tab.children()[1:]:
                    label, wdg = parent_wdg.children()[1:]
                    if isinstance(wdg, QtWidgets.QDoubleSpinBox):
                        wdg.setValue(config.getfloat("Plot", label.text()))
                    elif isinstance(wdg, QtWidgets.QSpinBox):
                        if label.text() == "axes.formatter.limits":
                            # Convert saved value as tuple to integer
                            value = config.get("Plot", label.text())
                            value = -int(value.split(",")[0])
                            wdg.setValue(value)
                        else:
                            wdg.setValue(config.getint("Plot", label.text()))
                    elif isinstance(wdg, LineStyleCombo):
                        wdg.setCurrentText(config.get("Plot", label.text()))
                    elif isinstance(wdg, MarkerCombo):
                        wdg.setCurrentText(config.get("Plot", label.text()))
                    elif isinstance(wdg, QtWidgets.QComboBox):
                        wdg.setCurrentText(config.get("Plot", label.text()))
                    elif isinstance(wdg, QtWidgets.QCheckBox):
                        wdg.setChecked(config.getboolean("Plot", label.text()))
                    elif isinstance(wdg, ColorSelector):
                        wdg.setColor(config.get("Plot", label.text()))

        self.updatePlot()

    def updateStyle(self, textStyle):
        """Update rcParams with the new selected style and update plot"""
        rcdefaults()

        if textStyle != "default":
            rcParams.update(style.library[textStyle])

        self._setRcParams(rcParams)
        self.updatePlot()

    def updatePlot(self):
        """Update plot with the new configuration, global style or only a
        rcParams changed"""
        if self.customize.isChecked():
            textStyle = self._getRcParams()

        else:
            textStyle = self.style.currentText()

        self.plot(textStyle)

    def plot(self, textStyle):
        """Do the sample plot using the specified style confuguration, the
        style configuration can be a named matplotlib style or a dictionary
        with rc_params"""
        with style.context(textStyle):
            x = arange(-2, 8, .2)
            y = .1 * x ** 3 - x ** 2 + 3 * x + 2
            z = .01 * x ** 3 + x ** 2 - 3 * x - 2

            plot = PlotWidget(width=1, height=1)
            plot.plot(x, y, label="f(x)")
            plot.plot(x, z, label="g(x)")
            plot.ax.set_title("Function representation")
            plot.ax.set_xlabel("x")
            plot.ax.set_ylabel(r"$y=f\left(x\right)$")
            plot.ax.legend()
            self.layout().addWidget(plot, 9, 1, 1, 3)

    def _getRcParams(self):
        """Get rc parameters from widgets"""
        kw = {}
        for idx in range(self.tabRcParams.count()):
            tab = self.tabRcParams.widget(idx).widget()

            # The first children is always the layout used so we discard it
            for parent_wdg in tab.children()[1:]:
                label, wdg = parent_wdg.children()[1:]

                if isinstance(wdg, QtWidgets.QDoubleSpinBox):
                    # Correct value for axes.titley
                    if label.text() == "axes.titley" and wdg.value() == 1.:
                        value = None
                    else:
                        value = wdg.value()
                elif isinstance(wdg, QtWidgets.QSpinBox):
                    # Correct value to matplotlib expected format
                    if label.text() == "axes.formatter.limits":
                        value = f"-{wdg.value()}, {wdg.value()}"
                    else:
                        value = wdg.value()
                elif isinstance(wdg, (LineStyleCombo, MarkerCombo)):
                    value = wdg.currentValue()
                elif isinstance(wdg, QtWidgets.QComboBox):
                    value = wdg.currentText()
                elif isinstance(wdg, QtWidgets.QCheckBox):
                    value = wdg.isChecked()
                elif isinstance(wdg, ColorSelector):
                    value = wdg.color.name()
                else:
                    value = None

                kw[label.text()] = value
        return kw

    def _setRcParams(self, rcparams):
        """Populate widgets with the dict rcparams given as parameter"""
        for idx in range(self.tabRcParams.count()):
            tab = self.tabRcParams.widget(idx).widget()
            for parent_wdg in tab.children()[1:]:
                label, wdg = parent_wdg.children()[1:]
                wdg.blockSignals(True)
                value = rcparams.get(label.text(), None)

                if isinstance(wdg, QtWidgets.QDoubleSpinBox):
                    # Correct default value
                    if label.text() == "axes.titley" and value is None:
                        value = 1.
                    if label.text() == "legend.framealpha" and value is None:
                        value = 0.8
                    wdg.setValue(value)
                elif isinstance(wdg, QtWidgets.QSpinBox):
                    # Correct default value
                    if label.text() == "savefig.dpi" and value == "figure":
                        value = rcparams["figure.dpi"]
                    if label.text() == "axes.formatter.limits":
                        value = -value[0]
                    wdg.setValue(int(value))
                elif isinstance(wdg, LineStyleCombo):
                    wdg.setCurrentValue(value)
                elif isinstance(wdg, MarkerCombo):
                    wdg.setCurrentValue(value)
                elif isinstance(wdg, QtWidgets.QComboBox):
                    # Do text manipulation for joinstyle enum
                    if label.text() in (
                            "lines.dash_joinstyle", "lines.solid_joinstyle"):
                        value = value.split(".")[-1]
                    # Do text manipulation for capstyle enum
                    if label.text() in (
                            "lines.dash_capstyle", "lines.solid_capstyle"):
                        value = value.split(".")[-1]
                    # Do text manipulation for font.family
                    if label.text() == "font.family":
                        value = value[0]
                    wdg.setCurrentText(str(value))
                elif isinstance(wdg, QtWidgets.QCheckBox):
                    wdg.setChecked(value)
                elif isinstance(wdg, ColorSelector):

                    # Correct auto and inherit values
                    if label.text() == "axes.titlecolor" and value == "auto":
                        value = rcparams["text.color"]
                    if label.text() == "legend.facecolor" \
                            and value == "inherit":
                        value = rcparams["axes.facecolor"]
                    if label.text() == "legend.edgecolor" \
                            and value == "inherit":
                        value = rcparams["axes.edgecolor"]
                    if label.text() == "xtick.labelcolor" \
                            and value == "inherit":
                        value = rcparams["xtick.color"]
                    if label.text() == "ytick.labelcolor" \
                            and value == "inherit":
                        value = rcparams["ytick.color"]
                    if label.text() == "savefig.facecolor" and value == "auto":
                        value = rcparams["figure.facecolor"]
                    if label.text() == "savefig.edgecolor" and value == "auto":
                        value = rcparams["figure.edgecolor"]
                    if label.text() == "lines.markeredgecolor" \
                            and value == "auto":
                        value = "#000000"
                    if label.text() == "lines.markerfacecolor" \
                            and value == "auto":
                        value = "C2"

                    # Get the RGB code of color to avoid the several color
                    # definition in matplotlib and its styles
                    r, g, b = to_rgb(value)

                    r = int(float(r)*255)
                    g = int(float(g)*255)
                    b = int(float(b)*255)
                    value = f"#{r:x}{g:x}{b:x}"

                    wdg.setColor(value)

                wdg.blockSignals(False)

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("Plot"):
            config.add_section("Plot")
        config.set("Plot", "style", str(self.style.currentIndex()))
        config.set("Plot", "customize", str(self.customize.isChecked()))

        for k, val in self._getRcParams().items():
            config.set("Plot", k, str(val))

        return config


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
