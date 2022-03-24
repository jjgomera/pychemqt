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
# Tools to create a python shell with pychemqt libraries imported
# For now only work in linux with xterm as terminal
###############################################################################

import atexit

from PyQt5 import QtCore, QtWidgets

from tools.firstrun import which


class XTerm(QtCore.QProcess):
    """Gui container for terminal widget"""
    def __init__(self, config, parent=None):
        super(XTerm, self).__init__(parent)
        self.config = config
        atexit.register(self.kill)
        self.show_term()

    def sizeHint(self):
        size = QtCore.QSize(400, 300)
        return size.expandedTo(QtWidgets.QApplication.globalStrut())

    def show_term(self):
        term = self.config.get("Applications", 'Shell')

        args = [
            "-bg", self.config.get("Applications", "backgroundColor"),
            "-fg", self.config.get("Applications", "foregroundColor"),
            # blink cursor
            "-bc",
            # title
            "-T", QtWidgets.QApplication.translate(
                "pychemqt", "pychemqt python console")]

        if self.config.getboolean("Applications", "maximized"):
            args.append("-maximized")

        if self.config.getboolean("Applications", 'ipython') and \
                which("ipython"):
            args.append("ipython3")
        else:
            args.append("python3")

        self.start(term, args)
        if self.error() == QtCore.QProcess.FailedToStart:
            print("xterm not installed")


if __name__ == "__main__":
    import sys
    from configparser import ConfigParser
    import os
    app = QtWidgets.QApplication(sys.argv)
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    preferences = ConfigParser()
    preferences.read(conf_dir+"pychemqtrc")

    terminal = XTerm(preferences)
    app.exec_()
