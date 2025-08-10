.. image:: https://dl.circleci.com/status-badge/img/gh/jjgomera/pychemqt/tree/master.svg?style=svg
    :target: https://dl.circleci.com/status-badge/redirect/gh/jjgomera/pychemqt/tree/master
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/jjgomera/pychemqt/badge.svg?branch=master
    :target: https://coveralls.io/github/jjgomera/pychemqt?branch=master
    :alt: coveralls.io analysis

.. image:: https://codecov.io/gh/jjgomera/pychemqt/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/jjgomera/pychemqt
    :alt: codecov.io analysis

.. image:: https://app.codacy.com/project/badge/Grade/457297f080904ae5aa2ae52a4c1e7f9d
    :target: https://www.codacy.com/gh/jjgomera/pychemqt/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jjgomera/pychemqt&amp;utm_campaign=Badge_Grade
    :alt: Code Quality

.. image:: https://readthedocs.org/projects/pychemqt/badge/?version=latest
    :target: http://pychemqt.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


Introduction
============
pychemqt is intended to be a free software tool for calculation and design of unit operations in chemical engineering. The goal is to develop free software equivalent to CHEMCAD or hysys. It is written in python using qt as a graphics library, so that it is cross-platform.


Dependencies
============

* `python3 <https://www.python.org/>`__, version 3.x required
* `pyqt binding <https://riverbankcomputing.com/software/pyqt>`__, support both PyQt5 and PyQt6
* `Numpy <https://numpy.org/>`__: python library for numerical computation
* `Scipy <https://scipy.org/>`__: python library for mathematical computation
* `matplotlib <https://matplotlib.org/>`__: python library for plotting
* `numdifftools <https://github.com/pbrod/numdifftools>`__: python library for numerical differentiation
* `iapws <https://github.com/jjgomera/iapws/>`__: python library for thermodynamic properties of water by IAPWS standards

Optional applications that are required for pychemqt to work with full functionality:

* `freesteam <http://freesteam.sourceforge.net/>`__: package for calculating thermodynamic properties of water by IAPWS-IF97
* `coolprop <http://coolprop.org/>`__: package for calculating thermodynamic properties using multiparameter equations of state, required 6.x version.
* `python-refprop <https://github.com/BenThelen/python-refprop>`__: package for calculating thermodynamic properties using `refprop <http://www.nist.gov/srd/nist23.cfm>`__ NIST application
* `openbabel <http://openbabel.org/wiki/Main_Page>`__: used to show compound extended formula in database
* `ezodf <https://bitbucket.org/mozman/ezodf>`__: package for integration with OpenDocument spreadsheet (ods)
* `openpyxl <https://bitbucket.org/ericgazoni/openpyxl>`__: package for integration with Microsoft Excel 2007/2010 (xlsx)
* `xlwt <https://pypi.python.org/pypi/xlwt>`__: package for integration with Microsoft Excel 97/2000/XP/2003 (xls)
* `reportlab <https://bitbucket.org/rptlab/reportlab>`__: package to export pdf reports
* `PyQt5-Qscintilla <https://riverbankcomputing.com/software/qscintilla/intro>`__: Custom code viewer and editor with syntax highlight


Features
========

The development is slow, so the software in in pre-alpha status, with many bugs and with only a few features implemented:

* `UI <ui.html>`__ with support for flow diagrams
* `Database with 1000 predefined components <tools.UI_databank.html>`__
* Definition of custom compounds
	* `Petroleum fraction pseudocomponent <pseudocomponent.html>`__
	* `Group contribution methods <lib.newComponent.html>`__
* `Stream definition <UI.UI_corriente.html>`__ with temperature, pressure and composition
* Thermodynamic EoS:
	* Redlich-Kwong (RK)
	* Soave-Redlich-Kwong (SRK)
	* Modified Soave-Redlich-Kwong (MSRK)
	* Peng-Robinson (PR)
	* Peng-Robinson-Stryjek-Vera (PRSV)
	* Benedict-Webb-Rubin-Starling (BWRS)
	* Lee-Kesler
	* EoS multi-parameter type Setzmann-Wagner for several pure fluids
	* GERG EoS for mixtures (Partial)
* Equipment:
	* Divider
	* Mixer
	* Valve
	* Pipe (Partial)
	* Compressor
	* Expander
	* Pump
	* Generic Heat Exchanger (without design)
	* Double Pipe Heat Exchanger (Partial)
	* Shell and Tube Heat Exchanger (Partial)
	* Fired Heater Heat Exchanger
	* Flash LV
	* Distillation column (simple method FUG)
	* Ciclon
	* Gravity Chamber
	* ElectricPrecipitator
	* Baghouse
	* Spreadsheet equipment (ods,xlsx)
* Tools:
	* `Units converter <tools.UI_unitConverter.html>`__
	* Currency converter
	* `Periodic table of elements <tools.qtelemental.html>`__
	* `Psychometric chart <tools.UI_psychrometry.html>`__
	* `High quality calculation of properties using multi-parameter equations <tools.UI_Tables.html>`__
	* Introspection support with a python shell (Linux only)
	* Tiny plot
		* `Moody chart <plots.moody.html>`__
		* `Drag sphere chart <plots.drag.html>`__
		* `Standing-Katz chart <plots.standing.html>`__

* `Calculation method fully configurable <tools.wizard.html>`__: Units system, property correlation, EoS to use...
* `Complete set of preferences <tools.UI_Preferences.html>`__: For adjust the gui to the user preferences
* Basic `cost estimation <tools.costIndex.html>`__ utility
* Internationalization support: english, spanish.
* `Tools to show references <tools.doi.html>`__ for academic purposes


TODO
====

* Make more equations of state available
* Improve the gui
* Add more equipment: complete heat exchanger, distillation columns, reactors...
* Clean code and debug bugs
* Add testing of libraries
* Improve documentation

For any suggestions, comments, bugs ... you can contact me at `email <jjgomera@gmail.com>`__.
