.. image:: https://travis-ci.org/jjgomera/pychemqt.svg?branch=master
    :target: https://travis-ci.org/jjgomera/pychemqt
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/jjgomera/pychemqt/badge.svg?branch=master
    :target: https://coveralls.io/github/jjgomera/pychemqt?branch=master
    :alt: coveralls.io analysis

.. image:: https://codecov.io/gh/jjgomera/pychemqt/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/jjgomera/pychemqt
    :alt: codecov.io analysis

.. image:: https://readthedocs.org/projects/pychemqt/badge/?version=latest
    :target: http://pychemqt.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


Introduction
============
pychemqt is intended to be a free software tool for calculation and design of unit operations in chemical engineering. The goal is to develop free software equivalent to CHEMCAD or hysys. It is written in python using qt as a graphics library, so that it is cross-platform.


Dependencies
============

* `python3 <http://www.python.org/>`__, version 3.x required
* `pyqt binding <http://www.riverbankcomputing.co.uk/news>`__, support both PyQt5 and PyQt6
* `Numpy-scipy <http://scipy.org/Download>`__: python library for mathematical computation
* `matplotlib <http://matplotlib.sourceforge.net/>`__: python library for graphical representation of data
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

* UI with support for flow diagrams
* Database with 800 components
* Definition of custom compounds
* Stream definition with temperature, pressure and composition
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
* Tools
	* `Units converter <tools.UI_unitConverter.html>`__
	* Currency converter
	* `Periodic table of elements <tools.qtelemental.html>`__
	* `Psychometric chart <tools.UI_psychrometry.html>`__
	* `High quality calculation of properties using multi-parameter equations <tools.UITables.html>`__
	* Introspection support with a python shell (Linux only)

* Configurable: Units system, property correlation, EoS to use...
* Internationalization support: english, spanish.


TODO
====

* Make more equations of state available
* Improve the gui
* Add more equipment: complete heat exchanger, distillation columns, reactors...
* Clean code and debug bugs
* Add testing of libraries
* Improve documentation

For any suggestions, comments, bugs ... you can contact me at `email <jjgomera@gmail.com>`__.
