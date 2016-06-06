pychemqt is intended to be a free software tool for calculation and design of unit operations in chemical engineering. The goal is to develop free software equivalent to CHEMCAD or hysys. It is written in python using qt as a graphics library, so that it is cross-platform.


Dependencies
--------------------

* `python3 <http://www.python.org/>`__, version 3.x required
* `pyqt5 <http://www.riverbankcomputing.co.uk/news>`__, developed with version 5.3 
* `Numpy-scipy <http://scipy.org/Download>`__: python library for mathematical computation
* `matplotlib <http://matplotlib.sourceforge.net/>`__: python library for graphical representation of data

Optional applications that are required for pychemqt to work with full functionality:

* `freesteam <http://freesteam.sourceforge.net/>`__: package for calculating thermodynamic properties of water by IAPWS-IF97
* `coolprop <http://coolprop.org/>`__: package for calculating thermodynamic properties using multiparameter equations of state, required 6.x version.
* `python-refprop <https://github.com/BenThelen/python-refprop>`__: package for calculating thermodynamic properties using `refprop <http://www.nist.gov/srd/nist23.cfm>`__ NIST application
* `oasa <http://bkchem.zirael.org/oasa_en.html>`__: used to show compound extended formula in database
* `ezodf <https://bitbucket.org/mozman/ezodf>`__: package for integration with OpenDocument spreadsheet (ods)
* `openpyxl <https://bitbucket.org/ericgazoni/openpyxl>`__: package for integration with Microsoft Excel 2007/2010 (xlsx)
* `xlwt <https://pypi.python.org/pypi/xlwt>`__: package for integration with Microsoft Excel 97/2000/XP/2003 (xls)
* `reportlab <https://bitbucket.org/rptlab/reportlab>`__: package to export pdf reports
* `PyQt5-Qscintilla <https://riverbankcomputing.com/software/qscintilla/intro>`__: Custom code viewer and editor with syntax highlight


Features
--------------------

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
	* Periodic table of elements
	* Steam Tables
	* Psychometric chart
	* Units converter 
	* Currency converter
	* Introspection support with a python shell (Linux only)
	* High quality calculation of properties using multi-parameter equations


* Configurable: Units system, property correlation, EoS to use...
* Internationalization support: english, spanish.



TODO
--------------------

* Make more equations of state available
* Improve the gui
* Add more equipment: complete heat exchanger, distillation columns, reactors...
* Clean code and debug bugs
* Improve documentation

For any suggestions, comments, bugs ... you can contact me at `email <jjgomera@gmail.com>`__.
