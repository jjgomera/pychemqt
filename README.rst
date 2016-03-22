pychemqt intended as a free software tool for calculation and design of unit operations in chemical engineering. The goal is to obtain an equivalent free software to CHEMCAD or hysys. It is written in python using qt as graphics libraries, so is cross-platform.


Dependencies
--------------------

* `python3 <http://www.p7ython.org/>`__, version 3.x required
* `pyqt5 <http://www.riverbankcomputing.co.uk/news>`__, developed with version 5.3 
* `Numpy-scipy <http://scipy.org/Download>`__: python library for mathematical computation
* `matplotlib <http://matplotlib.sourceforge.net/>`__: python library for graphical representation of data
* `python-graph <http://code.google.com/p/python-graph/>`__:  python library for working with graphs

Optional applications, pychemqt work but some options will be disabled 

* `freesteam <http://freesteam.sourceforge.net/>`__: package for calculating thermodynamic properties of water by IAPWS-IF97
* `python-refprop <https://github.com/BenThelen/python-refprop>`__: package for calculating thermodynamic properties using `refprop <http://www.nist.gov/srd/nist23.cfm>`__ NIST application
* `coolprop <http://coolprop.org/>`__: package for calculating thermodynamic properties using multiparameter equation of state
* `oasa <http://bkchem.zirael.org/oasa_en.html>`__: used to show compound extended formula in database
* `ezodf <https://bitbucket.org/mozman/ezodf>`__: package to integration with OpenDocument spreadsheet (ods)
* `openpyxl <https://bitbucket.org/ericgazoni/openpyxl>`__: package to integration with Microsoft Excel 2007/2010 (xlsx)
* `xlwt <https://pypi.python.org/pypi/xlwt>`__: package to integration with Microsoft Excel 97/2000/XP/2003 (xls)


Features
--------------------

The development is slow, so the software in in pre-alpha status, with many bugs and with only a few features implemented:

* UI with support for flow diagram
* Databank with 800 components
* Let define hypothetical compound 
* Stream definition with temperature, pressure and composition
* Thermodinamic:
	* Redlich-Kwong (RK)
	* Soave-Redlich-Kwong (SRK)
	* Modificada Soave-Redlich-Kwong (MSRK)
	* Peng-Robinson (PR)
	* Peng-Robinson-Stryjek-Vera (PRSV)
	* Benedict-Webb-Rubin-Starling (BWRS)
	* Lee-Kesler
	* EoS multiparameter type Setzmann-Wagner for several pure fluids
	* GERG EoS for mix (Partial)
* Equipments:
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
	* Psychrometric chart
	* Units converter 
	* Currency converter
	* Introspection support with a python shell (Linux only)
	* High quality properties calculation using multiparameter equations 


* Configurable: Units system, property correlation, EoS to use...
* Support units systems.
* Internationalization support: english, spanish.



TODO
--------------------

* Enlarge equation of state available
* Improve gui
* Add more equipment: complete heat exchanger, distillation columns, reactors...
* Clean code and debug bugs
* Improve documentation

For any suggestions, comments, bug ... you can contact me at `email <jjgomera@gmail.com>`__.
