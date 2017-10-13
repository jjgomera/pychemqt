The procedure to install dependences in windows isn't so easy as in POSIX. I'll try to explain here the procedure.

Current versions:
    * python v3.6
    * Link checked at 12/03/17
    * Procedure checked in windows 7 and windows 10


Using Anaconda
--------------
If the disk space isn't problem, we can install a python distribution with full of packages to boost it, like anaconda:

    * x86: https://repo.continuum.io/archive/Anaconda3-4.3.1-Windows-x86.exe
    * amd64: https://repo.continuum.io/archive/Anaconda3-4.3.1-Windows-x86_64.exe

We need only install iapws run the command from a power shell:

    * pip install https://github.com/jjgomera/iapws/archive/master.zip


Install only python and neccessary library
------------------------------------------

First we need install python3, from his `page <https://www.python.org/downloads/release/python-360/>`__, it depends of your computer architecture, x86, x86-64:

    * x86: https://www.python.org/ftp/python/3.6.0/python-3.6.0.exe 
    * amd64: https://www.python.org/ftp/python/3.6.0/python-3.6.0-amd64.exe

Now using the pip command shipped with python we can install dependences easy. Open a cmd console and run:

    * py -m pip install PyQt5
    * py -m pip install matplotlib

Scipy is a problematic dependence in windows, if we try to install like matplotlib we have error at compilation. We need find a precompiled versions fine for our system, same architecture and same python version. Surfing in the web we can find this `page <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`__ with full of python modules precompiled. We install numpy too for this repositories. So first we download the appropiate wheels:

    * numpy-x86: numpy‑1.12.0+mkl‑cp36‑cp36m‑win32.whl
    * numpy-amd64: numpy‑1.12.0+mkl‑cp36‑cp36m‑win_amd64.whl
    * scipy-x86: scipy‑0.18.1‑cp36‑cp36m‑win32.whl
    * scipy-amd64: scipy‑0.18.1‑cp36‑cp36m‑win_amd64.whl

To install both we are going to use the cmd terminal and using the command:

    * py -m pip install C:\folder\file_numpy.whl
    * py -m pip install C:\folder\file_scipy.whl

Remerber it can use TAB to autocomplete command. Respect the order because scipy depend numpy, first install numpy wheel.

Optional dependences:
The spreedsheet support is easy to get, running the command:

    * py -m pip install ezodf lxml
    * py -m pip install openpyxl
    * py -m pip install xlwt
    * py -m pip install openbabel


Get pychemqt code
-----------------

To get pychemqt download the zip file from https://github.com/jjgomera/pychemqt/archive/master.zip and unzip whatever you want. To run doble click over the main script pychemqt.py.
