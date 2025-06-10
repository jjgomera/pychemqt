#!/usr/bin/sh
# -*- coding: utf-8 -*-


#Script for dependences instalation in a debian based distribution,
#tested in a debian jessie and stretch


########################################################################
#Mandatory dependences
########################################################################
apt-get install python3 python3-numpy python3-scipy python3-matplotlib 
apt-get install python3-pyqt5 



########################################################################
# IAPWS
# External package to implement the IAPWS standard
# Preferred the github option to get the last version
########################################################################

#pip3 install iapws
pip3 install https://github.com/jjgomera/iapws/archive/master.zip


########################################################################
# pyqt5-qscintilla
# pagkage to show code
########################################################################

apt-get python3-pyqt5.qsci


########################################################################
#Freesteam
#package for calculating thermodynamic properties of water by IAPWS-IF97
########################################################################

# We need scons and swig3.0 to compile
apt-get install scons swig3.0 libgsl0-dev
# Probably you have yet but anyway you need install the development files
# of current python3 version
apt-get install python-dev

# Proceed to download and compile freesteam
wget -Nc https://sourceforge.net/projects/freesteam/files/freesteam/2.1/freesteam-2.1.tar.bz2
tar -xf freesteam-2.1.tar.bz2
cd freesteam-2.1
scons

# scons only work for python2.x, so the python binding created work only for the python2.x
# installed in your system
# Check your python3 version installed:
#    python3.4 in jessie
#    python3.5 in stretch
# For create the python3 binding:
swig3.0 -o python/freesteam_wrap.c -python python/freesteam.i
gcc -o python/freesteam_wrap.os -c -Wall -W -Wconversion -Wimplicit -fPIC -I. -I/usr/include/python3.5 python/freesteam_wrap.c
gcc -o python/_freesteam.so -shared python/freesteam_wrap.os -L. -L/usr/libs -lfreesteam

# Finally move the python3 binding to a forder in path, need root user
cp python/_freesteam.so /usr/local/lib/python3.5/dist-packages/
cp python/freesteam.py /usr/local/lib/python3.5/dist-packages/
cp libfreesteam.so /usr/lib
cp libfreesteam.so.1 /usr/lib


########################################################################
#CoolProp
#package for calculating thermodynamic properties with advanced equations
########################################################################

apt-get install pip3
pip3 install CoolProp


########################################################################
#refProp
#package for calculating thermodynamic properties with advanced equations
########################################################################

# This library is a python wrapper to work with a installed copy of REFPROP.
# In order to legally obtain a copy of the REFPROP software and fortran source
# codes, use the following web page: http://www.nist.gov/srd/nist23.cfm.
# Once purchased, the software will be delivered as a windows self extracting
# installation file (.exe). In order to access the files required by this python module,
# you can install REFPROP using this self extracting installation file or extract files.

# First at all instal the dependences
apt-get install dos2unix gfortran

# Download the python library
git clone https://github.com/BenThelen/python-refprop.git

# Compile the REFPROP source files
cd python-refprop
./rp2so
# The script need root privileges to run
# The script ask for the refprop files folder and the destination folder for compiled files,

# We suppose the python3 library (for python3.2) is going to work in other python3 versions
cp python3.2/multiRP.py /usr/local/lib/python3.5/dist-packages
cp python3.2/refprop.py /usr/local/lib/python3.5/dist-packages
cp python3.2/rptest.py /usr/local/lib/python3.5/dist-packages


########################################################################
# spreedsheet supoort

########################################################################
#lxml is a required dependendence of ezodf don't configured in setup
#so we need install too
pip3 install ezodf lxml
pip3 install openpyxl
pip3 install xlwt
pip3 install openbabel  # work too with python3-openbabel in debian bullseye


###########################################################################
# Other optional dependences

###########################################################################

apt install python3-icu
