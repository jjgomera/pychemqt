#!/usr/bin/sh
# -*- coding: utf-8 -*-


#Script for dependences instalation in a debian based distribution, tested in a debian stable


########################################################################
#Mandatory dependences
########################################################################
apt-get install python3 python3-pyqt5 python3-numpy python3-scipy python3-matplotlib 


########################################################################
#Freesteam
#package for calculating thermodynamic properties of water by IAPWS-IF97
########################################################################

# We need scons and swig3.0 to compile
apt-get install scons swig3.0 python3.4-dev

# Proceed to download and compile freesteam
wget -Nc https://sourceforge.net/projects/freesteam/files/freesteam/2.1/freesteam-2.1.tar.bz2
tar -xf freesteam-2.1.tar.bz2
cd freesteam-2.1
scons

# scons only work for python2.x, so the python binding created work only for the python2.x
# installed in your system
# For create the python3 binding:
swig3.0 -o python/freesteam_wrap.c -python python/freesteam.i
gcc -o python/freesteam_wrap.os -c -Wall -W -Wconversion -Wimplicit -fPIC -I. -I/usr/include/python3.4 python/freesteam_wrap.c
gcc -o python/_freesteam.so -shared python/freesteam_wrap.os -L. -L/usr/libs -lfreesteam

# Finally move the python3 binding to a forder in path
cp python/_freesteam.so /usr/local/lib/python3.4/dist-packages/
cp python/freesteam.py /usr/local/lib/python3.4/dist-packages/


########################################################################
#CoolProp
#package for calculating thermodynamic properties with advanced equations
########################################################################

apt-get install pip3
pip3 install CoolProp


########################################################################
# IAPWS
# External package to implement the IAPWS standard
########################################################################

pip3 install iapws


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
cp python3.2/multiRP.py /usr/local/lib/python3.4/dist-packages
cp python3.2/refprop.py /usr/local/lib/python3.4/dist-packages
cp python3.2/rptest.py /usr/local/lib/python3.4/dist-packages





pip3 install ezodf

#OASA  used to show compound extended formula in database
#Luego de instalado bkchem
sudo python /usr/lib/bkchem/bkchem/oasa/setup.py install
