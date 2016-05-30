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

pip3 install ezodf
