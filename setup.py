#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


import io

from setuptools import setup, find_packages


with open("VERSION") as version_file:
    __version__ = version_file.read().strip()

with io.open('README.rst', encoding="utf8") as file:
    long_description = file.read()

setup(
    name='pychemqt',
    author='Juan José Gómez Romera',
    author_email='jjgomera@gmail.com',
    url='https://github.com/jjgomera/pychemqt',
    description='pychemqt is intended to be a free software tool for '
                'calculation and design of unit operations in chemical '
                'engineering.',
    long_description=long_description,
    license="gpl v3",
    version=__version__,

    packages=find_packages(exclude=["iapws"]),
    package_data={'': ['../README.rst', '../LICENSE'],
                  'images': ["*.png", "*/*.png", "*.jpg", "*/*.jpg",
                             "*/*.svg", "*/*.gif"],
                  'docs': ["*.rst"],
                  'dat': ["*"],
                  'i18n': ["*"],
                  'Samples': ["*.pcq"],
                  },
    exclude_package_data={'docs': ['*.mEoS.*.rst', "*_ref.rst"]},

    install_requires=['scipy>=1.14',
                      'numpy>=1.26',
                      'matplotlib>=3.8',
                      'iapws>=1.5',
                      'numdifftools'],
    extras_require={
        'CoolProp': ["CoolProp>=6.1.0"],
        'openbabel': ["openbabel>=2.4.1"],
        'spreadsheet': ["openpyxl>=2.3.0", "xlwt>=1.2.0", "ezodf>=0.3.2"],
        'icu': ["PyICU>=2.1"],
        'reportlab': ["reportlab>=3.5.8"]},

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: X11 Applications :: Qt",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics"]
)
