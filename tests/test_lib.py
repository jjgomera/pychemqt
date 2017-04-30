#!/usr/bin/python
# -*- coding: utf-8 -*-

from doctest import DocTestSuite

from lib import unidades
from lib import adimensional

TestUnidades = DocTestSuite(unidades)
TestAdimensional = DocTestSuite(adimensional)
