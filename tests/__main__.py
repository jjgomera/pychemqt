#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Main script for unit testing

from unittest import TextTestRunner, TestSuite

from tests.test_meos import TestMEOS

suite = TestSuite()
suite.addTest(TestMEOS)

runner = TextTestRunner()
results = runner.run(suite)
