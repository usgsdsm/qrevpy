# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='qrevpy',
    version='0.1.0',
    description='QRev port to python',
    long_description=readme,
    author='David S. Mueller',
    author_email='dmueller@usgs.gov',
    url='https://hydroacoustics.usgs.gov/movingboat/QRev.shtml',
    license=license,
    packages=find_packages(exclude=('MiscLibs',))
)