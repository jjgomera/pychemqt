from setuptools import setup, find_packages

import io


__version__ = "0.1.0"

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

    # packages=find_packages(exclude=["iapws"]),
    package_dir={'': '../'},
    packages=["pychemqt." + pkg for pkg in find_packages(exclude=["iapws"])],
    package_data={'': ['../README.rst', '../LICENSE'],
                  'images': ["*.png", "*/*.png", "*.jpg", "*/*.jpg",
                             "*/*.svg", "*/*.gif"],
                  'docs': ["*.rst"],
                  'dat': ["*"],
                  'i18n': ["*"],
                  'Samples': ["*.pcq"],
                  },
    exclude_package_data={'docs': ['*.mEoS.*.rst', "*_ref.rst"]},

    install_requires=['scipy>=0.14',
                      'numpy>=1.8',
                      'matplotlib>=1.4',
                      'iapws>=1.3'],
    extras_require={
        'CoolProp':  ["CoolProp>=6.1.0"],
        'openbabel':  ["openbabel>=2.4.1"],
        'spreadsheet':  ["openpyxl>=2.3.0", "xlwt>=1.2.0", "ezodf>=0.3.2"],
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
