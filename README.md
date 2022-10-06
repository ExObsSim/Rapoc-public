[![PyPI version](https://badge.fury.io/py/rapoc.svg)](https://badge.fury.io/py/rapoc)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/ExObsSim/Rapoc-public?color=gree&label=GitHub%20release)
[![Downloads](https://pepy.tech/badge/rapoc)](https://pepy.tech/project/rapoc)
[![Documentation Status](https://readthedocs.org/projects/rapoc-public/badge/?version=latest)](https://rapoc-public.readthedocs.io/en/latest/?badge=latest)
![GitHub](https://img.shields.io/github/license/ExObsSim/Rapoc-public)

# RAPOC: Rosseland And Planck Opacity Converter

The `RAPOC` code is written by Lorenzo V. Mugnai and Darius Modirrousta-Galian and is the product 
of a collaboration between Sapienza Università di Roma, Università degli Studi di Palermo and 
INAF - Osservatorio Astronomico di Palermo. It uses molecular absorption measurements 
(i.e. wavelength-dependent opacities) to calculate Rosseland and 
Planck mean opacities that are commonly used in atmospheric modelling.

`RAPOC` is designed to be simple, straightforward, and easily incorporated 
into other codes. It is completely written in Python and documented with docstrings. 
In addition, a Sphinx version of the documentation with a full user guide 
that includes examples is available in html format.

### Reports
`RAPOC` is under development, please report any issues or inaccuracies 
to the developers to support the implementation.

### Cite
If you use this code or its results, please cite `RAPOC: the Rosseland and Planck opacity converter` by Mugnai L. V. and Modirrousta-Galian D. (submitted).

## Installation
### Installing from Pypi
`RAPOC` can be installed from the Pypi repository with the following script::

    pip install rapoc

### Installing from git
`RAPOC` may also be cloned from the main git repository::

    git clone https://github.com/ExObsSim/Rapoc-public.git

The next step is to move into the `RAPOC` folder::

    cd /your_path/Rapoc

Then::

    pip install .

To check if one has the correct setup::

    python -c "import rapoc"


## Use
`RAPOC` is designed to be used on its own or in conjunction with other Python 
codes. Given an ExoMol file in the TauREx.h5 format, Rosseland and Planck mean opacities can be calculated. 
For example, in order to estimate the mean opacities at a temperature (T) of 1000 K with a pressure (P) of 
10,000 Pa in the wavelength range of  0.3-50 micron the following script is used,
    
    from rapoc import Rosseland, Planck

    r_model = Rosseland(input_data='exomol_file.TauREx.h5')
    opacity = r_model.estimate(P_input=10000 * u.Pa, T_input=1000 * u.K, band=(0.3 *u.um, 50*u.um))

    p_model = Planck(input_data='exomol_file.TauREx.h5')
    opacity = p_model.estimate(P_input=10000 * u.Pa, T_input=1000 * u.K, band=(0.3 *u.um, 50*u.um))

### Inputs
To run the code you need measured data. The supported file formats are:

- ExoMol opacities (downloadable [here](http://exomol.com/data/data-types/opacity)) with the `TauREx.h5` format.

## Documentation
The full documentation is available [here](https://rapoc-public.readthedocs.io/en/latest/) 

Alternatively, `RAPOC` accepts user-defined documentation by using `sphinx`. To install it run
    
    pip install sphinx sphinx_rtd_theme
    
From the `Rapoc/docs` folder running
    
    cd docs
    make html

This will create a html version of the documentation in `Rapoc/doc/build/html/index.html`.
 