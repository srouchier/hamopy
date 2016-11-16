======
hamopy
======

[![DOI](https://zenodo.org/badge/11513/srouchier/hamopy.svg)](https://zenodo.org/badge/latestdoi/11513/srouchier/hamopy)

Hamopy is a python package for the numerical simulation of one-dimensional heat, air and moisture (HAM) transfer in porous materials. Its principle is the finite-element resolution of the HAM conservation equations. The only requirement to run it is a fairly recent version of SciPy.

The original field of application is the hygrothermal modelling of building materials, although the code is not restricted to it.

Hamopy makes good use of its open-source nature, and gives users complete control over the simulation process. One can:

* add new materials and customise the equations defining their properties,

* account for water flow and storage in both liquid and vapor states,

* include time-dependent boundary conditions,

* work with fully coupled hygrothermal transfer, or with thermal transfer only (saves time),

* easily automate many simulations for sensitivity analyses, evolutionary algorithms and such.

The documentation is available here: http://srouchier.github.io/hamopy/

Any contribution into improving hamopy is welcome, as to make open-source HAM modelling available and understandable by all.
