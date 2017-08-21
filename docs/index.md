# Welcome to the hamopy documentation

Hamopy is a python package for the numerical simulation of one-dimensional heat, air and moisture (HAM) transfer in porous materials. Its principle is the finite-element resolution of the HAM conservation equations. The only requirement to run it is a fairly recent version of SciPy.

The original field of application is the hygrothermal modelling of building materials, although the code is not restricted to it.

Hamopy makes good use of its open-source nature, and gives users complete control over the simulation process. One can:

* add new materials and customise the equations defining their properties,
* account for water flow and storage in both liquid and vapor states,
* include time-dependent boundary conditions,
* work with fully coupled hygrothermal transfer, or with thermal transfer only (saves time),
* easily automate many simulations for sensitivity analyses, evolutionary algorithms and such.

Any contribution into improving hamopy is welcome, as to make open-source HAM modelling available and understandable by all.

## Installation

You can download and install hamopy by cloning the [GitHub repository](https://github.com/srouchier/hamopy), or download and unpack the ZIP file.

Then include the path to the download directory to your PATH or PYTHONPATH environment variable.

## How things work

Once `hamopy` is detected by your Python installation, this is basically how a simulation is run:

	from hamopy.algorithm import calcul
	results = calcul(mesh, clim, init, time)

The first line imports the main algorithm of hamopy, the second line runs the simulation under specified conditions (the `mesh`, `clim`, `init` and `time` objects) and stores `results` as a python dictionary.

Of course, some questions remain unanswered, which is what this page is for.

## User guide

* [Inputs](user-guide/inputs.md): how to define the conditions of simulation
	* [Materials](user-guide/materials.md): how to create a new material
	* [Boundary conditions](user-guide/boundary.md): how to set up boundary conditions
* [Simulation](user-guide/simulation.md): some options to customise and monitor the simulation
	* [Library](user-guide/library.md): small library of useful functions within hamopy
* [Post-processing](user-guide/post-processing.md): how to visualise results

## Examples

* [Hamstad BM 3](examples/Hamstad_BM3.md): 3rd benchmark exercise of the Hamstad package
* [Hamstad BM 5](examples/Hamstad_BM5.md): 5th benchmark exercise of the Hamstad package
