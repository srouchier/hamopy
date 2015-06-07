# Input definition

The hamopy algorithm is called like so:

	from hamopy.algorithm import calcul
	results = calcul(mesh, clim, init, time)

This means that all inputs of the simulation are contained within 4 objects that serve as arguments for the `calcul` method:

* `mesh` contains all information regarding the material definition and spatial discretisation of the domain
* `clim` stores boundary conditions
* `init` defines the initialisation of the domain
* `time` specifies the time of simulation, time step size, etc.

Below is the syntax for the definition of these objects.

## mesh: material definition and discretisation

The `mesh` argument is created by instantiating a new object of the `hamopy.classes.Mesh` class.

	from hamopy.classes import Mesh
	from hamopy.materials.standard import wood_fibre

	mesh = Mesh(materials    = [wood_fibre],
            	sizes        = [0.08],
            	nbr_elements = [16] )

The example above imports the wood fibre board from the `hamopy.materials.standard` material library. Then, it creates a wall of thickness 0.08 m, discretised with 16 finite elements.

Here is an example of a multi-layered wall:

	from hamopy.classes import Mesh
	from hamopy.materials.standard import concrete, wood_fibre

	mesh = Mesh(materials    = [concrete, wood_fibre],
	            sizes        = [0.10, 0.08],
	            nbr_elements = [24, 16] )

From left to right, this example creates a wall with a 10 cm concrete layer divided into 24 finite elements, and an 8 cm wood fibre layer with 16 finite elements.

The instantiation of the `Mesh` class always takes 3 arguments: `materials`, `sizes` and `nbr_elements`. These arguments are *always lists*, even when only one layer is involved as in the first example.

Materials must be imported before the mesh is created. They are objects of the `Material` class: the [materials](materials.md) page of this site shows how to define one.

## clim: boundary conditions

The `clim` argument sent to the simulation is a list of two objects of the `Boundary` class.

	from hamopy.classes import Boundary
	
	clim_01 = Boundary('Dirichlet',**{"T"  : 293.15,
	                                  "HR" : 0.7 })
	
	clim_02 = Boundary('Fourier',**{"T"   : 278.15,
	                                "HR"  : 0.9,
	                                "h_t" : 5 })
	
	clim = [clim_01, clim_02]

Two boundaries are created separately and stored in the `clim` list. The first one is a constant type one (Dirichlet) condition of 293.15 K and 0.7 relative humidity. The second one is a type three (Fourier) conditions with a surface transfer coefficient of 5 [W/(m2K)].

The instantiation of a boundary requires two arguments: a string to specify the boundary type, and a dictionary to specify numerical data.

The example above is very simple, but things can get a little more complicated: boundary conditions can be variable in time and more parameters can be specified (air pressure, surface transfer coefficients...) The [boundary conditions](boundary.md) page shows how to do all that.

## init: initialisation

The `init` object is simply a dictionary containing the initial value of the temperature under the `'T'` key, and of the relative humidity under the `'HR'` key:

	init = {'T'  : 293.15,
	        'HR' : 0.95}

Alternatively, the relative humidity can be replaced by a value of vapor pressure with the `'PV'` key or of capillary pressure with the `'PC'` key.

It is possible to define non-uniform initial conditions, by adding an `'x'` key to the dictionary:

	init = {'x'  : [0, 0.05, 0.10],
	        'T'  : [273, 287, 291],
	        'HR' : [0.7, 0.4, 0.6] }

In this example, the initial conditions will be a linear interpolation of the values given at the coordinates 0, 0.05, and 0.10 m.

## time: temporal discretisation

The object storing information on the discretisation in time is an occurence of the `Time` class. Its instantiation resembles that of the `Boundary` object:

	from hamopy.classes import Time
	time = Time('constant',**{"delta_t" : 600,
	                          "t_max"   : 7200 } )

The definition of the `time` requires two input arguments: a string indicating whether the time step size should be constant or variable, and a dictionary containing numerical data. In the example above, it is declared that the time step size `delta_t` is 600 seconds, and that the simulation time `t_max` is 7200 seconds.

It is possible to ask for the time step size to automatically adapt to eventual convergence difficulties. This is done by setting the method string to `'variable'`, and optionnaly adding some keys to the dictionary:

* `'iter_max'`: maximum number of iterations in a time step (default is 12)
* `'delta_min'`: minimum time step size (default is 1e-3)
* `'delta_max'`: maximum time step size (default is 900)

Once all 4 objects are defined, you can start the [simulation](simulation.md).

[Documentation main page](../index.md)
