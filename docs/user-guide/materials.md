# Creating materials

As shown on the [input definition](inputs.md) page, creating the `mesh` object for the simulation requires including a material object in the arguments of the `Mesh` instantiation.

	from hamopy.classes import Mesh
	from hamopy.materials.standard import wood_fibre

	mesh = Mesh(materials = [wood_fibre],
		    sizes     = [0.08],
		    nbr_elem  = [16] )

This means that the `hamopy\materials\standard.py` file contains the definition of an object called `wood_fibre`, which is of the `Material` class and includes all material data.

Material definition in hamopy is a bit tricky, as it was designed to give the user maximum freedom for defining fully customisable transport properties.

This page tells you how to

* define a new material,
* use its properties once they have been registered,
* modify a material in an existing mesh.

## Material definition

### Instantiation

For purposes of reusability, the definition of a new material should preferably be done in a separate file which you can then import into your script (like `hamopy.materials.standard` above)

Once there, you may start by instantiating a new `Material` object:

	from hamopy.classes import Material

	concrete = Material('concrete', rho = 2400., cp = 880.)

The instantiation may take up to three arguments: the material name, its dry density (kg/m3) and specific heat (J/(kg.K)).

Once created, the `Material` object contains the following set of methods to register all material properties.

### set_density()

If not provided at the instantiation, you can specify the dry density of the material like so:

	concrete.set_density(2430.)


### set_capacity()

In the current state of hamopy, the heat capacity (J/(kg.K)) of the dry material can be defined as a function of the temperature.

	concrete.set_capacity(cp_0 = 2430., cp_t = 27.)

When called, the value of the heat capacity will then be `cp = cp_0 + T(°C) * cp_t`

The default value for `cp_t` is 0.

### set_conduc()

The thermal conductivity (W/(m.K)) may be defined as a function of the moisture content and temperature.

	concrete.set_conduc(lambda_0 = 1.75, lambda_m = 4.5, lambda_t = 1e-4) 

When called, the value of the conductivity will then be `lambda = lambda_0 + w(kg/m3)/1000 * lambda_m + T(°C) * lambda_t`

The default value for `lambda_m` and `lambda_t` is 0.

### set_isotherm()

There are currently two methods for defining the sorption isotherm:

Either a 3rd degree polynomial interpolation, fitted on a list of measurement points:

	wood_fibre.set_isotherm('polynomial', **{"HR" : [0, 0.25, 0.5, 0.75],
			                         "W"  : [0, 6.2, 12.4, 20.9] })

where the list given in the `'W'` key are the values of the moisture content measured at the relative humidities `'HR'`.

The second method for defining the sorption isotherm is the van Genuchten mono- or multimodal law:

	lightweight.set_isotherm('vangenuchten', **{"w_sat" : 871,
		                                    "l"     : [0.41, 0.59],
		                                    "alpha" : [6.12e-7, 1.22e-6],
		                                    "m"     : [0.5981, 0.5816] })

### set_perm_vapor()

There are currently two methods for defining the water vapor permeability:

Either by interpolation between measurement points:

	concrete.set_perm_vapor('interp', **{"HR" : [0.25, 0.75],
		                             "dp" : [4.2e-12, 7.8e-12] } )

Ot with the Schirmer law:

	lightweight.set_perm_vapor('schirmer', **{"mu" : 5.6,
		                                  "p"  : 0.2 })

### set_perm_liquid()

There are currently two methods for defining the liquid permeability:

Either with an exponential law:

	lightweight.set_perm_liquid('exp', **{"a" : [-46.245, 294.506, -1439, 3249, -3370, 1305] } )


Or with the Durner multi-modal law:

	concrete.set_perm_liquid('durner', **{"K_sat" : 2.2182e-13,
		                              "tau"   : -4.6974,
		                              "l"     : [0.5062, 0.4938],
		                              "alpha" : [5.5383e-7, 2.2493e-8],
		                              "m"     : [0.6148, 0.1913] } )

Note that if this method is not used in the definition of a material, this permeability will be set to `0` and liquid transfer will not be considered in the calculation.

### set_perm_air()

The air permeability should be given in (m^2^). Only a constant value is expected:

	lightweight.set_perm_air(1.08e-10)

Calling this method is optional: the default air permeability is `0`.

## Calling properties

Once defined, all properties may be called as fonctions of the temperature `T` and capillary pressure `p_c`. The following methods are bound within the `Material` class for this purpose

* `material.rho(T)`: density
* `material.cp(T)`: heat capacity
* `material.conduc(p_c, T)`: heat conductivity (`p_c` and `T` are optional arguments)
* `material.w(p_c)`: moisture content (kg/m3)
* `material.c_v(p_c)`: moisture capacity (derivative of the moisture content)
* `material.delta_p(p_c, T)`: vapour permeability (`T` is an optional argument)
* `material.k_l(p_c, T)`: liquid permeability (`T` is an optional argument)
* `material.k_air()`: air permeability

**Important**: hamopy is written with the capillary pressure as the driving potential for moisture transfer. All variables depending on the humidity are therefore methods expecting `p_c` as input argument, rather than the value of relative humidity or vapor pressure. Some functions are however available in the [library](library.md) to easily switch from one another.

## Switch materials

Once a material has been integrated into a `Mesh` object, it is still possible to change some of its properties and tell the mesh of the modifications:

	from hamopy.materials.standard import concrete, wood_fibre
	from copy import deepcopy

	mesh = Mesh(**{"materials"    : [concrete, wood_fibre],
		       "sizes"        : [0.1, 0.08],
		       "nbr_elements" : [16, 12] })

	wood_fibre_2 = deepcopy(wood_fibre)
	wood_fibre_2.set_density(170.)

	mesh.replace_materials([concrete, wood_fibre_2])

This functionality is particularly interesting when running a large number of simulations with different material properties, like in case of a sensitivity analysis. The instantiation of the `Mesh` and `Material` objects is not repeated, which saves some time.

[Documentation main page](../index.md)
