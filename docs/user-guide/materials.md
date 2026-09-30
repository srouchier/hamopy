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

## STEICOinternal

The standard library includes the STEICOinternal wood fibre internal insulation
board:

```python
from hamopy.materials.standard import steico_internal
```

Its properties come from the manufacturer's
[technical data sheet (PDF)](https://www.steico.com/fileadmin/user_upload/importer/downloads/40281de298de2e31019976a323500188/STEICO_technical-data-sheet_internal_GBR_en_i.pdf),
GBR/IRL version 3, dated 2026-09-21, page 1 (accessed 2026-09-29).
The [current datasheet](https://www.steico.com/tds_steicointernal_gbr-irl_en)
is also available from STEICO.

| Property | Value |
| --- | --- |
| Density | Approximately 160 kg/m³ |
| Specific heat capacity | 2100 J/(kg·K) |
| Declared thermal conductivity | 0.038 W/(m·K) |
| Water vapour diffusion resistance factor µ | 5 |

Conductivity and heat capacity are modelled as constants. Vapour permeability
uses `interp_mu` with constant µ = 5 across relative humidity; Hamopy applies
its usual temperature conversion from µ to permeability.

**This is a datasheet-based entry for heat-only use, not a complete coupled HAM
material.** The datasheet does not supply a sorption isotherm, liquid permeability
curve, or moisture/temperature conductivity coefficients. No substitute curves
are assumed. Liquid and air transport retain Hamopy's disabled defaults; these
are not measured zero permeabilities. The reported airflow resistivity is only
a lower bound and is not registered as an exact air permeability.

For example, a 60 mm board can be used with `hamopy.algorithm.calcul_thermo`:

```python
from hamopy.classes import Mesh

mesh = Mesh(materials=[steico_internal], sizes=[0.06], nbr_elements=[12])
```

Before using this entry in a coupled moisture simulation, make a `deepcopy` and
define a measured sorption isotherm with `set_isotherm`, together with suitable
liquid transport data using `set_perm_liquid`. Calling moisture storage methods
on the supplied entry fails because no isotherm is defined.

## STEICO materials from WUFI XML

The `hamopy.materials.steico` catalog builds 19 coupled heat/moisture materials
from STEICO's local XML files. First [download and verify the manufacturer
bundle](../../data/wufi/steico/README.md). With the files in the source checkout's
`data/wufi/steico/xml` directory, use a named material directly:

```python
from hamopy.classes import Mesh
from hamopy.materials.steico import steico_therm

mesh = Mesh(materials=[steico_therm], sizes=[0.06], nbr_elements=[12])
```

Load the whole library or select a liquid-transport mode:

```python
from hamopy.materials.steico import load_material, load_materials

materials = load_materials()  # dict of 19 independent Material objects
board = materials['steico_protect_h_dry']
absorption_board = load_material('steico_therm', liquid_transport='suction')

# For an installed package or XML files stored elsewhere:
materials = load_materials('/path/to/steico/xml')
```

Alternatively set `HAMOPY_STEICO_XML_DIR` to that directory before using named
imports. `MATERIAL_FILES` lists all identifiers and their source filenames.
Every load/import creates a fresh material. Imports perform no downloads; missing
files produce an error with setup instructions. The manufacturer XMLs remain
local and are excluded from version control.

The importer preserves every sorption point using monotone cubic interpolation
(PCHIP). `c_w` is the derivative of that interpolant with respect to capillary
pressure. Density, dry heat capacity, constant vapour resistance, and
conductivity coefficients come from the XML. The conductivity conversion is
`lambda = lambda10 * (1 + b*w/rho) + temperature_slope*(T - 283.15)`, where `b`
is the XML supplement in % per mass-percent and `T` is kelvin. This uses Hamopy's
additive thermal model and reference-temperature moisture properties.

The zero entries in the generated liquid functions are placeholders. Following
the [WUFI help, pages 31–34](https://wufi.de/download/WUFI41_help_e.pdf), the
importer reconstructs positive coefficients at the reference moisture `w80` and
free saturation `wf` using `Dws(w) = 3.8*(A/wf)**2 * 1000**(w/wf - 1)`.
It interpolates logarithmically between those anchors, uses zero transport below
`w80`, and holds the final coefficient above `wf`. Redistribution shares the
`w80` coefficient and divides the `wf` coefficient by the XML drying factor.
Conversion to Hamopy's pressure-based permeability is `k_l = D_w * c_w`.

**The default liquid mode is `redistribution`.** Use `suction` when modelling
absorption. Hamopy uses the selected curve throughout a simulation and does not
automatically switch it with rain. These objects support Hamopy's coupled solver;
they do not reproduce every WUFI solver convention. In particular, PCHIP sorption
interpolation is a Hamopy modelling choice. Air transport retains its disabled
default because these XMLs provide no air-permeability measurement.

`material.source` retains the original name, source path and SHA-256, database
version/date, scalar values, and the complete XML `Info` text describing the
source models and their limitations. These source curves are not all measured
data. **STEICOinternal is absent from this bundle** and retains its separate
datasheet-based, heat-only entry above.

To load another file with the same supported XML profile:

```python
from hamopy.materials.wufi import load_wufi_material

material = load_wufi_material('/path/to/material.xml')
```

The reader accepts SI version-2 records with tabulated sorption, constant vapour
resistance and dry heat capacity, and generation-mode-0 liquid and conductivity
functions. Unsupported models, units, inconsistent tables, and missing properties
raise errors rather than producing an incomplete simulation material.

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

The following methods can define the sorption isotherm:

Either a 3rd degree polynomial interpolation, fitted on a list of measurement points:

	wood_fibre.set_isotherm('polynomial', **{"HR" : [0, 0.25, 0.5, 0.75],
			                         "W"  : [0, 6.2, 12.4, 20.9] })

where the list given in the `'W'` key are the values of the moisture content measured at the relative humidities `'HR'`.

The second method for defining the sorption isotherm is the van Genuchten mono- or multimodal law:

	lightweight.set_isotherm('vangenuchten', **{"w_sat" : 871,
		                                    "l"     : [0.41, 0.59],
		                                    "alpha" : [6.12e-7, 1.22e-6],
		                                    "m"     : [0.5981, 0.5816] })

For a full moisture-storage table, use shape-preserving monotone interpolation:

```python
material.set_isotherm('tabulated', HR=[0, 0.5, 0.8, 1], W=[0, 5, 10, 100])
```

`HR` must increase strictly from 0 to 1, and `W` must be finite, nondecreasing,
and start at zero. The final value is the saturation moisture content. This
method reproduces the supplied points without overshoot and computes moisture
capacity from the interpolant's derivative. Positive capillary pressures are
clamped to saturation, with zero capacity beyond that endpoint.

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
