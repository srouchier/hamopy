#Boundary conditions

## Constant conditions

The `clim` argument sent to the simulation is a list of two objects of the `Boundary` class. The instantiation of a boundary requires two arguments: a string to specify the boundary type, and a dictionary to specify numerical data. The [input definition](inputs.md) page shows a simple example of how boundary conditions can be defined:

	from hamopy.classes import Boundary

	clim_01 = Boundary('Dirichlet',**{"T"  : 293.15,
		                          "HR" : 0.7 })

	clim_02 = Boundary('Fourier',**{"T"   : 278.15,
		                        "HR"  : 0.9,
		                        "h_t" : 5 })

	clim = [clim_01, clim_02]

Two boundaries are created separately and stored in the `clim` list. The first one is a constant type one (Dirichlet) condition of 293.15 K and 0.7 relative humidity. The second one is a type three (Fourier) conditions with a surface transfer coefficient of 5 [W/(m2K)].

## Variable conditions

The example above is very simple, but things can get a little more complicated: boundary conditions can be variable in time and more parameters can be specified (air pressure, surface transfer coefficients...)

Let's say one of your boundaries has a variable temperature and air pressure, which are stored in a .txt file **with headers**. In this case, this is how you instantiate a new boundary with the second column as the temperature and the third column as the air pressure:

	clim_file = 'path_to_your_file/file_name.txt'

	clim_01 = Boundary('Fourier',**{"file"  : clim_file,
		                        "time"  : "Time (s)",
		                        "T"     : "Temp. (K)",
		                        "P_air" : "Air pressure",
		                        "HR"    : 0.7,
		                        "h_t"   : 10,
		                        "h_m"   : 2e-7 })

where "Time (s)", "Temp. (K)" and "Air pressure" are headers of the text file, on top of each variable of interest.

A boundary can be assigned both constant and variable values at the same time:

* If a key of the dictionary points to a string (like `'T'` and `'P_air'` above), then `Boundary` will attempt to read it under the corresponding header in the .txt file.
* If a key points to a numerical value (like `'HR'`, `'h_t'` and `'h_m'` above), `Boundary` will keep it constant during the entire simulation time.

If any boundary value varies, the dictionary **must** contain a `'file'` key pointing to the file location on your drive, and a `'time'` key pointing to the column containing time coordinates. Note that you can use a single file to store data for both of your domain boundaries, like so:

	clim_file = 'path_to_your_file/file_name.txt'

	clim_01 = Boundary('Fourier',**{"file"  : clim_file,
		                        "time"  : "Time (s)",
		                        "T"     : "Temp. (K)",
		                        "P_air" : "Air pressure",
		                        "HR"    : 0.7,
		                        "h_t"   : 10,
		                        "h_m"   : 2e-7 })

	clim_02 = Boundary('Dirichlet',**{"file" : clim_file,
		                         "time"  : "Time (s)",
		                         "T"     : "Temp. 2 (K)",
		                         "HR"    : 0.85 })

	clim = [clim_01, clim_02]

This formulation aims at giving maximum flexibility for the definition of the boundary conditions.

### Keys

The dictionary which you give as argument for the instantiation of a new `'Boundary'` object may specify more content than in the examples above. The following is a list of the keys that may (or should) be included in it.

* `'T'`: temperature in °C or K
* `'HR'`: relative humidity (dimensionless)
* `'p_v'`: water vapor pressure (Pa), only needed if HR is not given
* `'h_t'`: surface heat transfer coefficient (W/(m2.K)). Optional, default value is 5 W/(m2.K)
* `'h_m'`: surface moisture transfer coefficient (s/m). Optional, default is `7.45e-9 * h_t`
* `'T_eq'`: equivalent temperature, accounting for effects of solar radiation. Optional, it is set equal to `'T'` if not given.
* `'g_l'`: liquid water income caused by rain (kg/(m2.s)). Optional, default value is zero
* `'P_air'`: air pressure (Pa), impacting eventual air transfer in the wall. Optional, default value is zero

Any of these values can either be constant, or read from a .txt file. In this case, these additional keys should be given in the dictionary:

* `'file'`: file location on the hard drive
* `'delimiter'`: delimiter used in the text file. Default is tab.
* `'time'`: location of the time data in the text file

[Documentation main page](../index.md)
