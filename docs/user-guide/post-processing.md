# Post-processing

Let's say you have correctly defined all simulation [inputs](inputs.md) and successfully performed a [simulation run](simulation.md). The last part of your script may look like this:

	from hamopy.algorithm import calcul

	results = calcul(mesh, clim, init, time)

Now you want to display and analyse some results, which is why you use this program in the first place.

## Simulation outcome

By default, the main algorithm of hamopy returns a dictionary containing all simulation results and everything one needs to interpret them. Here is a list of all keys and values stored within this dictionary.

* `'x'`: Mesh node coordinates; size: `(N,)`
* `'t'`: All time coordinates of the simulation; size: `(M,)`
* `'T'`: Temperature; size: `(M,N)`
* `'PC'`: Capillary pressure; size: `(M,N)`
* `'HR'`: Relative humidity; size: `(M,N)`
* `'PV'`: Vapor pressure; size: `(M,N)`

All values are numpy arrays. This is for instance how to read the temperature of the i^th node at the j^th time of simulation:

	results['T'][j,i]


As this data is a bit raw, two methods are available to easily extract data at user-defined times and locations without having to directly manipulate elements of the `results` dictionary.

### evolution()

The `evolution()` method of the `hamopy.postpro` module helps extract the temporal evolution of a variable at a specific location.

	from hamopy.postpro import evolution
	import numpy as np

	x_out = 0.05
	t_out = np.array([0, 60, 120, 180, 240, 300, 360])
	T_out = evolution(results, 'T', x_out, t_out)

This example returns the evolution of the temperature over time, at the point given by `x_out`, with the temporal discretisation given by `t_out`. The function may take 4 input arguments:

* the dictionary of results, provided by the simulation
* a string denoting which variable to extract (it must be one of the keys of `results`)
* the location of the point (preferably a single value)
* the time scale on which to extract the data (numpy array)

The last argument is optional: if not given, `evolution()` will take all time coordinates in `results['t']` (this is not advised if the simulation time step size was adaptative).

### distribution()

The `distribution()` method of the `hamopy.postpro` module helps extract the spatial distribution of a variable at a specific time.

	from hamopy.postpro import distribution
	import numpy as np

	x_out  = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10])
	t_out  = 3600
	HR_out = distribution(results, 'HR', x_out, t_out)

This example returns the distribution of relative humidity, at the time given by `t_out`, over the spatial discretisation given by `x_out`. The function may take 4 input arguments:

* the dictionary of results, provided by the simulation
* a string denoting which variable to extract (it must be one of the keys of `results`)
* the coordinates on which the distribution spans (numpy array)
* the time of the distribution (preferably a single value)

The third argument is optional: if not given, `distribution()` will take all mesh node coordinates in `results['x']`.

[Documentation main page](../index.md)
