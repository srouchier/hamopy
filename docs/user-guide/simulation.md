#Running a simulation

Let's say you have correctly defined all conditions of the simulation by reading [this helpful page](inputs.md). All you need to do now is import the main algorithm and run it

	from hamopy.algorithm import calcul

	results = calcul(mesh, clim, init, time)

The rest is [post-processing](post-processing.md).

## Options

The syntax for calling the simulation may be slightly modified for three purposes: choosing the type of output, keeping a record of the simulation process, or simulate thermal transfer only.

### Output type

The `calcul()` method normally returns a dictionary storing all necessary data for [post-processing](post-processing.md). Another available option is to save this data into a file in the uncompressed .npz format

	from hamopy.algorithm import calcul

	calcul(mesh, clim, init, time, output_type = 'file')

This will create a file called `hamopy_output.npz` in the working directory. This method of output selection can be customised by modifying the `algorithm` module.

### Diary

Should you encounter convergence difficulties, you may want the algorithm to store its progress somewhere so you may locate the problems. The `calcul()` method may take an additional argument for this purpose:

	from hamopy.algorithm import calcul

	diary = 'file_name'

	results = calcul(mesh, clim, init, time, logfile = diary)

This will create a file at the specified location, on which the record of the convergence criteria will be saved.

### Heat transfer only

The code is initially designed for coupled heat and moisture transfer, but the user may want to skip the humidity and only calculate heat. This is done through an alternate version of the algorithm:

	from hamopy.algorithm import calcul_thermo

	results = calcul_thermo(mesh, clim, init, time)

This should run about 4 times faster than coupled heat and mass transfer.

[Documentation main page](../index.md)
