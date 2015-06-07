# Library of functions

All material-independant values, universal constants and state equations are stored within the `hamopy/ham_library` file. In addition to constant properties such as the latent heat of evaporation, the universal gas constant and the thermal conductivity of air, the file includes 5 useful methods:

* `p_sat(T)`: Water vapor saturation pressure (Pa)
* `D_va(T)`: Water vapor diffusivity in air (m2/s)
* `p_v(p_c, T)`: Water vapor pressure (Pa)
* `HR(p_c, T)`: Relative humidity
* `p_c(HR, T)`: Capillary pressure (Pa)

The last three methods are equivalent formulations of the Clausius-Clapeyron equation. It is used to calculate the relative humidity from values of the capillary pressure, and inversely.

Here is an example of use:

	from hamopy import ham_library as ham

	ham.cp_liq
	ham.p_sat(280)
	ham.p_c(0.92, 299)

This script first returns the value of the liquid water specific heat, then the water vapor saturation pressure at 280 K, then the capillary pressure at 92 %RH and 299 K.

[Documentation main page](../index.md)
