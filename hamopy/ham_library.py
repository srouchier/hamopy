# Copyright (C) 2014  Simon Rouchier
# 
# This file is part of Hamopy.
#
# Hamopy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Hamopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Hamopy.  If not, see <http://www.gnu.org/licenses/>.
# 
# To contact the author : s.rouchier@gmail.com

from pylab import exp, log

# Properties regarding temperature
R       = 8.314     # ideal gas constant
T_0     = 273.15    # zero Celsius
T_ref   = 293.15    # a reference temperature [K]

# Properties of air
p_atm       = 101325.       # atmospheric pressure [Pa]
lambda_air  = 0.0262        # thermal conductivity [W.m-1.K-1]
alpha_air   = 2.216e-5      # thermal diffusivity [m2.s-1]
rho_air     = 1.2           # density [kg.m-3]
mu_air      = 1.8e-5        # dynamic viscosity [Pa.s]
cp_air      = 1004.         # specific heat [J.kg-1.K-1]

# Properties of water
rho_liq     = 998.          # density [kg.m-3]
eta_liq     = 1e-3          # dynamic viscosity [Pa.s]
cp_liq      = 4180.         # heat capacity [J.kg-1.K-1]
lambda_liq  = 0.6           # thermal conductivity [W.m-1.K-1]
tension     = 72.7e-3       # water/air surface tension [N.m-1]

# Properties of water vapour
cp_vap = 1850.              # specific heat [J.kg-1]
l_lv   = 2.5e6              # latent heat of evaporation at 0 C [J.kg-1]
Rv     = R/18 * 1e3         # specific gas constant [J.kg-1.K-1]

# FONCTIONS

def p_sat(T = 293.15):
    """
    Water vapor saturation pressure [Pa]
    
    input : temperature T [K]
    """
    
    T_celsius = T - 273.15
    log_p_sat = 2.7858 + 7.5*T_celsius / (237.3+T_celsius)
    
    return 10**log_p_sat
    
    
def D_va(T):
    """
    Water vapor diffusivity in air [m2/s]
    
    input : temperature T [K]
    """
    return 2.306e-5 * (T/273.15)**1.81;

def p_v(p_c, T):
    """
    Clausius-Clapeyron formula: returns the value of water vapor pressure as a
    function of the capillary pressure and temperature
    
    input : capillary pressure p_c [Pa], temperature T [K]
    """
    return p_sat(T) * exp( p_c / (rho_liq*Rv*T) )
    
def HR(p_c, T = 293.15):
    """
    Clausius-Clapeyron formula: returns the value of relative humidity as a
    function of the capillary pressure and temperature
    
    input : capillary pressure p_c [Pa], temperature T [K]
    """
    return exp( p_c / (rho_liq*Rv*T) )

def p_c(HR, T = 293.15):
    """
    Clausius-Clapeyron formula: returns the value of capillary pressure as a
    function of the relative humidity and temperature
    
    input : capillary pressure p_c [Pa], temperature T [K]
    """
    return rho_liq * Rv * T * log(HR)
