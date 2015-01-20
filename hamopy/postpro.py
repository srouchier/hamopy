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

from scipy.interpolate import interp1d
import numpy  as np
#from . import ham_library as ham

def distribution(result, var, x = None, t = None):
    """
    Distribution of a variable at a given time
    
    :result: dictionary of results provided by the simulation
    :var: string, which variable to extract from result
    :x: coordinates on which the distribution spans
    :t: float, time of the distribution
    """

    if x is None:
        x = result['x']
    
    if t is None:
        t = result['t'][-1]

    try:
        f_t = interp1d(result['t'], result[var], axis=0)
        
    except KeyError:
        
        print("The key %s is not in the dictionary of results" % var)
        return np.array([])
    
    # Distribution of 'var' at the time 't' on all nodes of the mesh
    distr_nodes = f_t(t)
    
    # Interpolate from the node coordinates to the output discretisation
    f_x = interp1d(result['x'], distr_nodes, 'cubic', bounds_error=False)
    
    return f_x(x)


def evolution(result, var, x, t = None):
    """
    Temporal evolution of a variable at a specific location
    
    :result: dictionary of results provided by the simulation
    :var: string, which variable to extract from result
    :x: float, location of the point
    :t: time scale on which to extract the data (optional)
    """

    if t is None:
        t = result['t']
    
    try:
        f_x = interp1d(result['x'], result[var], 'cubic', axis=1)
        
    except KeyError:
        
        print("The key %s is not in the dictionary of results" % var)
        return np.array([])
    
    # Evolution of 'var' at the location x at all time steps of the simulation
    evol_temps_simul = f_x(x)
    
    # Interpolate from the time of simulation to the output time
    f_t = interp1d(result['t'], evol_temps_simul, bounds_error=False)
    
    return f_t(t)


def surface_heat_flow(result, mesh, clim, side, t = None, total = False):
    """
    Temporal evolution of the heat flow from the air towards the boundary. It
    is calculated with the temperature gradient inside the material and is
    appropriate in case of a Dirichlet boundary condition.
    
    The direction of the flow is from the ambiance to the material
    
    :result: dictionary of results provided by the simulation
    :mesh: the mesh
    :clim: the clim list
    :side: 0 or 1, indicates which side of the wall to consider
    :t: time scale on which to extract the data
    :total: bool, False for conduction only (default), True for coupled effects
    
    Returns:
    
    :q_conv: heat flow caused by sole thermal conduction
    :q_total: total heat flow including convection, evaporation and condensation
    """
    
    if side == 0 or side == 'left':
        x  = result['x'][0]
        x2 = result['x'][1]
        i  = 0
    elif side == 1 or side == 'right':
        x  = result['x'][-1]
        x2 = result['x'][-2]
        i  = 1
    
    # Select the material on the chosen side of the wall
    m = mesh.materials[-i]
    
    # Temperature and pressure at the surface node and the next material node
    T_surf,  T_surf2  = [evolution(result, 'T', _, t)  for _ in [x, x2]]
    PC_surf, PC_surf2 = [evolution(result, 'PC', _, t) for _ in [x, x2]]
    PV_surf, PV_surf2 = [evolution(result, 'PV', _, t) for _ in [x, x2]]
    
    # Heat flow caused by sole air flow
    g_air = mesh.C_air * (clim[0].P_air(t) - clim[1].P_air(t))
    q_air = g_air * ham.cp_air * ((T_surf+T_surf2)/2 - 273.15)
    
    # Heat flow caused by sole conduction
    lambda_ = (m.conduc(PC_surf, T_surf)+m.conduc(PC_surf2, T_surf2)) /2
    q_cond  = lambda_ * (T_surf-T_surf2) / np.abs(x-x2)
    
    # In case of pure thermal calculation, we stop here
    if not result.has_key('PV'):
        
        q_conv = 0.
        
    else:
    
        # Otherwise, let's add moisture-induced heat flow
        # Average conductivities between two nodes closest to the surface
        lambda_ = (m.conduc(PC_surf, T_surf)+m.conduc(PC_surf2, T_surf2)) /2
        delta_  = (m.delta_p(PC_surf) + m.delta_p(PC_surf2)) /2
        k_l     = (m.k_l(PC_surf) + m.k_l(PC_surf2)) /2
        
        # Heat flow caused by evaporation or condensation
        q_conv = ham.l_lv * delta_ * (PV_surf-PV_surf2) / np.abs(x-x2)
        q_conv += ham.cp_liq * ((T_surf+T_surf2)/2 - 273.15) \
            * k_l * (PC_surf-PC_surf2) / np.abs(x-x2)
    
    # Output : either conduction only, or all heat transfer with coupling
    if not total:
        return q_cond
    else:
        return q_cond + q_conv + q_air


def heat_flow(result, mesh, clim, x, t = None, total = False):
    """
    Temporal evolution of the heat flow in some point inside the material
    
    The direction of the flow is from left to right
    
    :result: dictionary of results provided by the simulation
    :mesh: the mesh
    :x: where
    :t: time scale on which to extract the data
    :total: bool, False for conduction only (default), True for coupled effects
    
    Returns:
    
    :q_conv: heat flow caused by sole thermal conduction
    :q_total: total heat flow including convection, evaporation and condensation
    """
    
    x1 = x - 0.001
    x2 = x + 0.001
    
    if x1 < 0:
        return surface_heat_flow(result, mesh, clim, 0, t, total)
    elif x2 > sum(mesh.sizes):
        return -surface_heat_flow(result, mesh, clim, 1, t, total)
    
    m1, m2 = [mesh.find_material(_) for _ in [x1, x2]]
    T1,  T2  = [evolution(result, 'T',  _, t) for _ in [x1, x2]]
    PC1, PC2 = [evolution(result, 'PC', _, t) for _ in [x1, x2]]
    PV1, PV2 = [evolution(result, 'PV', _, t) for _ in [x1, x2]]
    
    # Heat flow caused by sole air flow
    g_air = mesh.C_air * (clim[0].P_air(t) - clim[1].P_air(t))
    q_air = g_air * ham.cp_air * ((T1+T2)/2 - 273.15)
    
    # Heat flow caused by sole conduction
    lambda_ = (m1.conduc(PC1, T1) + m2.conduc(PC2, T2)) /2
    q_cond  = lambda_ * (T1-T2) / np.abs(x1-x2)
    
    # In case of pure thermal calculation, we stop here
    if not result.has_key('PV'):
        
        q_conv = 0.
        
    else:
    
        # Otherwise, let's add moisture-induced heat flow       
        # Average conductivities between two nodes closest to the surface
        delta_ = (m1.delta_p(PC1) + m2.delta_p(PC2)) /2
        k_l    = (m1.k_l(PC1) + m2.k_l(PC2)) /2
        
        # Heat flow caused by evaporation or condensation
        q_conv = ham.l_lv * delta_ * (PV1-PV2) / np.abs(x1-x2)
        q_conv += ham.cp_liq * ((T1+T2)/2 - 273.15) \
            * k_l * (PC1-PC2) / np.abs(x1-x2)
    
    # Output : either conduction only, or all heat transfer with coupling
    if not total:
        return q_cond
    else:
        return q_cond + q_conv + q_air


# Deprecated function

def surface_heat_flow_out(result, mesh, clim, side, t = None):
    """
    Temporal evolution of the heat flow from the air towards the boundary. It
    is calculated with the surface transfer coefficient and is NOT appropriate
    in case of a Dirichlet boundary condition.
    
    :result: dictionary of results provided by the simulation
    :mesh: the mesh
    :clim: the list of boundaries
    :side: 0 or 1, indicates which side of the wall to consider
    :t: time scale on which to extract the data
    
    Returns:
    
    :q_conv: heat flow caused by sole thermal conduction
    :q_total: total heat flow including convection, evaporation and condensation
    """
    
    if side == 0 or side == 'left':
        x  = result['x'][0]
        i  = 0
    elif side == 1 or side == 'right':
        x  = result['x'][-1]
        i  = 1
    
    T_surf  = evolution(result, 'T', x, t)
    
    # Convective heat flow caused by air transfer
    g_air  = mesh.C_air * (clim[0].P_air(t) - clim[1].P_air(t))
    q_wind = (-1)**i * ham.cp_air * (clim[i].T(t)-T_surf) * g_air
    
    # Heat flow caused by sole thermal convection
    T_surf  = evolution(result, 'T', x, t)
    q_cond  = clim[i].h_t(t) * ( clim[i].T_eq(t) - T_surf )
    
    # In case of pure thermal calculation, we stop here
    if not result.has_key('PV'):
        
        q_evap = 0
        q_rain = 0
        
    else:
        
        # Otherwise, let's add moisture-related heat flow
        PV_surf = evolution(result, 'PV', x, t)
        q_rain = clim[i].g_l(t) * ham.cp_liq * ( clim[i].T(t) - 273.15 )
        q_evap = clim[i].h_m(t) * ham.l_lv * ( clim[i].p_v(t) - PV_surf )
        
        # Let's modify the wind-related term as well
        q_wind += ham.l_lv * g_air * (-1)**i * (clim[i].p_v(t)/clim[i].T(t) - PV_surf/T_surf ) / (ham.Rv * ham.rho_air)
          
    return (q_cond, q_cond + q_rain + q_evap + q_wind)
