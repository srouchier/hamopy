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

from __future__ import division

import numpy as np
import ham_library as ham
import copy

from scipy.sparse.linalg import spsolve
from scipy.interpolate   import interp1d


def initialisation(init, mesh, clim, thermo = False):
    """
    Initialisation of the array of unknown variables U
    
    If the init dictionary has an 'x' key, the initial conditions are a linear
    interpolation of the values of P and T provided by the user.

    This function is called by hamopy.hamopy1D.calcul only
    """
    
    if init.has_key('x'):
        # The domain is initialised with non-uniform values of P and T
        
        f_t = interp1d(init['x'], init['T'], 'linear')
        T   = f_t(mesh.x)
        
        if init.has_key('HR'):
            f_p = interp1d(init['x'], init['HR'], 'linear')
            HR  = f_p(mesh.x)
            P   = ham.p_c(HR, T)
        elif init.has_key('PV'):
            f_p = interp1d(init['x'], init['PV'], 'linear')
            HR  = f_p(mesh.x) / ham.p_sat(T)
            P   = ham.p_c(HR, T)
        elif init.has_key('PC'):
            f_p = interp1d(init['x'], np.log10(-init['PC']), 'linear')
            P   = 10 ** f_p(mesh.x)
        
    else:
        # The domain is initialised with uniform values of P and T
        
        T = init['T'] * np.ones(mesh.nbr_nodes)
        
        if init.has_key('HR'):
            P = ham.p_c(init['HR'], init['T']) * np.ones(mesh.nbr_nodes)
        elif init.has_key('PV'):
            HR = init['PV'] / ham.p_sat(init['T'])
            P  = ham.p_c(HR, init['T']) * np.ones(mesh.nbr_nodes)
        elif init.has_key('PC'):
            P = init['PC'] * np.ones(mesh.nbr_nodes)
    
    if thermo == True:
        # In case of thermal transfer only
        return T
    else:
        return np.concatenate((P,T))
    
    
def iteration(P, T, S_old, time, mesh, clim):
    """
    This is where things happen.
    
    This function calculates a new iteration of P and T, given the conditions
    provided by the time, mesh and clim objects.
    """
    
    # Construction of the matrices for the linearised equations system
    [C, K]    = mesh.system_matrices(P, T)               # global matrices
    [F, dFdU] = mesh.system_boundary(P, T, clim, time.t) # boundary conditions
    S         = mesh.system_conserv(P, T)                # conservative term
    Q         = mesh.system_airflux(P, T, clim, time.t)  # air transfer source term
    
    # Assembly of the system
    A = C + time.delta*K - time.delta*dFdU
    b = time.delta*(F-Q) - time.delta*K*np.concatenate((P,T)) - (S-S_old)
    
    # Solving
    DELTA_U = spsolve(A, b)
    P_new = P + DELTA_U[:mesh.nbr_nodes]
    T_new = T + DELTA_U[mesh.nbr_nodes:]
    
    return P_new, T_new
    
    
def solution_not_acceptable(P = -5e8, T = 293.15):
    """
    This function raises a flag if the newly calculated values of P or T are
    problematic (either complex, positive or not calculable)
    """
    a = np.any(np.isnan(P)) or np.any(np.iscomplex(P)) or np.any(P>0)
    b = np.any(np.isnan(T)) or np.any(np.iscomplex(T))
    return a or b
    
def calcul(mesh_in, clim_in, init_in, time_in, output_type='dict', logfile = None):
    """
    This is the main algorithm of hamopy.
    
    Required input arguments are the following objects : mesh, clim, init, time
    
    :mesh: A :class:`Mesh` object containing all information on material
           properties and mesh coordinates, along with most methods called for
           the calculation of the new values of P and T at each iteration
    :clim: A list of 2 :class:`Boundary` objects including methods to call the
           value of all field variables as functions of the time of simulation
    :init: A :class:`dict` of initial conditions
    :time: A :class:`Time` object storing information regarding the duration of
           the simulation and the discretisation in time
    
    Optional arguments:
    
    :output_type: String, either 'dict' or 'file'
    :logfile: Specify to save the process of the simulation in a diary
    """
    
    mesh = copy.deepcopy( mesh_in )
    clim = copy.deepcopy( clim_in )
    init = copy.deepcopy( init_in )
    time = copy.deepcopy( time_in )
    
    # U is the array of interest : a concatenation of P and T at each time step
    U = initialisation(init, mesh, clim)
    U = np.asarray([U])
    
    # Initialisation of the simulation diary, if asked by the user
    if logfile is None:
        logobj = None
    else:
        logobj = open(logfile,'w')
        logobj.write("Conv. P \t Conv. T \n")
    
    # Beginning of the time loop
    while time.t <= time.end:
        
        # Initialisation of the new time step
        if logobj is not None:
            logobj.write("t = %s \n" % time.t)
        time.vec.append(time.t)
        nbr_iter = 0
        conv_t   = 1
        conv_p   = 1
        U = np.concatenate( (U, [U[-1,:]]) )
        U_old = copy.deepcopy( U[-1,:] )
        S_old   = mesh.system_conserv(U_old[:mesh.nbr_nodes],U_old[mesh.nbr_nodes:])
        
        # Dirichlet boundary conditions are manually updated
        for i in range(2):
            # Each side of the domain is treated separately
            if clim[i].type == 'Dirichlet':
                ind_P =  i * (mesh.nbr_nodes - 1)
                ind_T =  i * (mesh.nbr_nodes - 1) + mesh.nbr_nodes
                P = ham.p_c(clim[i].HR(time.t), clim[i].T(time.t))
                T = clim[i].T(time.t)
                U[-1, ind_P] = P
                U[-1, ind_T] = T
        
        # Iteration loop of the current time step
        while (conv_t > 1e-4) or (conv_p > 1e-4) or (nbr_iter <= 1) :
            
            # Number of iterations since the time step began
            nbr_iter += 1
            
            # Check that the maximum number of iterations has not been exceeded
            if time.method == 'variable':
                if nbr_iter > time.iter_max:
                    time.try_again(2, logobj)
                    U[-1,:]  = U_old
                    nbr_iter = 0
                    continue
            
            # Key step : calculation of the new P and T profiles
            P = U[-1, :mesh.nbr_nodes ]
            T = U[-1,  mesh.nbr_nodes:]
            [P_new, T_new] = iteration(P, T, S_old, time, mesh, clim)
            
            # Check if the new values of P and T are within acceptable bounds
            if solution_not_acceptable(P = P_new, T = T_new):
                # Divergence : divide time step size by 10 and try again
                time.try_again(10, logobj)
                if time.stop:
                    break
                else:
                    U[-1,:]  = U_old
                    nbr_iter = 0
                    continue
            else:
                conv_p  = np.sum( (P_new-P)**2 ) / np.sum( P**2 )
                conv_t  = np.sum( (T_new-T)**2 ) / np.sum( T**2 )
                U[-1,:] = np.concatenate((P_new,T_new))
                if logobj is not None:
                    logobj.write( "%.6e \t %.6e \n" % (conv_p, conv_t) )
        
        # Time step is over: prepare the next one
        if time.stop:
            break
        
        time.next_step(nbr_iter)
        if time.t > time.end and time.vec[-1] < time.end:
            time.t = time.end
    
    # If the simulation has been interrupted for some reason, return nan
    if time.stop:
        if logobj is not None:
            logobj.write("Convergence not reached: simulation stopped")
            logobj.close()
        print("Convergence not reached: simulation stopped")
        return np.nan
    
    # Otherwise, the simulation has been completed without interruption
    if logobj is not None:
        logobj.close()
    
    # The results are now saved
    if output_type == 'file':
        
        np.savez('hamopy_output',
                 x  = mesh.x,
                 t  = time.vec,
                 T  = U[:,mesh.nbr_nodes:],
                 U  = U,
                 PC = U[:,:mesh.nbr_nodes],
                 HR = ham.HR(U[:,:mesh.nbr_nodes], U[:,mesh.nbr_nodes:]),
                 PV = ham.p_v(U[:,:mesh.nbr_nodes], U[:,mesh.nbr_nodes:])
                 )
                 
    elif output_type == 'dict':

        results = {'x'  : mesh.x,
                   't'  : np.array(time.vec),
                   'U'  : U,
                   'T'  : U[:,mesh.nbr_nodes:],
                   'PC' : U[:,:mesh.nbr_nodes],
                   'HR' : ham.HR(U[:,:mesh.nbr_nodes], U[:,mesh.nbr_nodes:]),
                   'PV' : ham.p_v(U[:,:mesh.nbr_nodes], U[:,mesh.nbr_nodes:])
                    }
                    
        return results


def iteration_thermo(T, T_old, time, mesh, clim):
    """
    This is where things happen.
    
    This function calculates a new iteration of T, given the conditions
    provided by the time, mesh and clim objects.
    """
    
    # Construction of the matrices for the linearised equations system
    [C, K]    = mesh.system_matrices_thermo(T)               # matrices
    [F, dFdU] = mesh.system_boundary_thermo(T, clim, time.t) # boundary
    
    # Newton Raphson
    A = C + time.delta*K - time.delta*dFdU
    b = time.delta*F - time.delta*K*T - C*(T-T_old)
    DELTA_T = spsolve(A, b)
    T_new   = T + DELTA_T

    return T_new

def calcul_thermo(mesh_in, clim_in, init_in, time_in, output_type='dict', logfile = None):
    """
    This is an alternate version of the hamopy algorithm, for heat transfer
    only (without moisture)
    
    Required input arguments are the following objects : mesh, clim, init, time
    
    :mesh: A :class:`Mesh` object containing all information on material
           properties and mesh coordinates, along with most methods called for
           the calculation of the new values of P and T at each iteration
    :clim: A list of 2 :class:`Boundary` objects including methods to call the
           value of all field variables as functions of the time of simulation
    :init: A :class:`dict` of initial conditions
    :time: A :class:`Time` object storing information regarding the duration of
           the simulation and the discretisation in time
    
    Optional arguments:
    
    :output_type: String, either 'dict' or 'file'
    :logfile: Specify to save the process of the simulation in a diary
    """

    mesh = copy.deepcopy( mesh_in )
    clim = copy.deepcopy( clim_in )
    init = copy.deepcopy( init_in )
    time = copy.deepcopy( time_in )
    
    # Initialisation de la matrice des solutions
    T = initialisation(init, mesh, clim, thermo = True)
    T = np.asarray([T])
    
    # Initialisation of the simulation diary, if asked by the user
    if logfile is None:
        logobj = None
    else:
        logobj = open(logfile,'w')
        logobj.write("Conv. T \n")
    
    # Beginning of the time loop
    while time.t <= time.end:
        
        # Initialisation of the new time step
        if logobj is not None:
            logobj.write("t = %s \n" % time.t)
        time.vec.append(time.t)
        nbr_iter = 0
        conv_t   = 1
        T = np.concatenate( (T, [T[-1]]) )
        T_old = copy.deepcopy( T[-1] )
        
        # Dirichlet boundary conditions are manually updated
        for i in range(2):
            # Each side of the domain is treated separately
            if clim[i].type == 'Dirichlet':
                ind =  i * (mesh.nbr_nodes - 1)
                T[-1, ind] = clim[i].T(time.t)
        
        # Iteration loop of the current time step
        while (conv_t > 1e-4) or (nbr_iter <= 1) :
            
            # Number of iterations since the time step began
            nbr_iter += 1
            
            # Check that the maximum number of iterations has not been exceeded
            if time.method == 'variable':
                if nbr_iter > time.iter_max:
                    time.try_again(2, logobj)
                    T[-1]  = T_old
                    nbr_iter = 0
                    continue
            
            # Key step : calculation of the new T profile
            T_new = iteration_thermo(T[-1], T_old, time, mesh, clim)
            
            # Check if the new values of T are within acceptable bounds
            if solution_not_acceptable(P = -5e8, T = T_new):
                # Divergence : divide time step size by 10 and try again
                time.try_again(10, logobj)
                if time.stop:
                    break
                else:
                    T[-1]  = T_old
                    nbr_iter = 0
                    continue
                
            else:
                conv_t = np.sum( (T_new-T[-1])**2 ) / np.sum( T[-1]**2 )
                T[-1]  = T_new
                if logobj is not None:
                    logobj.write( "%.6e \n" % conv_t )
        
        # Time step is over: prepare the next one
        if time.stop:
            break
        
        time.next_step(nbr_iter)
        if time.t > time.end and time.vec[-1] < time.end:
            time.t = time.end
    
    # If the simulation has been interrupted for some reason, return nan
    if time.stop:
        if logobj is not None:
            logobj.write("Convergence not reached: simulation stopped")
            logobj.close()
        print("Convergence not reached: simulation stopped")
        return np.nan
        
    # Otherwise, the simulation has been completed without interruption
    if logobj is not None:
        logobj.close()
    
    # The results are now saved
    if output_type == 'file':
        
        np.savez('hamopy_output',
                 x  = mesh.x,
                 t  = time.vec,
                 T  = T )
                 
    elif output_type == 'dict':

        results = {'x'  : mesh.x,
                   't'  : np.array(time.vec),
                   'T'  : T }
                    
        return results