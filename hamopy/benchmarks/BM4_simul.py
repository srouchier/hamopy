"""
4th exercise of the Hamstad benchmark package
"""

from hamopy.classes import Mesh, Boundary, Time
import hamopy.ham_library as ham

# Choice of materials
from hamopy.materials.hamstad import BM4_load, BM4_finishing

# Add the liquid permeability to the load material
import pandas as pd
import numpy  as np
k_file = 'BM4 Perm.txt'
dataK = pd.read_csv(k_file, delimiter='\t')
PC = -10**( np.array(dataK['log(Psuc)']) )
KL =  10**( np.array(dataK['log(K)']) )

BM4_load.set_perm_liquid('interp', **{"PC" : PC[::-1],
                                      "KL" : KL[::-1]})

# Geometry
mesh = Mesh(**{"materials"    : [BM4_load, BM4_load, BM4_load, BM4_finishing, BM4_finishing],
               "sizes"        : [0.01, 0.08, 0.01, 0.005, 0.015],
               "nbr_elements" : [100, 100, 100, 100, 100] })

# Boundary conditions
clim_file = 'BM4 Climate.txt'

clim1 = Boundary('Fourier',**{"file" : clim_file,
                              "time" : "time (s)",
                              "T"    : "Ta,e",
                              "T_eq" : "Teq,e",
                              "p_v"  : "pa,e",
                              "g_l"  : "gl (kg/m2s)",
                              "h_t"  : 25,
                              "h_m"  : 2e-7 })

clim2 = Boundary('Fourier',**{"file" : clim_file,
                              "time" : "time (s)",
                              "T"    : "Teq,i",
                              "p_v"  : "pa,i",
                              "h_t"  : 8,
                              "h_m"  : 3e-8 })

clim = [clim1, clim2]

# Initial conditions
init = {'T'  : 293.15,
        'PC' : -120738829}

time = Time('variable',**{"delta_t"  : 600,
                          "t_max"    : 432000,
                          "iter_max" : 12,
                          "delta_min": 1e-3,
                          "delta_max": 600 } )

if __name__ == "__main__":
    
    diary = 'hamopy_log'
    
    # Calculation
    from hamopy.algorithm import calcul
    result = calcul(mesh, clim, init, time, logfile = diary)
    
    from hamopy.postpro import evolution
    data0 = pd.read_csv(clim_file, delimiter='\t')
    t_plot = result['t']
    x_plot = [0., 0.1, 0.12]
    
    Temperature = np.column_stack([evolution(result, 'T', _, t_plot) for _ in x_plot]) - 273.15
    Humidite    = np.column_stack([evolution(result, 'HR', _, t_plot) for _ in x_plot])
    TeneurEnEau1 = BM4_load.w(ham.p_c(Humidite[:,0], Temperature[:,0]+273.15), Temperature[:,0]+273.15)
    TeneurEnEau2 = BM4_load.w(ham.p_c(Humidite[:,1], Temperature[:,1]+273.15), Temperature[:,1]+273.15)
    TeneurEnEau3 = BM4_finishing.w(ham.p_c(Humidite[:,1], Temperature[:,1]+273.15), Temperature[:,1]+273.15)
    TeneurEnEau4 = BM4_finishing.w(ham.p_c(Humidite[:,2], Temperature[:,2]+273.15), Temperature[:,2]+273.15)

    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(t_plot/3600, Temperature[:,0])
    plt.figure()
    plt.plot(t_plot/3600, TeneurEnEau1, '-r')
    plt.plot(t_plot/3600, TeneurEnEau2, '--r')
    plt.plot(t_plot/3600, TeneurEnEau3, '-b')
    plt.plot(t_plot/3600, TeneurEnEau4, '--b')
    plt.plot([24, 24], [0, 200], '--k')