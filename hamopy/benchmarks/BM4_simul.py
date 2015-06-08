"""
4th exercise of the Hamstad benchmark package
"""

from hamopy.classes import Mesh, Boundary, Time

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
mesh = Mesh(**{"materials"    : [BM4_load, BM4_finishing],
               "sizes"        : [0.1, 0.02],
               "nbr_elements" : [200, 200] })

# Boundary conditions
clim_file = 'D:\MCF\Simulation\Python\Hamstad/BM4 Climate.txt'

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
    
    diary = 'D:\MCF\Simulation\Python\hamopy_log'
    
    # Calculation
    from hamopy.algorithm import calcul
    result = calcul(mesh, clim, init, time, logfile = diary)
