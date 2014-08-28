"""
3rd exercise of the Hamstad benchmark package
"""

from hamopy.classes import Mesh, Boundary, Time

# Choice of materials and geometry
from hamopy.materials.hamstad import BM3

mesh = Mesh(**{"materials"    : [BM3],
               "sizes"        : [0.2],
               "nbr_elements" : [40] })

# Boundary conditions
clim_file = 'D:\MCF\Simulation\Python\Hamstad/BM3 climate.txt'

clim1 = Boundary('Fourier',**{"file"      : clim_file,
                              "time"      : "Time (s)",
                              "T"         : 293.15,
                              "HR"        : 0.7,
                              "h_t"       : 10,
                              "h_m"       : 2e-7,
                              "P_air"     : "DeltaP"})

clim2 = Boundary('Fourier',**{"T"         : 275.15,
                              "HR"        : 0.8,
                              "h_t"       : 10,
                              "h_m"       : 7.38e-12 })
clim = [clim1, clim2]

# Initial conditions
init = {'T'  : 293.15,
        'HR' : 0.95}

# Time step control
time = Time('variable',**{"delta_t"  : 900,
                          "t_max"    : 8640000,
                          "iter_max" : 12,
                          "delta_min": 1e-3,
                          "delta_max": 900 } )

if __name__ == "__main__":
    
    import numpy  as np
    import pandas as pd
    
    # Calculation
    from hamopy.algorithm import calcul
    result = calcul(mesh, clim, init, time)
    
    # Post processing
    from hamopy.postpro import evolution
    data0 = pd.read_csv(clim_file, delimiter='\t')
    t_plot = np.array( data0['Time (s)'] )
    x_plot = [0.05, 0.1, 0.15, 0.17, 0.19]
    
    from hamopy import ham_library as ham
    
    Temperature = np.column_stack([evolution(result, 'T', _, t_plot) for _ in x_plot]) - 273.15
    Humidite    = np.column_stack([evolution(result, 'HR', _, t_plot) for _ in x_plot])
    Eau = BM3.w(ham.p_c(Humidite, Temperature+273.15), Temperature+273.15)
