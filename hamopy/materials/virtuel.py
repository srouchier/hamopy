# -*- coding: utf-8 -*-
"""
Librairie de matériaux virtuels créés pour un benchmark numérique
"""
import copy
from hamopy.classes import Material

# Isolant hygroscopique proche de la fibre de bois

isolant = Material('isolant', rho = 200., cp = 2000.)


isolant.set_conduc(lambda_0 = 0.05,
                   lambda_m = 0.5,
                   lambda_t = 0.1e-3)

isolant.set_isotherm('slope', **{"HR" : [0.25, 0.5, 0.75],
                                 "XI" : [17, 19, 47] })
                                      
isolant.set_perm_vapor('interp', **{"HR" : [0.25, 0.75],
                                    "dp" : [5e-11, 1e-10] } )

# Mur porteur (béton ou brique)

porteur = Material('porteur', rho = 2300., cp = 900.)

porteur.set_conduc(lambda_0 = 1.5,
                   lambda_m = 2.,
                   lambda_t = 0.1e-3)

porteur.set_isotherm('vangenuchten',**{"w_sat" : 120,
                                       "l"     : 1,
                                       "alpha" : 6.2e-7,
                                       "m"     : 0.22} )
                                      

porteur.set_perm_vapor('interp', **{"HR" : [0.25, 0.75],
                                    "dp" : [5e-12, 1e-11] } )

porteur2 = copy.deepcopy(porteur)

porteur2.set_isotherm('polynomial', **{"HR" : [0, 0.25, 0.5, 0.75],
                                       "W"  : [0, 31.3, 38.1, 48.7] })