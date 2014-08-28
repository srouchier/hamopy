# -*- coding: utf-8 -*-
"""
Library of some standard materials

@author: Rouchier
"""

from hamopy.classes import Material

################
### Concrete ###
################

concrete = Material('concrete', rho = 2450., cp = 880.)

concrete.set_conduc(lambda_0 = 1.75,
                    lambda_m = 4.5)
                 
concrete.set_isotherm('vangenuchten',**{"w_sat" : 101.44,
                                        "l"     : 1.,
                                        "alpha" : 6.165e-7,
                                        "m"     : 0.219 } )
                                      
concrete.set_perm_vapor('schirmer',**{"mu" : 30.,
                                      "p"  : 0.497 } )
                                    
concrete.set_perm_liquid('durner', **{"K_sat" : 2.2182e-13,
                                      "tau"   : -4.6974,
                                      "l"     : [0.5062, 0.4938],
                                      "alpha" : [5.5383e-7, 2.2493e-8],
                                      "m"     : [0.6148, 0.1913] } )

#############################                       
### Wood fibre insulation ###
#############################

wood_fibre = Material('wood_fibre', rho = 146., cp = 1103.1)

wood_fibre.set_capacity(cp_0 = 1103.1,
                        cp_t = 11.271)

wood_fibre.set_conduc(lambda_0 = 0.038,
                      lambda_m = 0.192,
                      lambda_t = 0.108e-3)
"""                        
wood_fibre.set_isotherm('polynomial', **{"HR" : [0, 0.25, 0.5, 0.75],
                                         "W"  : [0, 6.981, 11.133, 19.299] })
wood_fibre.set_isotherm('slope', **{"HR" : [0.4, 0.6, 0.8],
                                    "XI" : [15.841, 28.686, 59.049] })
"""
wood_fibre.set_isotherm('slope', **{"HR" : [0.25, 0.5, 0.75],
                                    "XI" : [17.704, 20.074, 49.816] })


wood_fibre.set_perm_vapor('interp', **{"HR" : [0.25, 0.75],
                                       "dp" : [3.75e-11, 6.59e-11] } )

wood_fibre.set_perm_air(4.16e-13)