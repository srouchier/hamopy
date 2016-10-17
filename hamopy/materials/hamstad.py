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
"""

This library contains the definition of materials involved in the exercises of
the Hamstad benchmark package, formatted so they can be imported into hamopy

An user-defined material library for hamopy may typically look like this file

"""

from hamopy.classes import Material

###################
### Benchmark 1 ###
###################

BM1_load = Material('BM1_load', rho = 2000, cp = 912)

BM1_load.set_conduc(lambda_0 = 1.5, lambda_m = 15.8)

BM1_load.set_isotherm('vangenuchten', **{"w_sat" : 146,
                                         "l"     : 1.,
                                         "alpha" : 8e-8,
                                         "m"     : 0.375 })
                                         
BM1_load.set_perm_vapor('schirmer', **{"mu" : 200,
                                       "p" : 0.497 })
                                       
BM1_load.set_perm_liquid('exp2', **{"a" : [-39.2619, 0.0704, -1.742e-4, -2.7953e-6, -1.1566e-7, 2.5969e-9],
                                    "w0" : 73 })

BM1_insulation = Material('BM1_insulation', rho = 100, cp = 739)

BM1_insulation.set_conduc(lambda_0 = 0.033, lambda_m = 0.59)

BM1_insulation.set_isotherm('vangenuchten', **{"w_sat" : 900,
                                               "l"     : 1.,
                                               "alpha" : 2e-4,
                                               "m"     : 0.5 })

BM1_insulation.set_perm_vapor('schirmer', **{"mu" : 9.6,
                                             "p" : 0.497 })

BM1_seal= Material('BM1_seal', rho=100, cp=100)

BM1_seal.set_conduc(lambda_0 = 100)

BM1_seal.set_isotherm('polynomial', **{"HR" : [0, 0.25, 0.5, 0.75],
                                       "W"  : [0, 1, 2, 3] })


###################
### Benchmark 3 ###
###################

BM3 = Material('BM3', rho = 212, cp = 1000)

BM3.set_conduc(lambda_0 = 0.06,
               lambda_m = 0.56)

BM3.set_isotherm('vangenuchten', **{"w_sat"    : 871,
                                    "l"        : [0.41, 0.59],
                                    "alpha"    : [6.12e-7, 1.22e-6],
                                    "m"        : [0.5981, 0.5816]  })

BM3.set_perm_vapor('schirmer', **{"mu" : 5.6,
                                  "p"  : 0.2 })

BM3.set_perm_liquid('exp', **{"a" : [-46.245, 294.506, -1439, 3249, -3370, 1305] } )

BM3.set_perm_air(3e-5*1.8e-5*0.2)

###################
### Benchmark 4 ###
###################

BM4_load = Material('BM4_load', rho = 2005., cp = 840.)

BM4_load.set_conduc(lambda_0 = 0.5,
                    lambda_m = 4.5)
                    
BM4_load.set_isotherm('vangenuchten',**{"w_sat" : 157.,
                                        "l"     : [0.3, 0.7],
                                        "alpha" : [1.25e-5, 1.8e-5],
                                        "m"     : [0.394, 0.833] })
                                         
BM4_load.set_perm_vapor('schirmer',**{"mu" : 30.,
                                      "p"  : 0.497 } )


BM4_finishing = Material('BM4_finishing', rho = 790., cp = 870.)

BM4_finishing.set_conduc(lambda_0 = 0.2,
                         lambda_m = 4.5)

BM4_finishing.set_isotherm('vangenuchten',**{"w_sat" : 209.,
                                             "l"     : 1.,
                                             "alpha" : 2e-6,
                                             "m"     : 0.2126 })

BM4_finishing.set_perm_vapor('schirmer',**{"mu" : 3.,
                                           "p"  : 0.497 } )

BM4_finishing.set_perm_liquid('exp2', **{"a" : [-33, 0.0704, -1.742e-4, -2.7953e-6, -1.1566e-7, 2.5969e-8],
                                         "w0" : 120 })

###################
### Benchmark 5 ###
###################

BM5_brick = Material('BM5_brick', rho = 1600, cp = 1000)

BM5_brick.set_conduc(lambda_0 = 0.682)

BM5_brick.set_isotherm('vangenuchten',**{"w_sat" : 373.5,
                                         "l"     : [0.46, 0.54],
                                         "alpha" : [4.796e-5, 2.041e-5],
                                         "m"     : [0.333, 0.737] })

BM5_brick.set_perm_vapor('schirmer',**{"mu" : 7.5,
                                       "p"  : 0.2 } )

BM5_brick.set_perm_liquid('exp', **{"a" : [-36.484, 461.325, -5240, 2.907e4, -7.41e4, 6.997e4] } )


BM5_mortar = Material('BM5_mortar', rho = 230, cp = 920)

BM5_mortar.set_conduc(lambda_0 = 0.6, lambda_m = 0.56)

BM5_mortar.set_isotherm('vangenuchten',**{"w_sat" : 700,
                                          "l"     : [0.2, 0.8],
                                          "alpha" : [5.102e-5, 4.082e-7],
                                          "m"     : [0.333, 0.737] })

BM5_mortar.set_perm_vapor('schirmer',**{"mu" : 50,
                                        "p"  : 0.2 } )

BM5_mortar.set_perm_liquid('exp', **{"a" : [-40.425, 83.319, -175.961, 123.863] } )


BM5_insulation = Material('BM5_insulation', rho = 212, cp = 1000)

BM5_insulation.set_conduc(lambda_0 = 0.06, lambda_m = 0.56)

BM5_insulation.set_isotherm('vangenuchten',**{"w_sat" : 871,
                                              "l"     : [0.41, 0.59],
                                              "alpha" : [6.122e-7, 1.224e-6],
                                              "m"     : [0.6, 0.5833] })

BM5_insulation.set_perm_vapor('schirmer',**{"mu" : 5.6,
                                            "p"  : 0.2 } )

BM5_insulation.set_perm_liquid('exp', **{"a" : [-46.245, 294.506, -1439, 3249, -3370, 1305] } )




