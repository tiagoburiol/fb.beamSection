#=======================================================================
# MAT-fem 1.0  - MAT-fem is a learning tool for undestanding 
#                the Finite Element Method with MATLAB and GiD
#=======================================================================
# PROBLEM TITLE = Titulo del problema
import numpy as np

malha = dict(
#
#  Material Properties
#
young =   210000000000.00000 ,
poiss =              0.20000 ,


#
# Coordinates
#

coordinates = np.array([
[         -5.00000   ,        -0.00000  ],
[         -2.62048   ,         0.85166  ],
[         -2.62048   ,        -0.85166  ],
[         -1.31024   ,         0.42583  ],
[         -1.31024   ,        -0.42583  ],
[          0.00000   ,        -0.00000  ],
[          0.00000   ,         1.00000  ],
[          0.00000   ,        -1.00000  ],
[          1.31024   ,        -0.42583  ],
[          1.31024   ,         0.42583  ],
[          2.62048   ,        -0.85166  ],
[          2.62048   ,         0.85166  ],
[          5.00000   ,        -0.00000  ]]),

#
# Elements
#

elements = np.array([


[       2   ,      3   ,      6   ,      1   ,      5   ,      4   ],





[       6   ,     11   ,     12   ,      9   ,     13   ,     10   ],





[       2   ,      6   ,     12   ,      4   ,     10   ,      7   ],





[       6   ,      3   ,     11   ,      5   ,      8   ,      9   ]]),





) # Fim da malha


