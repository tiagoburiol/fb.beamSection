#=======================================================================
# MAT-fem 1.0  - MAT-fem is a learning tool for undestanding 
#                the Finite Element Method with MATLAB and GiD
#=======================================================================
# PROBLEM TITLE = Titulo del problema
#
#  Material Properties
#
import numpy as np

malha = dict(
young =   210000000000.00000,
poiss =              0.20000,
denss = 0.00,
pstrs =  1,
thick =              0.10000,
#
# Coordinates
#
coordinates = np.array([
[         -5.00000   ,         0.00000  ],
[         -0.57354   ,         0.99340  ],
[         -0.57354   ,        -0.99340  ],
[          0.00000   ,        -0.00000  ],
[          0.57354   ,        -0.99340  ],
[          0.57354   ,         0.99340  ],
[          5.00000   ,        -0.00000  ]]),
#
# Elements
#
elements = np.array([
[       1   ,      3   ,      2   ],
[       3   ,      4   ,      2   ],
[       4   ,      5   ,      6   ],
[       5   ,      7   ,      6   ],
[       6   ,      2   ,      4   ],
[       4   ,      3   ,      5   ]]),

)# fim malha
