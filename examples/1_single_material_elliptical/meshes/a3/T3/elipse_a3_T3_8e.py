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
[          3.00000   ,         0.00000  ],
[          1.65908   ,        -0.83316  ],
[          1.65908   ,         0.83316  ],
[          0.00000   ,         0.00000  ],
[          0.00000   ,        -1.00000  ],
[          0.00000   ,         1.00000  ],
[         -1.65907   ,        -0.83317  ],
[         -1.65907   ,         0.83317  ],
[         -2.99999   ,         0.00000  ]]),
#
# Elements
#
elements = np.array([
[       8   ,      9   ,      4   ],
[       6   ,      8   ,      4   ],
[       3   ,      6   ,      4   ],
[       1   ,      3   ,      4   ],
[       2   ,      1   ,      4   ],
[       5   ,      2   ,      4   ],
[       7   ,      5   ,      4   ],
[       9   ,      7   ,      4   ]]),

)# fim malha
