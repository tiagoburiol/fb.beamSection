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
[          0.00000   ,        -0.00000  ],
[          0.00000   ,         1.00000  ],
[          0.00000   ,        -1.00000  ],
[          2.62048   ,        -0.85166  ],
[          2.62048   ,         0.85166  ],
[          5.00000   ,         0.00000  ]]),

#
# Elements
#

elements = np.array([



[       1   ,      3   ,      4   ,      2  ],





[       5   ,      2   ,      4   ,      8  ],





[       6   ,      7   ,      4   ,      3  ],





[       8   ,      4   ,      7   ,      9  ]]),




) # Fim da malha


