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
[         -4.00000   ,        -0.00000  ],
[         -2.13691   ,        -0.84534  ],
[         -2.13691   ,         0.84534  ],
[          0.00000   ,         0.00000  ],
[          0.00000   ,        -1.00000  ],
[          0.00000   ,         1.00000  ],
[          2.13691   ,        -0.84534  ],
[          2.13691   ,         0.84534  ],
[          4.00000   ,         0.00000  ]]),

#
# Elements
#

elements = np.array([



[       1   ,      2   ,      4   ,      3  ],





[       6   ,      3   ,      4   ,      8  ],





[       8   ,      4   ,      7   ,      9  ],





[       5   ,      7   ,      4   ,      2  ]]),




) # Fim da malha


