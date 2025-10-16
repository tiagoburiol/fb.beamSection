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
[         -2.00000   ,         0.00000  ],
[         -0.89443   ,        -0.89443  ],
[         -0.89443   ,         0.89443  ],
[          0.00000   ,         0.00000  ],
[          0.00000   ,         1.00000  ],
[          0.00000   ,        -1.00000  ],
[          0.89443   ,         0.89443  ],
[          0.89443   ,        -0.89443  ],
[          2.00000   ,        -0.00000  ]]),

#
# Elements
#

elements = np.array([



[       4   ,      7   ,      5   ,      3  ],





[       8   ,      9   ,      7   ,      4  ],





[       2   ,      6   ,      8   ,      4  ],





[       1   ,      2   ,      4   ,      3  ]]),




) # Fim da malha


