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
[         -3.00000   ,         0.00000  ],
[         -0.94868   ,        -0.94868  ],
[         -0.94868   ,         0.94868  ],
[          0.00000   ,         0.00000  ],
[          0.00000   ,         1.00000  ],
[          0.00000   ,        -1.00000  ],
[          0.94868   ,         0.94868  ],
[          0.94868   ,        -0.94868  ],
[          3.00000   ,         0.00000  ]]),

#
# Elements
#

elements = np.array([



[       1   ,      2   ,      4   ,      3  ],





[       5   ,      3   ,      4   ,      7  ],





[       6   ,      8   ,      4   ,      2  ],





[       7   ,      4   ,      8   ,      9  ]]),




) # Fim da malha


