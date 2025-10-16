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
[         -1.00000   ,         0.00000  ],
[         -0.70711   ,         0.70711  ],
[         -0.70711   ,        -0.70711  ],
[          0.00000   ,         0.00000  ],
[          0.00000   ,         1.00000  ],
[          0.00000   ,        -1.00000  ],
[          0.70711   ,         0.70711  ],
[          0.70711   ,        -0.70711  ],
[          1.00000   ,         0.00000  ]]),

#
# Elements
#

elements = np.array([



[       2   ,      1   ,      4   ,      5  ],





[       7   ,      5   ,      4   ,      9  ],





[       8   ,      9   ,      4   ,      6  ],





[       3   ,      6   ,      4   ,      1  ]]),




) # Fim da malha


