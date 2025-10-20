*loop elems
*if(elemsmat==0)
*messagebox - An element without any material was found, some posibilities are that you forget to asign a material to the problem or there is a repited entity over another - Process finished with error
*endif
*end
*loop nodes
*if(NodesCoord(3)==0)
*else
*messagebox - This is a plane problem so z coordinate must be =0 for all nodes - Process finished with error
*endif
*end
#=======================================================================
# MAT-fem 1.0  - MAT-fem is a learning tool for undestanding 
#                the Finite Element Method with MATLAB and GiD
#=======================================================================
import numpy as np

mesh_data = dict(
#
#  Material Properties
#  [young, nu, dens]
#
materials = np.array([
*loop elems
*format " %6i  %6i  %6i "
*set var material=0
*loop materials
*set var material=Operation(material(int)+1)
*if(elemsmat == material)
*if(ElemsNum == nelem)
[ *MatProp(1) , *MatProp(2) , *MatProp(3) ] ]),
*else
[ *MatProp(1) , *MatProp(2) , *MatProp(3) ],
*endif
*endif
*end 
*end 
#
# Coordinates
#
coordinates = np.array([
*loop nodes
*format " %15.5f  %15.5f "
*if(NodesNum == npoin)
[ 0, *NodesCoord(1) , *NodesCoord(2) ]]),
*else
[ 0, *NodesCoord(1) , *NodesCoord(2) ],
*endif
*end nodes
#
# Elements
#
elements = np.array([
*loop elems
*format " %6i  %6i  %6i  %6i  %6i  %6i  %6i  %6i  %6i "

*if(nnode == 3)
*if(ElemsNum == nelem)
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) ]]) -1,
*else
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) ],
*endif
*endif

*if(nnode == 6)
*if(ElemsNum == nelem)
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4) , *ElemsConec(5) , *ElemsConec(6) ]]) -1,
*else
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4) , *ElemsConec(5) , *ElemsConec(6) ],
*endif
*endif

*if(nnode == 4)
*if(ElemsNum == nelem)
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4)]]) -1,
*else
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4)],
*endif
*endif

*if(nnode == 8)
*if(ElemsNum == nelem)
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4) , *ElemsConec(5) , *ElemsConec(6) , *ElemsConec(7) , *ElemsConec(8) ]]) -1,
*else
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4) , *ElemsConec(5) , *ElemsConec(6) , *ElemsConec(7) , *ElemsConec(8) ], 
*endif
*endif

*if(nnode == 9)
*if(ElemsNum == nelem)
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4) , *ElemsConec(5) , *ElemsConec(6) , *ElemsConec(7) , *ElemsConec(8) , *ElemsConec(9) ]]) -1,
*else
[ *ElemsConec(1) , *ElemsConec(2) , *ElemsConec(3) , *ElemsConec(4) , *ElemsConec(5) , *ElemsConec(6) , *ElemsConec(7) , *ElemsConec(8) , *ElemsConec(9) ], 
*endif
*endif
*end elems


) # Fim da malha


