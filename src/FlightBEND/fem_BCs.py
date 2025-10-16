'''Module containing objects related to boundary conditions, 
like displacements, forces, moments, etc.'''
from .fem_geometry  import *
from .fem_BCs       import *
import numpy as np

## ================================================================ ##
##                      FORCE CLASSES                               ##
## ---------------------------------------------------------------- ##
class Force:
    '''Basic force class'''
    def __init__(self, direction):
        # Direction of application
        self._direction = direction

    ## --------------------- Setters --------------------- ##
    def setDirection(self, direction:tuple[float,float,float]):
        '''Changes the force direction'''
        self._direction = direction
        
    ## --------------------- Getters --------------------- ##
    def getDirection(self) -> tuple[float,float,float]:
        '''Gets the force direction'''
        return self._direction

class MomentCouple(Force):
    '''Subclass for moment couple. Applied at a node'''
    def __init__(self, components:tuple[float,float,float], nodeIndex:int):
        ## Inheritance (null direction)
        Force.__init__(self, [0,0,0])
        # Components
        self._components = components
        # Magnitude
        self._magnitude = np.linalg.norm(components)
        # Updating the direction
        direction = np.array(components)/self._magnitude
        self._direction = direction
        # Node where it is applied
        self._nodeIndex = nodeIndex
    def __repr__(self):
        return f'Couple,\n  Node: {self.getNode()},\n  Mag.: {self.getMagnitude()},\n  Dir.: {self.getDirection()}'
    def __str__(self):
        return 'MomentCouple'
    
    ## --------------------- Setters --------------------- ##
    def setMagnitude(self, magnitude:float):
        '''Changes the magnitude of the force'''
        self._magnitude = magnitude
    def setNode(self, node:Node):
        '''Changes the node where the force is applied'''
        self._node = node
    def setComponents(self, comps:tuple[float,float,float]):
        '''Changes the components and direction of the force'''
        # Updating the components
        self._components = np.hstack([comps, [0,0,0]])

        # Unit vector in the direction of the component
        v   = np.array(comps)
        unit = v/np.linalg.norm(v)

        # Updating the direction of the displacement
        self._direction = unit
    ## --------------------- Getters --------------------- ##
    def getMagnitude(self) -> float:
        '''Gets the force magnitude'''
        return self._magnitude
    def getNode(self) -> Node:
        '''Gets the node where the force is applied'''
        return self._nodeIndex
    def getComponents(self):
        return self._components

## Point Load
class PointLoad(Force):
    '''
    Subclass for point load. Applied at a node
    PARAMETERS:
        components: Force components (in the local axes of the reference element)
        pos: Application point on the cross-section (relative to the section geometry's reference frame [Z_arb,Z_arb])
        nodeIndex: Index of the Node
        refElement: reference element to get the local axes where the comp. are written
    
    '''
    def __init__(self, components:tuple[float,float,float], pos:tuple[float,float], nodeIndex:int, refElement:Elem1D):
        ## Inheritance (null direction)
        Force.__init__(self, [0,0,0])
        # Magnitude
        self._magnitude = np.linalg.norm(components)
        # Application point on the cross-section
        self._pos = np.array(pos)
        # Reference element
        self._refElement = refElement
        # Direction cosine matrix of the ref element
        R = refElement._localAxes
        # Components in the global reference frame
        F_local           = components              # Components in the local system of the ref element
        F_global          = R@F_local                 # Components in the global system
        self._globalDir   = F_global/self._magnitude  # Direction in the global system
        self._globalComp = F_global
        self._localComp   = components                # Local components
        
        # Node where it is applied
        self._nodeIndex = nodeIndex
        
        # Updating the direction
        self._direction = np.array(components)/self._magnitude
        
        
        
    def __repr__(self):
        return f'Point,\n  Node: {self.getNode()},\n  Mag.: {self.getMagnitude()},\n  Dir.: {self.getDirection()}'
    def __str__(self):
        return 'PointLoad'
    
    ## --------------------- Setters --------------------- ##
    def setMagnitude(self, magnitude:float):
        '''Changes the magnitude of the force'''
        self._magnitude = magnitude
    def setNode(self, node:Node):
        '''Changes the node where the force is applied'''
        self._node = node
    def setComponents(self, comps:tuple[float,float,float]):
        '''Changes the components and direction of the force'''
        # Updating the components
        self._components = np.hstack([comps, [0,0,0]])

        # Unit vector in the direction of the component
        v   = np.array(comps)
        unit = v/np.linalg.norm(v)

        # Updating the direction of the displacement
        self._direction = unit
    ## --------------------- Getters --------------------- ##
    def getMagnitude(self) -> float:
        '''Gets the force magnitude'''
        return self._magnitude
    def getNode(self) -> Node:
        '''Gets the node where the force is applied'''
        return self._nodeIndex
    def getComponents(self):
        return self._components

## Distributed Load
class DistributedLoad(Force):
    '''Subclass for distributed load. Applied on an element, meaning the magnitudes applied
    to each node are provided'''
    def __init__(self, direction, magnitudes:list[float,float], elem:Elem):
        ## Inheritance
        Force.__init__(self, direction)
        # Force magnitude at each node
        self._magnitudes = magnitudes if type(magnitudes) == list else list(magnitudes)
        # Element where it is applied
        self._element = elem
    def __repr__(self):
        return f'Distributed,\n  Elem.: {self.getElem()},\n  Mags.: {self.getMagnitude()},\n  Dir.: {self.getDirection()}'
    def __str__(self):
        return 'DistributedLoad'

    ## --------------------- Setters --------------------- ##
    def setMagnitude(self, magnitudes:list[float,float]):
        '''Changes the force magnitude'''
        self._magnitudes = magnitudes
    def setElem(self, elem:Elem):
        '''Changes the element where the force is applied'''
        self._element = elem
    ## --------------------- Getters --------------------- ##
    def getMagnitude(self) -> list[float,float]:
        '''Gets the force magnitude'''
        return self._magnitudes
    def getElem(self) -> Elem:
        '''Gets the element where the force is applied'''
        return self._element

## ================================================================ ##
##                   PRESCRIBED DISPLACEMENT CLASS                  ##
## ---------------------------------------------------------------- ##
class Displacement:
    '''
    Generic class for prescribed displacement. Applied at a node, has magnitude and direction.
    '''
    def __init__(self, direction:tuple[float,float,float], magnitude:float, nodeIndex:int):
        # Force magnitude at each node
        self._magnitude = magnitude
        # Node where it is applied
        self._nodeIndex = nodeIndex
        # Direction of application
        self._direction = direction
        # # Prescribed nodes
        # self.prescribedDof = (6*nodeIndex) + np.array([0,1,2,3,4,5])

    ## --------------------- Setters --------------------- ##
    def setDirection(self, direction:tuple[float,float,float]):
        '''Changes the displacement direction'''
        self._direction = direction
    def setMagnitude(self, magnitude:float):
        '''Changes the displacement magnitude'''
        self._magnitude = magnitude
    def setNode(self, node:Node):
        '''Changes the node where the displacement is applied'''
        self._node = node
    
    ## --------------------- Getters --------------------- ##
    def getDirection(self) -> tuple[float,float,float]:
        '''Gets the displacement direction'''
        return self._direction
    def getMagnitude(self) -> float:
        '''Gets the displacement magnitude'''
        return self._magnitude
    def getNode(self) -> Node:
        '''Gets the node where the displacement is applied'''
        return self._node

class LinearDisplacement(Displacement):
    '''Class for prescribed linear displacement. Applied at a node, has magnitude and direction.'''
    def __init__(self, direction:tuple[float,float,float], magnitude:float, nodeIndex:int):
        Displacement.__init__(self, direction, magnitude, nodeIndex)
        u = magnitude*np.array(direction)
        self._components = np.hstack([u, [0,0,0]])
    
    def __repr__(self):
        return f'Ang. Disp.: \n\t Node: {self.getNode()}, \n\t Mag.: {self.getMagnitude()}, \n\t Comps: {self.getComponents()}'
    def __str__(self):
        return 'LinearDisplacement'
    ## --------------------- Setters --------------------- ##
    def setComponents(self, comps:tuple[float,float,float]):
        '''Changes the displacement components and direction'''
        # Updating the components
        self._components = np.hstack([[0,0,0], comps])

        # Unit vector in the direction of the component
        v   = np.array(comps)
        unit = v/np.linalg.norm(v)

        # Updating the displacement direction
        self._direction = unit
    ## --------------------- Getters --------------------- ##
    def getComponents(self):
        return self._components

class AngularDisplacement(Displacement):
    '''Class for prescribed angular displacement. Applied at a node, has magnitude and direction.'''
    def __init__(self, direction:tuple[float,float,float], magnitude:float, nodeIndex:int):
        Displacement.__init__(self, direction, magnitude, nodeIndex)
        u = magnitude*np.array(direction)
        self._components = np.hstack([[0,0,0], u])

    def __repr__(self):
        return f'Ang. Disp.: \n\t Node: {self.getNode()}, \n\t Mag.: {self.getMagnitude()}, \n\t Comps: {self.getComponents()}'
    def __str__(self):
        return 'AngularDisplacement'
    ## --------------------- Setters --------------------- ##
    def setComponents(self, comps:tuple[float,float,float]):
        '''Changes the displacement components and direction'''
        # Updating the components
        self._components = np.hstack([[0,0,0], comps])

        # Unit vector in the direction of the component
        v   = np.array(comps)
        unit = v/np.linalg.norm(v)

        # Updating the displacement direction
        self._direction = unit
    ## --------------------- Getters --------------------- ##
    def getComponents(self):
        return self._components
    
class FixedSupport(Displacement):
    ''''Class for fixed support, zero displacements and rotations'''
    def __init__(self, nodeIndex:int):
        Displacement.__init__(self, [0,0,0], 0, nodeIndex)
        self._components = np.array([0,0,0,0,0,0])
    def __str__(self):
        return 'FixedSupport'