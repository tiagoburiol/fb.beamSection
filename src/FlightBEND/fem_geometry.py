'''Module containing classes related to mesh geometry, 
like nodes, elements, meshes, etc. '''
import numpy as np
from itertools import product                       # To avoid using nested loops
from nptyping import NDArray, Int, Shape, Float     # Use numpy type annotations
from abc import ABC, abstractmethod                 # To define abstract methods in parent classes

__all__ = ['Elem','Elem1D','Elem2D','Mesh','Mesh1D','Node','SectionElem']

## ================================================================ ##
##                            Node Class                            ##
## ================================================================ ##
class Node:
    '''Class for nodes
    '''
    ## Class attributes
    __nodeCount = 0     # Counter for the current number of class instances
    
    ## --------------------- Constructor --------------------- ##
    def __init__(self, coords:np.ndarray, index:int = None):
        if type(coords) in [int, float]:
            raise Exception('Coodinates must be an array or list')
        
        # Converts input to an array
        if type(coords) is not np.ndarray:
            self._coords = np.array(coords)
        else:
            self._coords = coords

        # Associates the node count with the node index
        if index is None:
            self._index = Node.__nodeCount
        else:
            # In case an index is given
            self._index = index

        # Increments the node counter
        Node.__nodeCount += 1     
    def __repr__(self):
        return str(np.round(self.getCoords(),3))

    ## --------------------- Getters --------------------- ##
    def getCoords(self) -> np.ndarray:
        # If the node is in 2D space, set x = 0
        if len(self._coords) == 2:
            return np.insert(self._coords,0,0) 
        return self._coords
    def getIndex(self) -> int:
        return self._index


## ================================================================ ##
##                          ELEMENT CLASSES                         ##
## ---------------------------------------------------------------- ##
class Elem:
    '''Generic class for elements'''
    ## --------------------- Constructor --------------------- ##
    def __init__(self, nodes:tuple[Node,Node], index:int):
        self._nodes = nodes    # Element nodes (connectivity)
        self._index = index    # Element index
    def __repr__(self):
        nodes  = self.getNodeIndices()
        index = self.getIndex()
        return 'Index: {}, Nodes: {}'.format(index,nodes)

    ## --------------------- Getters --------------------- ##
    def getNodes(self)          -> tuple[Node, Node]:
        '''Returns a tuple containing the element's nodes'''
        return self._nodes
    def getNodeIndices(self)    -> tuple[Node, Node]:
        '''Returns a tuple containing the indices of the element's nodes'''
        return [node._index for node in self._nodes]
    def getNodeCoords(self)     -> list:
        '''Returns a list containing the coordinates of the element's nodes'''            
        return list(node.getCoords() for node in self.getNodes())
    def getNCs(self)            -> list:
        '''Wrapper for the ``getNodeCoords()`` method.
        Returns a list containing the coordinates of the element's nodes. 
        '''
        # Keeping 'getNCs' as it might be an abbreviation (Node Coordinates) used elsewhere, 
        # but renaming the method in the docstring for clarity.
        # return list(node.getCoords() for node in self.getNodes())
        return self.getNodeCoords()
    def getIndex(self)          -> int:
        '''Get the element's index'''
        return self._index
    
## ---------------------------------------------------------------- ##
##                            1D Elements                           ##
## ---------------------------------------------------------------- ##
class Elem1D(Elem):
    '''
    General class for 1D (line) finite elements.  \n
    
    More specific element types like bars, or beams extend this class 
    with more specific methods for their stiffness matrices. 
    '''
    ## --------------------- Constructor --------------------- ##
    def __init__(self, nodes:tuple[Node,Node], index:int, rot_x:float=0):
        ## Inheritance from the element class
        Elem.__init__(self, nodes, index)

        
        ## ___ Initializations
        self._rot_x     = rot_x                          # [rad] Rotation around x
        self._length    = np.nan
        self._localAxes = np.nan * np.ones([3,3])    # vector of 3 vectors: [ex, ey, ez]
        self._midPoint  = np.nan


        ## ___ Local Variables
        p1, p2 = self.getNodeCoords()  # getting the node coordinates
        #p1, p2 = nos[0], nos[1]

        ## ___ Computations
        self._length    = self.__calcLength__(p1,p2)     # computing the length 
        self._localAxes = self.__calcLocalAxes__(p1,p2)  # computing the local axes
        self._midPoint  = self.__calcMidPoint__(p1,p2)   # computing the midpoint

        # If rotation around x is prescribed
        if rot_x != 0:
            self.setRotation(rot_x)            
    def __repr__(self):
        nodes  = self.getNodeIndices()
        length = self.getLength()
        rot    = self.getRotation()
        return 'Nodes: {}, Length: {:.4e}, x rotation: {:.4e}'.format(nodes,length,rot)

    ## --------------------- Calculations --------------------- ##
    def __calcLength__(self,p1,p2)      -> float:
        """
        Calculates the element length.

        Parameters
        ----------
        p1 : array
            Coordinates of node 1.
        p2 : array
            Coordinates of node 2.

        Returns
        -------
        L : float
            Element length.
        """
        # Length 
        L  = np.linalg.norm(p2-p1)
        return L
    def __calcLocalAxes__(self,p1,p2)   -> np.ndarray:
        """
        Calculates the element's local coordinate system axes.

        Parameters
        ----------
        p1 : array
            Coordinates of node 1.
        p2 : array
            Coordinates of node 2.

        Returns
        -------
        local_axes : np.ndarray
            A 3x3 matrix where rows are the local unit vectors [i, j, k].
        """
        ## Define global direction vectors
        # I = np.array([1,0,0])
        J = np.array([0,1,0])
        K = np.array([0,0,1])
        
        ## Define local x direction (i)
        i = (p2 - p1) / self._length
        # Direction cosines of the local x-axis
        # cos_11 = np.sum(i*I)
        cos_12 = np.sum(i*J)
        cos_13 = np.sum(i*K)
        # Direction angles of the local x-axis
        # lamb_11 = np.arccos(cos_11)
        lamb_12 = np.arccos(cos_12)
        lamb_13 = np.arccos(cos_13)

        ## Define local y direction as K x i
        prod  = np.cross(K,i)
        norma = np.linalg.norm(prod)
        # Check if e2 is parallel to the global Z axis
        if norma < 1e-6 or np.isnan(norma):
            # If it is, define local y as J x i
            j = np.cross(J,i)/np.sin(lamb_12)
        else:
            j = prod/np.sin(lamb_13)

        ## Define local z direction as i x j
        k = np.cross(i,j)
        return np.array([i,j,k])  
    def __calcMidPoint__(self,p1,p2)    -> np.ndarray:
        """
        Calculates the coordinates of the element's midpoint.

        Parameters
        ----------
        p1 : array
            Coordinates of node 1.
        p2 : array
            Coordinates of node 2.

        Returns
        -------
        midpoint : np.ndarray
            Coordinates of the midpoint.
        """
        # Midpoint
        dp = p2-p1
        return p1 + np.array(dp)/2
        

    ## --------------------- Getters --------------------- ##
    def getLength(self)     -> float:
        '''
        Returns
        -------
        L : float
            Element length.
        '''
        return self._length
    def getLocalAxes(self)  -> list[np.ndarray,np.ndarray,np.ndarray]:
        '''
        Returns
        -------
        local_axes : list[np.ndarray,np.ndarray,np.ndarray]
            List of direction vectors for the element's local axes [ex, ey, ez].
        '''
        return self._localAxes
    def getMidPoint(self)   -> np.ndarray:
        '''
        Returns
        -------
        midpoint : np.ndarray
            Coordinates of the element's midpoint [mx, my, mx].
        '''
        return self._midPoint
    def getJacob(self)      -> list[float,np.ndarray]:
        '''
        Calculates the element's Jacobian.

        Returns
        -------
        J : float
            Jacobian.
        Jacob : np.ndarray
            Jacobian matrix.
        '''
        J = self._length / 2
        Jacob = np.array(J)
        return J, Jacob
    def getRotation(self)   -> float:
        '''
        Gets the element's rotation angle in radians.

        Returns
        -------
        rot_x : float
            Rotation angle in radians.
        '''
        return self._rot_x

    ## --------------------- Setters --------------------- ##
    def setRotation(self, angle):
        '''
        Rotates the element around its local x-axis.

        Parameters
        ----------
        angle : float
            The rotation angle in radians.
        '''
        C = np.cos(angle)
        S = np.sin(angle)
        # Rotation about local 'x' axis
        Rx = np.array([[1,0,0],[0,C,S],[0,-S,C]])
        # Updates the local axes
        self._localAxes = Rx @ self._localAxes
        self._rot_x     = angle

## ---------------------------------------------------------------- ##
##                           2D Elements                            ##
## ---------------------------------------------------------------- ##
class Elem2D(Elem):
    '''
    General class for 2D finite elements.\n
    
    More specific elements like 2D sections, plane strain or plane stress 
    extend this class with more specific methods for their stiffness matrices. 
    '''
    __elementCount = 0
    ## --------------------- Constructor --------------------- ##
    def __init__(self, nodes:tuple[Node,...], index=None):
        # If no index is given, use the class variable 'elementCount'
        if index is None:
            index = __elementCount
            __elementCount += 1
        ## Inheritance from the element class
        Elem.__init__(self, nodes, index)

        self.nNodes = len(nodes)

        # Element type
        self.elemType    = ''
        if self.nNodes in [3,6]:
            self.elemType  = 'Tri'
        elif self.nNodes in [4,8,9]:
            self.elemType  = 'Quad' 
        else:
            raise Exception(f'Element type not implemented: {self.elemType}')
    

    ## JACOBIAN AND PHYSICAL DERIVATIVES
    def getJacob(self, zetas:np.ndarray) -> list[float, np.ndarray, np.ndarray]:
        '''
        Calculates the Jacobian and physical derivatives for 2D elements.

        Parameters
        ----------
        zetas : np.ndarray
            Natural position of the point inside the element (3d-array for Tri, 2d-array for Quad).

        Returns
        -------
        J : float
            Jacobian.
        dNdy : np.ndarray
            Derivatives of shape functions with respect to y.
        dNdz : np.ndarray
            Derivatives of shape functions with respect to z.
        '''
        # Zetas (Necessary for higher-order elements)
        node_coords = np.array(self.getNodeCoords())
        ys  = node_coords[:,1]
        zs  = node_coords[:,2]

        match self.nNodes:
            ## T3 ___________________________________________
            case 3:
                # Natural coordinates of the point
                zeta1, zeta2, zeta3 = zetas

                # Coordinates of the element nodes
                y1,y2,y3 = ys
                z1,z2,z3 = zs

                Jy1 = y1
                Jy2 = y2
                Jy3 = y3

                Jz1 = z1
                Jz2 = z2
                Jz3 = z3

                # Jacobian matrix
                jacobianMatrix = np.array([[  1,  1,  1],
                                           [Jy1,Jy2,Jy3],
                                           [Jz1,Jz2,Jz3]])

                # Jacobian
                J = .5*np.linalg.det(jacobianMatrix) # Fellipa p. 364

                # Vector of partial derivatives
                derivs = 1/(2*J)*np.array([[z2-z3, z3-z1, z1-z2], 
                                           [y3-y2, y1-y3, y2-y1]])
                
                dNdy = derivs[0,:]
                dNdz = derivs[1,:]

            ## T6 __________________________________________
            case 6:
                # Natural coordinates of the point
                zeta1, zeta2, zeta3 = zetas
                
                # Coordinates of the element nodes
                y1,y2,y3,y4,y5,y6 = ys
                z1,z2,z3,z4,z5,z6 = zs
                
                # Derivatives
                dNdzeta1 = np.array([4*zeta1-1,      0,      0, 4*zeta2,      0, 4*zeta3])
                dNdzeta2 = np.array([      0, 4*zeta2-1,      0, 4*zeta1, 4*zeta3,      0])
                dNdzeta3 = np.array([      0,      0, 4*zeta3-1,      0, 4*zeta2, 4*zeta1])

                # Components of the Jacobian matrix
                Jy1 = np.dot(ys,dNdzeta1)
                Jy2 = np.dot(ys,dNdzeta2)
                Jy3 = np.dot(ys,dNdzeta3)

                Jz1 = np.dot(zs,dNdzeta1)
                Jz2 = np.dot(zs,dNdzeta2)
                Jz3 = np.dot(zs,dNdzeta3)

                # Jacobian matrix
                jacobianMatrix = np.array([[ 1, 1, 1],
                                           [Jy1,Jy2,Jy3],
                                           [Jz1,Jz2,Jz3]])
                
                # Jacobian
                J = .5*np.linalg.det(jacobianMatrix)

                # Vector of partial derivatives
                P = 1/(2*J)*np.array([[Jz2-Jz3, Jy3-Jy2],
                                      [Jz3-Jz1, Jy1-Jy3],
                                      [Jz1-Jz2, Jy2-Jy1]])
                derivs = P.T@np.array([dNdzeta1, 
                                       dNdzeta2,
                                       dNdzeta3])
                
                dNdy = derivs[0,:]
                dNdz = derivs[1,:]

            ## Q4 __________________________________________
            case 4:
                # Natural coordinates of the point
                eta, zeta = zetas
                
                # Derivatives
                dNdeta  = 1/4*np.array([zeta-1,  1-zeta, 1+zeta, -(1+zeta)])
                dNdzeta = 1/4*np.array([ eta-1, -(1+eta),  1+eta,     1-eta])

                # Components of the Jacobian matrix
                Jy_eta  = np.dot(ys,dNdeta)
                Jy_zeta = np.dot(ys,dNdzeta)

                Jz_eta  = np.dot(zs,dNdeta)
                Jz_zeta = np.dot(zs,dNdzeta)

                # Jacobian matrix
                jacobianMatrix = np.array([[ Jy_eta, Jz_eta],   # Component Jyi multiplies the d_eta derivative
                                           [Jy_zeta,Jz_zeta]]) # Component Jzi multiplies the d_zeta derivative

                # Jacobian
                J = np.linalg.det(jacobianMatrix)

                # Inverse Jacobian
                jacobInv = np.linalg.inv(jacobianMatrix)

                # Physical Derivatives
                dNdksi = np.vstack((dNdeta, dNdzeta))#.reshape(2,-1)
                dNdy = jacobInv[0,:]@dNdksi
                dNdz = jacobInv[1,:]@dNdksi

            ## Q8 __________________________________________
            case 8:
                # Natural coordinates of the point
                eta, zeta = zetas
                
                # Derivatives
                dNdeta  = 1/4*np.array([ zeta*(2*eta-1)*(zeta-1),  zeta*(2*eta+1)*(zeta-1),  zeta*(2*eta+1)*(zeta+1),  zeta*(2*eta-1)*(zeta+1),
                                             4*eta*zeta*(1-zeta), -2*(2*eta+1)*(zeta**2-1),     -4*eta*zeta*(zeta+1),  2*(1-2*eta)*(zeta**2-1)])
                dNdzeta = 1/4*np.array([ eta*(eta-1)*(2*zeta-1),    eta*(eta+1)*(2*zeta-1),    eta*(eta+1)*(2*zeta+1),  eta*(eta-1)*(2*zeta+1),
                                          2*(1-2*zeta)*(eta**2-1),      -4*eta*zeta*(eta+1), -2*(eta**2-1)*(2*zeta+1),      4*eta*zeta*(1-eta)])

                # Components of the Jacobian matrix
                Jy_eta  = np.dot(ys,dNdeta)
                Jy_zeta = np.dot(ys,dNdzeta)

                Jz_eta  = np.dot(zs,dNdeta)
                Jz_zeta = np.dot(zs,dNdzeta)

                # Jacobian Matrix
                jacobianMatrix = np.array([[ Jy_eta, Jz_eta],   # Component Jyi multiplies the d_eta derivative
                                           [Jy_zeta,Jz_zeta]]) # Component Jzi multiplies the d_zeta derivative
                
                # Jacobian
                J = np.linalg.det(jacobianMatrix)
                
                # Inverse Jacobian
                jacobInv = np.linalg.inv(jacobianMatrix)

                # Physical Derivatives
                dNdksi = np.vstack((dNdeta, dNdzeta))#.reshape(2,-1)
                dNdy = jacobInv[0,:]@dNdksi
                dNdz = jacobInv[1,:]@dNdksi
            
            ## Q9 __________________________________________
            case 9:
                # Natural coordinates of the point
                eta, zeta = zetas
                
                # Derivatives
                dNdeta  = 1/4*np.array([ zeta*(2*eta-1)*(zeta-1),  zeta*(2*eta+1)*(zeta-1),  zeta*(2*eta+1)*(zeta+1),  zeta*(2*eta-1)*(zeta+1),
                                             4*eta*zeta*(1-zeta), -2*(2*eta+1)*(zeta**2-1),     -4*eta*zeta*(zeta+1),  2*(1-2*eta)*(zeta**2-1),
                                             8*eta*(zeta**2-1)])
                dNdzeta = 1/4*np.array([ eta*(eta-1)*(2*zeta-1),    eta*(eta+1)*(2*zeta-1),    eta*(eta+1)*(2*zeta+1),  eta*(eta-1)*(2*zeta+1),
                                          2*(1-2*zeta)*(eta**2-1),      -4*eta*zeta*(eta+1), -2*(eta**2-1)*(2*zeta+1),      4*eta*zeta*(1-eta),
                                             8*zeta*(eta**2-1)])

                # Components of the Jacobian matrix
                Jy_eta  = np.dot(ys,dNdeta)
                Jy_zeta = np.dot(ys,dNdzeta)

                Jz_eta  = np.dot(zs,dNdeta)
                Jz_zeta = np.dot(zs,dNdzeta)

                # Jacobian Matrix
                jacobianMatrix = np.array([[ Jy_eta, Jz_eta],   # Component Jyi multiplies the d_eta derivative
                                           [Jy_zeta,Jz_zeta]]) # Component Jzi multiplies the d_zeta derivative
                
                # Jacobian
                J = np.linalg.det(jacobianMatrix)
                
                # Inverse Jacobian
                jacobInv = np.linalg.inv(jacobianMatrix)

                # Physical Derivatives
                dNdksi = np.vstack((dNdeta, dNdzeta))#.reshape(2,-1)
                dNdy = jacobInv[0,:]@dNdksi
                dNdz = jacobInv[1,:]@dNdksi

            case _:
                raise Exception(f'Invalid number of nodes: {self.nodes}')
            
        return J, dNdy.reshape(1,-1), dNdz.reshape(1,-1)
    ## SHAPE FUNCTIONS
    def getShapeFun(self, zetas:np.ndarray) -> np.ndarray:
        '''
        Determines the value of the shape functions at a point in a 2D element.

        Parameters
        ----------
        zetas : np.ndarray
            Natural position of the point inside the element (3d-array for Tri, 2d-array for Quad).

        Returns
        -------
        N : np.ndarray
            Vector with the value of the shape functions evaluated at 'zetas'.
        '''
        match self.nNodes:
            ## T3 - LINEAR ------------------------------------------------------ ##
            case 3:
                # Natural coordinates of the point
                zeta1,zeta2,zeta3 = zetas
                N1 = zeta1
                N2 = zeta2
                N3 = zeta3
                N  = np.array([N1,N2,N3])

            ## T6 - BILINEAR ---------------------------------------------------- ##
            case 6:
                # Natural coordinates of the point
                zeta1,zeta2,zeta3 = zetas
                N1 = zeta1*(2*zeta1-1)
                N2 = zeta2*(2*zeta2-1)
                N3 = zeta3*(2*zeta3-1)
                N4 = 4*zeta1*zeta2
                N5 = 4*zeta2*zeta3
                N6 = 4*zeta3*zeta1
                N  = np.array([N1,N2,N3,N4,N5,N6])

            ## Q4 - BILINEAR ---------------------------------------------------- ##
            case 4:
                # Natural coordinates of the point
                eta,zeta = zetas
                N1 = 1/4*(1-eta)*(1-zeta)
                N2 = 1/4*(1+eta)*(1-zeta)
                N3 = 1/4*(1+eta)*(1+zeta)
                N4 = 1/4*(1-eta)*(1+zeta)
                N  = np.array([N1,N2,N3,N4])

            ## Q8 - QUADRATIC -------------------------------------------------- ##
            case 8:
                # Natural coordinates of the point
                eta,zeta = zetas
                N1 = 1/4*eta*zeta*(eta-1)*(zeta-1)
                N2 = 1/4*eta*zeta*(eta+1)*(zeta-1)
                N3 = 1/4*eta*zeta*(eta+1)*(zeta+1)
                N4 = 1/4*eta*zeta*(eta-1)*(zeta+1)
                N5 = 1/2*zeta*(1-eta**2)*(  zeta-1)
                N6 = 1/2* eta*(  eta+1)*(1-zeta**2)
                N7 = 1/2*zeta*(1-eta**2)*(  zeta+1)
                N8 = 1/2* eta*(  eta-1)*(1-zeta**2)
                N  = np.array([N1,N2,N3,N4,N5,N6,N7,N8])

            ## Q9 - BIQUADRATIC ------------------------------------------------ ##
            case 9:
                # Natural coordinates of the point
                eta,zeta = zetas
                N1 = 1/4*eta*zeta*(eta-1)*(zeta-1)
                N2 = 1/4*eta*zeta*(eta+1)*(zeta-1)
                N3 = 1/4*eta*zeta*(eta+1)*(zeta+1)
                N4 = 1/4*eta*zeta*(eta-1)*(zeta+1)
                N5 = 1/2*zeta*(1-eta**2)*(  zeta-1)
                N6 = 1/2* eta*(  eta+1)*(1-zeta**2)
                N7 = 1/2*zeta*(1-eta**2)*(  zeta+1)
                N8 = 1/2* eta*(  eta-1)*(1-zeta**2)
                N9 = (1-eta**2)*(1-zeta**2)
                N  = np.array([N1,N2,N3,N4,N5,N6,N7,N8,N9])

            case _:
                raise Exception(f'Invalid number of nodes {self.nNodes}')

        return N.reshape(1,-1)

    ## ------------------------------------------------------- ##
    ##                   QUADRATURE FUNCTIONS                  ##
    ## ------------------------------------------------------- ##
    # Triangular elements
    def __quadTri(self, degree):
        '''
        Returns the points (in triangular coordinates) and weights for quadrature on triangles.

        Parameters
        ----------
        degree : int
            Quadrature degree (1 to 5).

        Returns
        -------
        pnts : np.ndarray (n, 3)
            The quadrature points in triangular coordinates.
        weights : np.ndarray (1, n)
            The corresponding weights.
        '''

        # Internal Points
        match degree:
            case 1:
                pnts = [[1/3,1/3,1/3]]
                weights = [1]
            case 2:
                w  = 1/3
                p1 = 2/3; p2 = 1/6
                pnts = [[p1,p2,p2],
                        [p2,p1,p2],
                        [p2,p2,p1]]
                weights = [w,w,w]
            case 3 | 4:
                # Quadrature of degree 4, there is no stable quadrature of degree 3
                p1  = (8 - np.sqrt(10) + np.sqrt(38 - 44*np.sqrt(2/5)))/18
                p11 = 1 - 2*p1
                p2  = (8 - np.sqrt(10) - np.sqrt(38 - 44*np.sqrt(2/5)))/18
                p22 = 1 - 2*p2
                w1 = (620 + np.sqrt(213125 - 53320*np.sqrt(10)))/3720
                w2 = (620 - np.sqrt(213125 - 53320*np.sqrt(10)))/3720
                pnts = [[p11,  p1,  p1],
                        [ p1, p11,  p1],
                        [ p1,  p1, p11],
                        [p22,  p2,  p2],
                        [ p2, p22,  p2],
                        [ p2,  p2, p22]]
                weights = [w1,w1,w1,w2,w2,w2]
            case 5:
                p1  = 1/3

                p2  = (6 - np.sqrt(15))/21
                p22 = (9 + 2*np.sqrt(15))/21

                p3  = (6 + np.sqrt(15))/21
                p33 = (9 - 2*np.sqrt(15))/21

                w1  = 9/40
                w2  = (155 - np.sqrt(15))/1200
                w3  = (155 + np.sqrt(15))/1200

                pnts = [[  p1,  p1,  p1],
                        [ p22,  p2,  p2],
                        [  p2, p22,  p2],
                        [  p2,  p2, p22],
                        [ p33,  p3,  p3],
                        [  p3, p33,  p3],
                        [  p3,  p3, p33]]
                weights = [w1,w2,w2,w2,w3,w3,w3]
            case _:
                raise NotImplementedError(f'{degree} is not a valid or implemented degree')
        return np.array(pnts), np.array(weights)
    # Quadrilateral elements
    def __quadQuad(self, degree):
        '''
        Returns the points (in natural coordinates) and weights for 1D quadrature.

        Parameters
        ----------
        degree : int
            Quadrature degree.

        Returns
        -------
        pnts : np.ndarray (1, n)
            The quadrature points in natural coordinates.
        weights : np.ndarray (1, n)
            The corresponding weights.
        '''
        # Necessary Number of Points
        n_gp = np.ceil((degree+1)/2)

        # Internal Points
        # Data from Fish (2007)
        match n_gp:
            case 1:
                pnts = [0]
                weights = [2]
            case 2:
                p  = 1/np.sqrt(3)
                
                w  = 1

                pnts  = [-p, p]
                weights = [w,w]
            case 3 :
                p1  = 0.7745966692
                p2  = 0

                w1  = 0.5555555556
                w2  = 0.8888888889
                
                pnts  = [-p1, p2, p1]
                weights = [ w1, w2, w1]
            case 4:
                p1  = 0.8611363116
                p2  = 0.3399810436

                w1  = 0.3478548451
                w2  = 0.6521451549

                pnts  = [-p1, -p2, p2, p1]
                weights = [ w1,  w2 ,w2, w1]
            case 5:
                p1  = 0.9061798459
                p2  = 0.5384693101
                p3  = 0

                w1  = 0.2369268851
                w2  = 0.4786286705
                w3  = 0.5688888889

                pnts  = [-p1, -p2, p3, p2, p1]
                weights = [ w1,  w2, w3, w2, w1]
            case 6:
                p1  = 0.9324695142
                p2  = 0.6612093865
                p3  = 0.2386191861

                w1  = 0.1713244924
                w2  = 0.3607615730
                w3  = 0.4679139346

                pnts  = [-p1, -p2, -p3, p3, p2, p1]
                weights = [ w1,  w2,  w3, w3, w2, w1]
            case _:
                raise NotImplementedError(f'{degree} is not a valid or implemented degree')
        return np.array(pnts), np.array(weights)
    
    ## QUADRATURE 
    # Choose quadrature according to element type
    def getQuadrature(self, degree):
        '''
        Selects the appropriate Gaussian quadrature method based on the element type.

        Parameters
        ----------
        degree : int
            The degree of the quadrature to use.

        Returns
        -------
        pnts : np.ndarray
            The quadrature points (Gaussian points).
        weights : np.ndarray
            The corresponding weights.
        '''
        # Element type
        if self.elemType == 'Tri':
            return self.__quadTri(degree)
        elif self.elemType == 'Quad':
            return self.__quadQuad(degree)
        else:
            raise Exception('Element type not implemented')

    ## ------------------------------------------------------- ##
    ##                  FORMULATION-SPECIFIC                   ##
    ##                   ABSTRACT FUNCTIONS                    ##
    ## ------------------------------------------------------- ##
    # Specific elements MUST define the methods below
    ## LOCAL STIFFNESS MATRIX 
    @abstractmethod
    def getKlocal(self,degree:int) -> np.ndarray:
        '''
        Calculates the local stiffness matrix of the element.

        Parameters
        ----------
        degree : int
            Quadrature degree for numerical integration.

        Returns
        -------
        sum_val : np.ndarray
            Local stiffness matrix.
        '''
        pass
    
    ## LOCAL FORCE VECTOR
    @abstractmethod
    def getFlocal(self,degree:int) -> np.ndarray: 
        '''
        Calculates the local force vector of the element.

        Parameters
        ----------
        degree : int
            Quadrature degree for numerical integration.

        Returns
        -------
        sum_val : np.ndarray
            Local force vector.
        '''
        pass

## ------------------------------------------------------- ##
##                   2D Section Element                    ##
## ------------------------------------------------------- ##
class SectionElem(Elem2D):
    '''
    Class for 2D section elements. \n
    
    This element was designed for solving the Saint-Venant torsion 
    problem for the warping of a general multi-material cross-section. \n
    
    This class defines the element's local stiffness matrix and force vector 
    to be used inside a global assembly routine.  
    '''
    def __init__(self, nodes:tuple[Node,...], index=None):
        ## Inheritance from the element class
        Elem2D.__init__(self, nodes, index)
        
    ## LOCAL STIFFNESS MATRIX 
    def getKlocal(self,degree:int) -> np.ndarray:
        zetas, weights = self.getQuadrature(degree)

        sum_val = 0
        # Quadrature for the element; Iterates through Gaussian points
        if self.elemType == 'Tri':
            for i, coords in enumerate(zetas):
                #      i: index
                # coord: natural coordinates of the point (zeta1,2,3)

                J, By, Bz = self.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point

                # Integrand Term
                BB = By.T @ By + Bz.T @ Bz
                #BB = By[...,None] @ By + Bz[...,None] @ Bz


                sum_val += weights[i]*J*BB
        elif self.elemType == 'Quad':
            for i, zeta_i in enumerate(zetas):
                for j, zeta_j in enumerate(zetas):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)
                    coords = [zeta_i, zeta_j]
                    J, By, Bz = self.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point

                    # Integrand Term
                    BB = By.T @ By + Bz.T @ Bz
                    #BB = By[...,None] @ By + Bz[...,None] @ Bz


                    sum_val += weights[i]*weights[j]*J*BB

        return sum_val


    ## LOCAL FORCE VECTOR
    def getFlocal(self,degree:int) -> np.ndarray:         
        points, weights = self.getQuadrature(degree)
        node_coords = np.array(self.getNodeCoords())
        ys  = node_coords[:,1]
        zs  = node_coords[:,2]

        sum_val = 0
        # Quadrature for the element; Iterates through Gaussian points

        if self.elemType == 'Tri':
            for i, coords in enumerate(points):
                #      i: index
                # coord: natural coordinates of the point (zeta1,2,3)

                # Integrand Term
                N           = self.getShapeFun(coords)
                J, By, Bz   = self.getJacob(coords) 
                F           = By.T @ N @ zs.reshape(-1,1) - Bz.T @ N @ ys.reshape(-1,1)

                sum_val += weights[i]*J*F
        elif self.elemType == 'Quad':
            for j, zeta_j in enumerate(points):
                for i, zeta_i in enumerate(points):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)
                    coords = [zeta_i, zeta_j]
                    J, By, Bz = self.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point

                    # Integrand Term
                    N = self.getShapeFun(coords)
                    F = By.T @ N @ zs.reshape(-1,1) - Bz.T @ N @ ys.reshape(-1,1)

                    sum_val += weights[i]*weights[j]*J*F

        return sum_val




## ================================================================ ##
##                            MESH CLASSES                          ##
## ---------------------------------------------------------------- ##
class Mesh:
    '''
    Generic class for meshes.
    '''
    def __init__(self, elemType: Elem1D | Elem2D, 
                 elems:list[Elem]=None, 
                 coordinates:list[list[float,float,float]]=None, 
                 connectivity:list[list[int,int]]=None):
        
        ## ---------- Element Initialization ---------- ##
        # If only a list of elements is provided
        if elems is not None:
            self._elements = list(elems)

            elementsNodes = []
            coords = []
            for elem in elems:
                elementsNodes.append([node._index for node in elem._nodes]) # Connectivity
            self._elementsNodes = elementsNodes
            # print(self._elementsNodes)
            
            # Ordered list of nodes
            elemNodes = np.array(elementsNodes)
            elemNodes = set(elemNodes.flatten())
            nodeIndexes = list(elemNodes)
            indices   = []
            nodes     = []
            for i in nodeIndexes:
                for elem in elems:
                    for node in elem._nodes:
                        if node._index == i and node._index not in indices:
                            nodes.append(node)
                            coords.append(node._coords)
                            indices.append(node._index)

            self._nodeCoords = np.array(coords)
            self._nodes      = nodes
            # print(self._nodeCoords)

        # If lists of coordinates and connectivity are provided
        elif (coordinates is not None) and (connectivity is not None):
            # Creation of the list of node instances from coordinate list
            nodes   = [Node(nodeCoords,indNode) for indNode,nodeCoords in enumerate(coordinates)]

            # Elements list
            elems = []
                        
            # Iterating through the element vector to create element objects
            for e,elemNodes in enumerate(connectivity):
                elem_nodes = []
                mesh_nodes = []
                
                # Iterate through element node list
                for node_index in elemNodes:
                    # Searching for the element's node instance and adding it to a list
                    elem_nodes.append(nodes[node_index])
                    # Creating a list of mesh nodes
                    if nodes[node_index] not in mesh_nodes:
                        mesh_nodes.append(nodes[node_index])
                        
                # Passing the list of node instances to an element
                if elemType == Elem1D:
                    elems.append(Elem1D(elem_nodes, index=e)) # Added index to Elem1D constructor
                elif elemType == Elem2D:
                    elems.append(Elem2D(elem_nodes, index=e)) # Elem2D constructor doesn't take index
                elif elemType == SectionElem:
                    elems.append(SectionElem(elem_nodes, index=e)) # Elem2D constructor doesn't take index
            
            # Save lists to object    
            self._elements  = elems
            self._nodes     = mesh_nodes
            
            self._elementsNodes = connectivity
            
            # Check if the coordinates are 2D
            if np.shape(coordinates)[1] == 2:
                # Create an column array of zeros
                nElems = np.shape(coordinates)[0]
                zerosX = np.zeros(nElems).reshape(-1,1)
                # Add to coordinates array
                coordinates = np.hstack([zerosX, coordinates])
            self._nodeCoords    = coordinates
            
            ## Generate degree of freedom and valid nodes list
            dofs = []
            validNodes = []
            # validCoords = []
            
            if elemType == SectionElem:
                # Get a sorted list of nodes related to elements
                # The node's DoF will be the node number index on the sorted list
                sorted_nodes = np.sort(np.unique(connectivity),axis=None)
                # Number of Dofs
                nDofs = len(sorted_nodes)
                
                dofCount = 0
                # Iterate through coordinates list
                for n,coords in enumerate(coordinates):
                    # If the node is part of the mesh, add its index to the dof list 
                    # and increment dof count
                    if n in sorted_nodes:
                        dofs.append(dofCount)
                        # Add coordinates to valid mesh nodes
                        # validCoords.append(coords)
                        validNodes.append(n)
                        dofCount += 1
                    # if not, add a nan to the dof list
                    else:
                        dofs.append(np.nan)
                    
            elif elemType == Elem2D:
                raise NotImplementedError(f'Element type {elemType} not yet implemented.')
            else:
                raise Exception('Mesh element type invalid')        
            # Store the dofs list
            self._dofs = dofs
            # Store the number of dofs
            self._nDofs = nDofs
            # Store valid nodes (nodes that belong to at least one element)
            self._validNodes = np.array(validNodes)
            # Store coordinates of valid nodes
            # self._validCoords = np.array(validCoords)
                        
            
                
        
        # Raise error if signature doesn't follow either convention
        else:
            raise Exception('''Element list and coordinates may be provided at the same time. 
                            \n * If the element list exists, use Mesh(elems=elementList). 
                            \n * If the coordinates list and connectivity list exist, use Mesh(coordinates=coordList,connectivity=connectList)  ''')
        
        
    def __repr__(self):
        lista = [f'Element: {elem._index}, {elem}' for elem in self._elements]
        return str(lista)

    ## --------------------- Getters --------------------- ##
    def getElems(self):
        """
        Gets the list of elements in the mesh.

        Returns
        -------
        _elements : list
            List of element instances.
        """
        return self._elements
    def getElemIndex(self, elem:Elem):
        '''
        Get an elements index.

        Parameters
        ----------
        elem : Elem
            The element instance to search for.
        '''
        return self._elements.index(elem)
    def getElem(self, index:int):
        '''
        Get element from index.

        Parameters
        ----------
        index : int
            The element's index.
        '''
        return self._elements[index]
    def getElemDof(self, elem):
        '''
        Get element degrees of freedom.
        
        Parameters
        ----------
        elem : Elem
            The element instance to search for.

        Returns
        -------
        elementDof : np.Ndarray of int
            Array of element degrees of freedom.
        '''
        # Get element node indices
        nodesIndices = np.array(elem.getNodeIndices(),dtype=int)
        # Get element DoFs
        elementDof   = np.array(self._dofs)[nodesIndices]
        # Cast as an array of ints
        elementDof   = np.array(elementDof,dtype=int)
        return elementDof

    ## --------------------- Methods --------------------- ##
    def add(self, elem:Elem):
        '''
        Add an element to the mesh.

        Parameters
        ----------
        elem : Elem
            The element instance to add.
        '''
        self._elements.append(elem)
    def delete(self, index:int):
        '''
        Remove an element from the mesh by its position index in the internal list.

        Parameters
        ----------
        index : int
            The index of the element to remove.
        '''
        self._elements.remove(self._elements[index])
    def dataframe(self, precision:int = 4):
        '''
        Shows the mesh information as a pandas DataFrame.

        Parameters
        ----------
        precision : int, optional
            The number of decimal places for coordinate/rotation values. The default is 4.

        Returns
        -------
        df : pandas.DataFrame
            DataFrame containing element data.
        '''
        import pandas as pd

        col0 = 'Element'
        # Key for node indices
        col1 = 'Nodes'
        # Key for coordinates
        col2 = 'Coordinates'

        dic  = {col0:[], col1:[], col2:[]}
        dfIndex = []
        for elem in self._elements: 
            elemIndex = elem._index
            nodes     = elem._nodes
            nodeIndices = [node._index  for node in nodes]
            coords    = [np.round(node._coords, precision) for node in nodes]
            rot       = np.round(elem._rot_x, precision)
            i,j,k     = np.round(elem.getLocalAxes(),3)
            # class_name  = elem.__str__()
            
            # Building dictionary
            dic[col0].append(elemIndex)
            dic[col1].append(nodeIndices)
            dic[col2].append(coords)
            # dic[col4].append(class_name)
            dfIndex.append('')             # to inhibit the index column
        df = pd.DataFrame(data=dic,index=dfIndex)
        # df.index.name = ''
        return df
    

## ---------------------------------------------------------------- ##
##                              1D Mesh                             ##
## ---------------------------------------------------------------- ##
class Mesh1D(Mesh):
    '''
    Class for reticulated meshes, i.e., of rectilinear elements in space.
    Inherits from the generic Mesh class.
    '''
    def __init__(self, elems:list[Elem]=None, coordinates:list[list[float,float,float]]=None, connectivity:list[list[int,int]]=None):
        ## Inheritance from the mesh class
        Mesh.__init__(self,Elem1D,elems,coordinates,connectivity)            
    
    def dataframe(self, precision:int = 4):
        '''
        Shows the mesh information as a pandas DataFrame.

        Parameters
        ----------
        precision : int, optional
            The number of decimal places for coordinate/rotation values. The default is 4.

        Returns
        -------
        df : pandas.DataFrame
            DataFrame containing element data.
        '''
        import pandas as pd

        col0 = 'Element'
        # Key for node indices
        col1 = 'Nodes'
        # Key for coordinates
        col2 = 'Coordinates'
        col21 = 'Length'
        col3 = 'Rot. about i_e [rad]'
        col4 = 'i_e'
        col5 = 'j_e'
        col6 = 'k_e'
        dic  = {col0:[], col1:[], col2:[], col21:[], col3:[], col4:[], col5:[], col6:[]}
        dfIndex = []
        for elem in self._elements: 
            elemIndex = elem._index
            nodes     = elem._nodes
            nodeIndices = [node._index  for node in nodes]
            coords    = [np.round(node._coords, precision) for node in nodes]
            rot       = np.round(elem._rot_x, precision)
            i,j,k     = np.round(elem.getLocalAxes(),precision)
            length    = np.round(elem.getLength(),precision)
            # class_name  = elem.__str__()
            
            # Building dictionary
            dic[col0].append(elemIndex)
            dic[col1].append(nodeIndices)
            dic[col2].append(coords)
            dic[col21].append(length)
            dic[col3].append(rot)
            dic[col4].append(i)
            dic[col5].append(j)
            dic[col6].append(k)
            # dic[col4].append(class_name)
            dfIndex.append('')             # to inhibit the index column
        df = pd.DataFrame(data=dic,index=dfIndex)
        # df.index.name = ''
        return df
    
    ## --------------------- Mesh Plotting --------------------- ##
    def plotMesh(self, 
                 showElemLabels:bool=True, 
                 showNodeLabels:bool=True, 
                 showLocalAxes:bool=True, 
                 showSections:bool=False, 
                 sectionScaleFactor:float=1,
                 localAxisScaleFactor:float=0.25,
                 **kwargs):
        '''
        Plot the reticulated mesh.

        Parameters
        ----------
        showElemLabels : bool, optional
            Show element index labels. The default is True.
        showNodeLabels : bool, optional
            Show node index labels. The default is True.
        showLocalAxes : bool, optional
            Show the local coordinate axes for each element. The default is True.
        showSections : bool, optional
            Show the cross-sections of the elements (requires SectionMesh). The default is False.
        sectionScaleFactor : float, optional
            Scaling factor for plotting cross-sections. The default is 1.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The 3D axes object with the plotted mesh.
        '''
        import matplotlib.pyplot as plt
        # Create figure
        fig = plt.figure(**kwargs)
        
        # Create axes
        elevation    = 30
        azimuth      = 45
        ax = fig.add_subplot(1, 1, 1, projection='3d',azim=azimuth, elev=elevation)
        # ax           = plt.axes(projection='3d',azim=azimuth, elev=elevation)
        
        # maxX, maxY, maxZ = -np.inf,-np.inf,-np.inf
        # minX, minY, minZ = np.inf,np.inf,np.inf

        # Loop for elements
        for elemIndex, elem in enumerate(self._elements): 
            ## ----------------------------------------------------------------- ##
            ##                           SHOW ELEMENTS                           ##
            ## ----------------------------------------------------------------- ##
            # Getting node coordinates and converting to numpy array
            nodes_coords = np.array( elem.getNodeCoords() ) 
            # nodes_coords = list([X1,Y1,Z1],[X2,Y2,Z3])
            X = nodes_coords[:,0]
            Y = nodes_coords[:,1]
            Z = nodes_coords[:,2]
            
            # # Calculate current maximum an minimum coordinates
            # maxX = np.max(np.append(X,maxX)) 
            # maxY = np.max(np.append(Y,maxY)) 
            # maxZ = np.max(np.append(Z,maxZ)) 
            
            # minX = np.min(np.append(X,minX)) 
            # minY = np.min(np.append(Y,minY)) 
            # minZ = np.min(np.append(Z,minZ)) 
            
            # Plot elements
            ax.plot3D(X,Y,Z,marker='o',lw=1,c='k',mfc='w',markersize=5,label=None)

            ## ----------------------------------------------------------------- ##
            ##                          SHOW ELEMENT LABELS                      ##
            ## ----------------------------------------------------------------- ##
            # Get the midpoint
            [Xm, Ym, Zm] = elem._midPoint

            # Element Labels
            if showElemLabels:
                label = '  ('+str(elemIndex)+')'    
                ax.text(Xm,Ym,Zm,label,zdir=None,c='k',va='top',ha='left',fontweight='bold',fontstyle='italic')
            
            # Node Labels
            if showNodeLabels:
                # Searching for nodes of each element 
                nodes     = elem._nodes
                nodeIndices = [node._index  for node in nodes]    # node indices
                coords    = [node._coords for node in nodes]    # node coordinates
                # Searching for the index and coordinates of each node
                for i,index in enumerate(nodeIndices):
                    label = str(index) + '  '              # label for the node
                    X,Y,Z = coords[i]                      # node position
                    ax.text(X,Y,Z,label,zdir=None,c='#04a301',va='top',ha='right',fontweight='bold')
            
            ## ----------------------------------------------------------------- ##
            ##                           SHOW LOCAL AXES                           ##
            ## ----------------------------------------------------------------- ##
            if showLocalAxes:
                ### Loop for each direction cosine
                for i in range(3):
                    # cDir = self.cosDir[nElem][i] #+ Xm
                    # pm   = [Xm,Ym,Zm] # Element midpoint
                    cDir = elem._localAxes[i]
                    pm   = elem._midPoint
                    p2   = pm + localAxisScaleFactor*cDir
                    
                    # Local axes
                    x = [pm[0], p2[0]] 
                    y = [pm[1], p2[1]]
                    z = [pm[2], p2[2]]
                    
                    # Determine the line color
                    colors = ['#ff7f7f' ,'#A6D400', '#7fb0ff']
                    # colors = ['#ff7f7f' ,'#f7f456', '#7fb0ff']
                    # colors = ['#ff7f7f' ,'#fffc7f', '#7fb0ff']
                    # colors = ['#ff7f7f' ,'#f87fff', '#7fb0ff']
                    # colors = ['#ff7f7f' ,'#fff47f', '#7fb0ff']
                    # colors = ['#ff8150', '#e59400' , '#1073d6']
                        
                    # Plot lines indicating local axes
                    ax.plot3D(x,y,z,marker=None,lw=2,c=colors[i])
            
            ## ----------------------------------------------------------------- ##
            ##                            SHOW SECTIONS                            ##
            ## ----------------------------------------------------------------- ##
            if showSections:
                y_cg        = elem.secao.propsArea['Y_CG']   # Centroid position
                z_cg        = elem.secao.propsArea['Z_CG']
                elemNodes   = elem.secao._elementsNodes      # Section mesh connectivity
                nNodes      = elem.secao.nNodes
                # Correcting the system origin to the centroid and element midpoint
                nodes2D             = elem.secao._nodeCoords
                section_node_coords = np.zeros(nodes2D.shape)    # Section mesh coordinates

                # Correcting the system origin to the centroid and applying the scale factor
                section_node_coords[:,1] = sectionScaleFactor*(nodes2D[:,1] - y_cg)  
                section_node_coords[:,2] = sectionScaleFactor*(nodes2D[:,2] - z_cg)
                
                # print('elemNodes: \n', elemNodes)
                # Determining the vertex order
                match nNodes:
                    case 9:
                        indices = [0,4,1,5,2,6,3,7]
                    case 6:
                        indices = [0,3,1,4,2,5]
                    case 4:
                        indices = [0,1,2,3]
                    case 3:
                        indices = [0,1,2]
                # Direction cosine matrix of the 1D element
                Ae          = elem.getLocalAxes()            
                for elem2D_conn in elemNodes:                        # Iterate through 2D elements in the cross-section
                    # print('elem2D: \n',elem2D_conn)
                    
                    if type(elem2D_conn) != np.ndarray :
                        elem2D_conn = np.array(elem2D_conn)
                    # Reordering nodes to maintain consistent geometry
                    elem2D_conn = elem2D_conn[indices]

                    for k, node in enumerate(elem2D_conn):          # Iterate through the element's nodes
                        # Indices of the 2D element
                        IND0 = elem2D_conn[k]
                        try:
                            IND1 = elem2D_conn[k+1]
                        except IndexError:
                            # If pos0 is the last one, the next will be the initial one
                            IND1 = elem2D_conn[0]

                        # Transformation of body axes
                        X0,Y0,Z0 = Ae.T@section_node_coords[IND0,:] # node 0
                        X1,Y1,Z1 = Ae.T@section_node_coords[IND1,:] # node 1

                        Xs = np.array([X0,X1])
                        Ys = np.array([Y0,Y1])
                        Zs = np.array([Z0,Z1])

                        # Determining the vertex order
                        match nNodes:
                            case 9:
                                inds = [0,4,1,5,2,6,3,7]
                            case 6:
                                inds = [0,3,1,4,2,5]
                            case 4:
                                inds = [0,1,2,3]
                            case 3:
                                inds = [0,1,2]
                        XX = Xs+Xm
                        YY = Ys+Ym
                        ZZ = Zs+Zm
                        ax.plot3D(XX,YY,ZZ,marker=None,ls='--',lw=.5,c='m',alpha=.75)


        ax.set_proj_type('ortho')
        ax.set_xlabel(r'$X$',fontsize=18,usetex=True);  # ax.set_xlim([np.min(X),np.max(X)])
        ax.set_ylabel(r'$Y$',fontsize=18,usetex=True);  # ax.set_ylim([np.min(Y),np.max(Y)])
        ax.set_zlabel(r'$Z$',fontsize=18,usetex=True);  # ax.set_zlim([np.min(Z),np.max(Z)])

        if showLocalAxes:
            ax.legend(['Mesh','x-local','y-local','z-local'])
        else:
            ax.legend(['Mesh'])
        
        
        # Make axis the same equal scale
        ax.set_box_aspect([1,1,1])
        # ax.set_aspect('equal')
        
        # Set limits to largest and smallest coordinates among all axis to preserve aspect
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        zlim = ax.get_zlim()
        
        mmax = np.max([xlim,ylim,zlim])
        mmin = np.min([xlim,ylim,zlim])
        ax.set_xlim(mmin,mmax)
        ax.set_ylim(mmin,mmax)
        ax.set_zlim(mmin,mmax)
        
        fig.tight_layout()

        return fig, ax
    


    