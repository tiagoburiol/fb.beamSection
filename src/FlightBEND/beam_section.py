# Get definitions from the geometry module
from .fem_geometry          import *

import numpy as np
import time
from itertools import product

# Plotting functions
import matplotlib.pyplot as plt
from matplotlib             import gridspec, tri
from matplotlib.patches     import Polygon
from matplotlib.axes._axes  import Axes

# Sparse matrix funtions  
from scipy.sparse           import lil_matrix, csr_matrix
from scipy.sparse.linalg    import spsolve

## ---------------------------------------------------------------- ##
##                             2D Section                           ##
## ---------------------------------------------------------------- ##
class BeamSection(Mesh, SectionElem):
    '''
    Class for creating a beam cross-section from a mesh of 2D elements. 
    Area properties such as modulus-weighted area, moments and
    product of inetia, as well as mass  properties (from density distribuition)
    are calculated by integration using gaussian quadrature.  
    
    The torsion formulation is based on Saint-Venant's theory of uniform torsion 
    and the warping displacement is solved by a 2D finite element routine and 
    from that the torsional constant and similar properties are calculated. 

    Parameters
    ----------
    elems : list of Elem2D objects, optional
        List of Elem2D objects.
    coordenadas : list of list of float, optional
        List of node coordinates.
    conectividade : list of list of int, optional
        List of element node connectivity.
    young : list of list of int, optional
        List of Young's moduli for the elements.
    nu : list of list of int, optional
        List of Poisson's ratios for the elements.
    rho : list of list of int, optional
        List of mass densities for the elements.
    intDegree: int, optional 
        Degree for quadrature integration, default = 4.
    noSimulation: bool, optional 
       Flag to not run the simulation, default = False.
    displayTimes: bool, optional 
        Print assembly, solving and integration times at the end, default = True.
    '''
    def __init__(self, 
                 elems:         list[Elem2D]          = None, 
                 coordinates:   list[list[float]]     = None, 
                 connectivity:  list[list[int]]       = None,
                 young:         list[list[int]]|float = None,
                 nu:            list[list[int]]|float = None,
                 rho:           list[list[int]]|float = None,
                 intDegree:     int                   = 4,
                 noSimulation: bool                   = False,
                 displayTimes:  bool                  = True,
                 ):

        ## Inheritance from the 2D mesh class
        Mesh.__init__(self,SectionElem,elems,coordinates,connectivity)
        
        # Material properties
        if np.shape(young) == ():
            self.E      = young*np.ones(np.shape(connectivity)[0])
        else:
            self.E      = young
            
        if np.shape(nu) == ():
            self.nu     = nu*np.ones(np.shape(connectivity)[0])
        else:
            self.nu     = nu
            
        if np.shape(rho) == ():
            self.rho    = young*np.ones(np.shape(connectivity)[0])
        else:
            self.rho    = rho
            
            
        self.G                  = self.E / (2*(1 + self.nu))
        self.nNodes             = len(connectivity[0])
        self.totalDofs          = self._nDofs
        self.totalNodes         = len(coordinates)
        self.totalElements      = len(connectivity)
        
        # Properties for the warping boundary value problem
        self.activeDof          = None
        self.prescribedDof      = None
        self.stiffness          = lil_matrix((self.totalDofs,self.totalDofs)) # List of Lists sparse matrix
        # self.stiffness          = np.zeros([self.totalDofs,self.totalDofs])
        self.displacements      = np.zeros([self.totalDofs,1])
        self.forces             = np.zeros([self.totalDofs,1])
        self.phiStar            = np.zeros([self.totalDofs,1])
        self.areaProperties     = {}
        self.massProperties  = {}
        self.intDegree          = intDegree
        
        # List for timing
        self.times              = []
        

        # Element type
        self.elemType           = ''
        if self.nNodes in [3,6]:
            self.elemType       = 'Tri'
        elif self.nNodes in [4,8,9]:
            self.elemType       = 'Quad' 
        else:
            raise Exception('Element type not implemented: {sel}')

        # Simulate if noSimulation flag is off
        if noSimulation:
            return
        else:
            self.simulate(displayTimes=displayTimes)

    
    
    ## =================================================== ## 
    ##                   Private Methods                   ## 
    ## =================================================== ## 
    
    ## TORSIONAL CONSTANT
    def __calcTorsionalConstant(self, degree=2):
        '''
        Calculates the torsional constant (J) of the section.

        Parameters
        ----------
        degree : int, optional
            Quadrature degree for numerical integration. The default is 2.
        '''
        points, weights = self._elements[0].getQuadrature(degree)
        E0 = self.E[0] # np.min(self.E)
        G0 = self.G[0] # np.min(G)
        
        # Reads element type
        elemType = self._elements[0].elemType
        
        # Initializing torsional constant
        J_phi           = 0
        J_phi_w_G = 0
        for e,elem in enumerate(self._elements):

            # Degrees of freedom for the element
            elementDof = self.getElemDof(elem)
            # elementDof = np.array(self._elementsNodes)[e, :]

            # Get coordinates w.r. to the Shear Center (CT)
            ys = self._nodeCoords[elementDof,1] - self.areaProperties['Y_CT']
            zs = self._nodeCoords[elementDof,2] - self.areaProperties['Z_CT']
            # ys = Ys[elementDof] - self.areaProperties['Y_CT']
            # zs = Zs[elementDof] - self.areaProperties['Z_CT']

            # Finding the element's warping vector
            phi = self.displacements[elementDof]#.reshape(-1,1)
            # print(f'{e = }')
            # print(f'{elementDof = }')
            # print(f'{self.displacements.shape = }')
            # print(f'{phi.shape = }')
            # print('')
            
            # Material properties
            G   = self.G[e]

            # Quadrature for the element; Iterates through Gaussian points
            if elemType == 'Tri':
                for i, coords in enumerate(points):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)

                    # Shape functions and Jacobian
                    N = elem.getShapeFun(coords)
                    J, By, Bz = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point


                    # Integrated Function _______________
                    Fy = N@ys.reshape(-1,1)
                    Fz = N@zs.reshape(-1,1)

                    FF = Fy@Fy + Fz@Fz + Fy@Bz@phi - Fz@By@phi

                    J_phi += weights[i]*J*FF
                    J_phi_w_G += weights[i]*J*FF * G/G0
                    

            elif elemType == 'Quad':
                for j, coord_j in enumerate(points):
                    for i, coord_i in enumerate(points):
                        #      i: index
                        # coord: natural coordinates of the point (zeta1,2,3)
                        coords = [coord_i, coord_j]
                        
                        # Shape functions and Jacobian
                        N = elem.getShapeFun(coords)
                        J, By, Bz = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point


                        # Integrated Function _______________
                        Fy = N@ys.reshape(-1,1)
                        Fz = N@zs.reshape(-1,1)
                        FF = Fy@Fy + Fz@Fz + Fy@Bz@phi - Fz@By@phi

                        J_phi += weights[i]*weights[j]*J*FF
                        J_phi_w_G += weights[i]*weights[j]*J*FF * G/G0
            else:
                raise Exception('Element type not recognized')
                
        self.areaProperties.update(J_phi= float(J_phi),
                                   J_phi_w_G = float(J_phi_w_G))

        return
    ## AREA PROPERTIES
    def __calcAreaProperties(self): 
        '''
        Calculates the area properties of a composite section.
        '''

        points, weights = self._elements[0].getQuadrature(2)
        A, Qy, Qz, Iyy, Izz, Iyz = 0,0,0,0,0,0
        A_w, Qy_w , Qz_w , Izz_w, Iyy_w, Iyz_w = 0,0,0,0,0,0
        A_wG, Qy_wG , Qz_wG , Izz_wG, Iyy_wG, Iyz_wG = 0,0,0,0,0,0

        # Reference moduli
        E0 = self.E[0]
        G0 = self.G[0]

        # Reads element type
        elemType = self._elements[0].elemType        

        # for e in range(self.totalElements):
        for e,elem in enumerate(self._elements):

            # Degrees of freedom for the element
            elementDof = self.getElemDof(elem)

            # Finding the node coordinates
            ys = self._nodeCoords[elementDof,1]; 
            zs = self._nodeCoords[elementDof,2]

            # Quadrature for the element; Iterates through Gaussian points
            if elemType == 'Tri':
                for i, coords in enumerate(points):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)

                    # Shape functions and Jacobian
                    N = elem.getShapeFun(coords); #print(N)
                    J, _, _ = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                    
                    # AREA _______________________________
                    F = 1
                    A   += weights[i]*J*F

                    # STATIC MOMENTS OF AREA _________
                    Fy = N@ys.reshape(-1,1)
                    Qy  += weights[i]*J*Fy

                    Fz = N@zs.reshape(-1,1)
                    Qz  += weights[i]*J*Fz

                    # AREA MOMENTS OF INERTIA ________
                    Iyy += weights[i]*J*(Fy@Fy)

                    Izz += weights[i]*J*(Fz@Fz)

                    # PRODUCT OF INERTIA _________________
                    Iyz += weights[i]*J*(Fy@Fz)
                    
                    # WEIGHTED PROPERTIES ____________                    
                    E_weight = self.E[e] / E0
                    G_weight = self.G[e] / G0
                    
                    ## WEIGHTING (BY E)
                    A_w      += weights[i]*J*F        * E_weight
                    Qy_w     += weights[i]*J*Fy       * E_weight
                    Qz_w     += weights[i]*J*Fz       * E_weight
                    Iyy_w    += weights[i]*J*(Fy@Fy)  * E_weight
                    Izz_w    += weights[i]*J*(Fz@Fz)  * E_weight
                    Iyz_w    += weights[i]*J*(Fy@Fz)  * E_weight

                    ## WEIGHTING (BY G)
                    A_wG     += weights[i]*J*F        * G_weight
                    Qy_wG    += weights[i]*J*Fy       * G_weight
                    Qz_wG    += weights[i]*J*Fz       * G_weight
                    Iyy_wG   += weights[i]*J*(Fy@Fy)  * G_weight
                    Izz_wG   += weights[i]*J*(Fz@Fz)  * G_weight
                    Iyz_wG   += weights[i]*J*(Fy@Fz)  * G_weight

            elif elemType == 'Quad':
                for j, coord_j in enumerate(points):
                    for i, coord_i in enumerate(points):
                        #      i: index
                        # coord: natural coordinates of the point (zeta1,2,3)
                        coords = [coord_i, coord_j]

                        # Shape functions and Jacobian
                        N = elem.getShapeFun(coords); #print(N)
                        J, _, _ = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                        
                        # AREA _______________________________
                        F    = 1
                        A   += weights[i]*weights[j]*J*F

                        # STATIC MOMENTS OF AREA _________
                        Fy   = N@ys.reshape(-1,1)
                        Qy  += weights[i]*weights[j]*J*Fy

                        Fz   = N@zs.reshape(-1,1)
                        Qz  += weights[i]*weights[j]*J*Fz

                        # AREA MOMENTS OF INERTIA ________
                        Iyy += weights[i]*weights[j]*J*(Fy@Fy)

                        Izz += weights[i]*weights[j]*J*(Fz@Fz)

                        # PRODUCT OF INERTIA _________________
                        Iyz += weights[i]*weights[j]*J*(Fy@Fz)

                        # WEIGHTED PROPERTIES ____________
                        # Elastic moduli
                        E_weight = self.E[e] / E0
                        G_weight = self.G[e] / G0

                        ## WEIGHTING (BY E)
                        A_w      += weights[i]*weights[j]*J*F        * E_weight
                        Qy_w     += weights[i]*weights[j]*J*Fy       * E_weight
                        Qz_w     += weights[i]*weights[j]*J*Fz       * E_weight
                        Iyy_w    += weights[i]*weights[j]*J*(Fy@Fy)  * E_weight
                        Izz_w    += weights[i]*weights[j]*J*(Fz@Fz)  * E_weight
                        Iyz_w    += weights[i]*weights[j]*J*(Fy@Fz)  * E_weight

                        ## WEIGHTING (BY G)
                        A_wG     += weights[i]*weights[j]*J*F        * G_weight
                        Qy_wG    += weights[i]*weights[j]*J*Fy       * G_weight
                        Qz_wG    += weights[i]*weights[j]*J*Fz       * G_weight
                        Iyy_wG   += weights[i]*weights[j]*J*(Fy@Fy)  * G_weight
                        Izz_wG   += weights[i]*weights[j]*J*(Fz@Fz)  * G_weight
                        Iyz_wG   += weights[i]*weights[j]*J*(Fy@Fz)  * G_weight

   
        ## Determine the weighted centroid (always by E)
        Y_CG_w = float(Qy_w / A_w)
        Z_CG_w = float(Qz_w / A_w)
        ## Bending properties at the Centroid (parallel axis theorem)
        Iyy_centroid = Iyy_w - A_w * Y_CG_w**2
        Izz_centroid = Izz_w - A_w * Z_CG_w**2
        Iyz_centroid = Iyz_w - A_w * Y_CG_w*Z_CG_w

        # Without weighting (original properties)
        # self.areaProperties.update(Area=A, Qy=float(Qy), Qz=float(Qz), Iyy=float(Iyy), Izz=float(Izz), Iyz=float(Iyz), Y_CG=Y_CG, Z_CG=Z_CG)

        # Weighting by E
        self.areaProperties.update(
                                 A_w        = float(A_w), 
                                 Qy_w       = float(Qy_w), 
                                 Qz_w       = float(Qz_w), 
                                 Iyy_w      = float(Iyy_w), 
                                 Izz_w      = float(Izz_w), 
                                 Iyz_w      = float(Iyz_w), 
                                 Y_CG_w     = Y_CG_w, 
                                 Z_CG_w     = Z_CG_w
                               )

        # Weighting by G (These are mostly for internal calculations, not typically section properties)
        # self.areaProperties.update(Area_G=A_wG, Qy_G=float(Qy_wG), Qz_G=float(Qz_wG), Iyy_G=float(Iyy_wG), Izz_G=float(Izz_wG), Iyz_G=float(Iyz_wG), Y_CG=Y_CG_w, Z_CG=Z_CG_w)
        # Bending properties
        self.areaProperties.update(
                                 Iyy_w_cent = float(Iyy_centroid), 
                                 Izz_w_cent = float(Izz_centroid), 
                                 Iyz_w_cent = float(Iyz_centroid)
                               )
        return
    ## INERTIA AREA PROPERTIES
    def __calcMassProperties(self): 
        '''
        Calculates the inertia properties (mass moments) of a composite section.
        '''

        points, weights = self._elements[0].getQuadrature(2)
        rho_l           , Qy_rho    , Qz_rho          = 0,0,0     #      ρ dA,     ρ*y dA,     ρ*z dA
        Izz_rho         , Iyy_rho   , Iyz_rho         = 0,0,0     #  ρ*z^2 dA,   ρ*y^2 dA,   ρ*y*z dA
        phi_rho         , phi_phi_rho               = 0,0       #    ρ*φ dA,   ρ*φ^2 dA
        phi_y_rho       , phi_z_rho                 = 0,0       #  ρ*φ*y dA,   ρ*φ*z dA
        phi_yy_rho      , phi_yz_rho, phi_zz_rho    = 0,0,0     # ρ*φ*y^2 dA, ρ*φ*z^2 dA, ρ*φ*y*z dA
        
        Qy_cent_w, Qz_cent_w, Iyz_cent_weighted = 0,0,0 # E/E0*y dA,  E/E0*y dA, E/E0*y*z dA
        
        # Getting the centroid position
        Y_CG = self.areaProperties['Y_CG_w']
        Z_CG = self.areaProperties['Z_CG_w']
        
        # Getting the shear center position
        Y_CT = self.areaProperties['Y_CT']
        Z_CT = self.areaProperties['Z_CT']
        
        # Reads element type
        elemType = self._elements[0].elemType
    
        # # Get valid nodes coordinates
        # vn = self._validNodes
        # Ys = self._nodeCoords[vn,1]
        # Zs = self._nodeCoords[vn,2]

        # for e in range(self.totalElements):
        for e,elem in enumerate(self._elements):

            # Degrees of freedom for the element
            elementDof = self.getElemDof(elem)
            # elementDof = np.array(self._elementsNodes)[e, :]
            # print(f'elementDof = {elementDof}')

            # Finding the node coordinates in the local reference frame (section centroid)
            ys = self._nodeCoords[elementDof,1] - Y_CG
            zs = self._nodeCoords[elementDof,2] - Z_CG
            
            # Coordinates relative to the CT
            yts = self._nodeCoords[elementDof,1] - Y_CT
            zts = self._nodeCoords[elementDof,2] - Z_CT
            
            # Determining the nodal values of the warping function
            phi = self.displacements[elementDof]

            # Element properties
            rho     = self.rho[e]
            E_weight  = self.E[e] / self.E[0]

            # Quadrature for the element; Iterates through Gaussian points
            if elemType == 'Tri':
                for i, coords in enumerate(points):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)

                    # Shape functions and Jacobian
                    N = elem.getShapeFun(coords); #print(N)
                    J, _, _ = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                    WJ = weights[i]*J

                    # INTERPOLATION OF Y, Z and PHI
                    Fy   = N@ys.reshape(-1,1)
                    Fz   = N@zs.reshape(-1,1)
                    Fty  = N@yts.reshape(-1,1)
                    Ftz  = N@zts.reshape(-1,1)
                    Fphi = N@phi.reshape(-1,1)

                    ## PROPERTIES 
                    rho_l           += WJ*1         * rho           #      ρ dA
                    Qy_rho          += WJ*Fy        * rho           #    ρ*y dA
                    Qz_rho          += WJ*Fz        * rho           #    ρ*z dA
                    Iyy_rho         += WJ*Fy@Fy     * rho           #  ρ*z^2 dA
                    Izz_rho         += WJ*Fz@Fz     * rho           #  ρ*y^2 dA
                    Iyz_rho         += WJ*Fy@Fz     * rho           #  ρ*y*z dA
                    
                    ## PHI-DEPENDENT PROPERTIES
                    phi_rho         += WJ*Fphi        * rho        #    ρ*φ dA
                    phi_phi_rho     += WJ*Fphi@Fphi   * rho        #  ρ*φ^2 dA
                    phi_y_rho       += WJ*Fty@Fphi    * rho        #  ρ*φ*y dA
                    phi_z_rho       += WJ*Ftz@Fphi    * rho        #  ρ*φ*z dA
                    phi_yy_rho      += WJ*Fty@Fty@Fphi * rho        # ρ*φ*y^2 dA
                    phi_zz_rho      += WJ*Ftz@Ftz@Fphi * rho        # ρ*φ*z^2 dA
                    phi_yz_rho      += WJ*Fty@Ftz@Fphi * rho        # ρ*φ*y*z dA
                    
                    ## UPDATE STATIC MOMENTS (AT WEIGHTED CENTROID)
                    Qy_cent_w    += WJ*Fy        * E_weight  #   E/E0*y dA
                    Qz_cent_w    += WJ*Fz        * E_weight  #   E/E0*z dA
                    Iyz_cent_weighted   += WJ*Fy@Fz     * E_weight  #     y*z dA

            elif elemType == 'Quad':
                for j, coord_j in enumerate(points):
                    for i, coord_i in enumerate(points):
                        #      i: index
                        # coord: natural coordinates of the point (zeta1,2,3)
                        coords = [coord_i, coord_j]

                        # Shape functions and Jacobian
                        N = elem.getShapeFun(coords); #print(N)
                        J, _, _ = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                        WJ = weights[i]*weights[j]*J  # Product of weights and Jacobian
                        
                        # INTERPOLATION OF Y, Z and PHI
                        Fy   = N@ys.reshape(-1,1)
                        Fz   = N@zs.reshape(-1,1)
                        Fty  = N@yts.reshape(-1,1)
                        Ftz  = N@zts.reshape(-1,1)
                        Fphi = N@phi.reshape(-1,1)

                        ## PROPERTIES 
                        rho_l           += WJ*1         * rho           #      ρ dA
                        Qy_rho          += WJ*Fy        * rho           #    ρ*y dA
                        Qz_rho          += WJ*Fz        * rho           #    ρ*z dA
                        Iyy_rho         += WJ*Fy@Fy     * rho           #  ρ*z^2 dA
                        Izz_rho         += WJ*Fz@Fz     * rho           #  ρ*y^2 dA
                        Iyz_rho         += WJ*Fy@Fz     * rho           #  ρ*y*z dA
                        
                        ## PHI-DEPENDENT PROPERTIES
                        phi_rho         += WJ*Fphi        * rho        #    ρ*φ dA
                        phi_phi_rho     += WJ*Fphi@Fphi   * rho        #  ρ*φ^2 dA
                        phi_y_rho       += WJ*Fty@Fphi    * rho        #  ρ*φ*y dA
                        phi_z_rho       += WJ*Ftz@Fphi    * rho        #  ρ*φ*z dA
                        phi_yy_rho      += WJ*Fty@Fty@Fphi * rho        # ρ*φ*y^2 dA
                        phi_zz_rho      += WJ*Ftz@Ftz@Fphi * rho        # ρ*φ*z^2 dA
                        phi_yz_rho      += WJ*Fty@Ftz@Fphi * rho        # ρ*φ*y*z dA
                        
                        ## UPDATE STATIC MOMENTS (AT WEIGHTED CENTROID)
                        Qy_cent_w    += WJ*Fy        * E_weight  #    E/E0*y dA
                        Qz_cent_w    += WJ*Fz        * E_weight  #    E/E0*z dA
                        Iyz_cent_weighted   += WJ*Fy@Fz     * E_weight  #  E/E0*y*z dA
                        
                        
        # UPDATING INERTIA PROPERTIES
        self.massProperties.update(   
                                         rho_l          = float(rho_l        ),
                                         Qy_rho         = float(Qy_rho       ),
                                         Qz_rho         = float(Qz_rho       ),
                                         Iyy_rho        = float(Iyy_rho      ),
                                         Izz_rho        = float(Izz_rho      ),
                                         Iyz_rho        = float(Iyz_rho      ),
                                         I_phi_rho      = float(phi_rho      ),
                                         I_phi_phi_rho  = float(phi_phi_rho  ),
                                         I_phi_y_rho    = float(phi_y_rho    ),
                                         I_phi_z_rho    = float(phi_z_rho    ),
                                         I_phi_yy_rho   = float(phi_yy_rho   ),
                                         I_phi_zz_rho   = float(phi_zz_rho   ),
                                         I_phi_yz_rho   = float(phi_yz_rho   )
                                     )
        # UPDATING AREA PROPERTIES
        self.areaProperties.update(
                                         Qy_cent_w = float(Qy_cent_w),
                                         Qz_cent_w = float(Qz_cent_w)
        )
        return
    ## PROCESSING 
    def __globalStiffness (self, degree=4):
        '''
        Assembles the global stiffness matrix.

        Parameters
        ----------
        degree : int, optional
            Quadrature degree for numerical integration. The default is 4.

        Notes
        -----
        Call this method on the problem instance.
        
        Example
        -------
        >>> elementNodes = numpy.array([[0,1,2], [2,3,0]])
        >>> nodeCoords   = numpy.array([[0,0],[1,0],[1,1],[0,1]])
        >>> problem      = Torsion(elementNodes, nodeCoords)
        >>> problem.globalStiffness()
        '''
        # Loops for the elements
        for e, elem in enumerate(self._elements):

            #print(f'Elem.: {e}')
    
            # Degrees of freedom for the element
            elementDof   = self.getElemDof(elem)
            
            # nodesIndices = np.array(elem.getNodeIndices(),dtype=int)
            # elementDof   = np.array(self._dofs)[nodesIndices]
            # elementDof   = np.array(elementDof,dtype=int)
            
            # Filter only non nans and convert to int
            # elementDof = elementDof[np.logical_not(np.isnan(elementDof))]
            # print(elementDof)
            
            # elementDof = np.array(elem.getNodeIndices())
            # elementDof = np.array(self._elementsNodes)[e, :]

            # Finding the node coordinates
            # ys = self._nodeCoords[elementDof,1]
            # zs = self._nodeCoords[elementDof,2]

            # Element local stiffness matrix
            kLocal = self.G[e]*elem.getKlocal(degree)
            
            # Assembly of the global stiffness matrix
            cols, lins       = np.meshgrid(range(self.nNodes),range(self.nNodes))
            dofCol, dofLin   = np.meshgrid(elementDof,elementDof)
            self.stiffness[dofLin,dofCol] += kLocal[lins,cols]

            # Assembly of the Global Force Vector
            fLocal = self.G[e]*elem.getFlocal(degree)

            self.forces[dofLin,:] += fLocal[lins,:]

        
        return
    # BOUNDARY CONDITIONS
    def __boundaryConditions (self, prescribedDof):
        '''
        Calculates the free degrees of freedom (DOFs), given the constrained DOFs.

        Parameters
        ----------
        prescribedDof : np.ndarray
            A numpy array with the constrained node indices (DOFs).

        Example
        -------
        >>> prescribedDof = np.array([[1, 3, 4]])
        >>> problem.boundaryConditions(prescribedDof)
        '''
        
        # Nodes with prescribed displacement
        self.prescribedDof = prescribedDof # adjusting for Python
        
        # Nodes with free displacement
        # Note: self._elementsNodes is a list of lists of indices, so np.setdiff1d expects a flattened array.
        all_dofs = np.arange(self.totalDofs)
        self.activeDof = np.setdiff1d(all_dofs, self.prescribedDof)
        return
    # SOLVE FOR THE WARPING DISTRIBUITION
    def __solveWarping(self, degree):
        '''
        Solves [K]{d} = {f} for the free nodes.

        Parameters
        ----------
        degree : int
            Quadrature degree for numerical integration.

        Notes
        -----
        Call this method on the problem instance.
        
        Example
        -------
        >>> problem.solveWarping()
        '''
        
        free_dofs = self.activeDof # nodes without restriction
        gr1,gr2 = np.meshgrid(free_dofs,free_dofs)

        stiffTemp  = self.stiffness[gr1,gr2]
        forcesTemp = self.forces[free_dofs]

        ################################################################
        # System solution
        # SCR is best for solving
        tic      = time.perf_counter()    
        dispTemp =  spsolve(stiffTemp.tocsr(), forcesTemp)
        toc      = time.perf_counter()    
        self.times.append(f'System solve time:                   {toc-tic:.3f} seconds')
        
        
        # Filling the displacement vector with the free displacements
        j = 0
        for i in self.activeDof:
            self.phiStar[i, 0] = dispTemp[j]
            j += 1
                
        # Solving the force vector (reaction forces)
        self.forces = self.stiffness@self.phiStar

        ### Solution for the Shear Center (CT)
        tic      = time.perf_counter()    
        self.__solveShearCenter(degree)
        toc      = time.perf_counter()    
        self.times.append(f'Shear center integration time:       {toc-tic:.3f} seconds')
        
        
        ### Solution correction (to get the actual warping function $\varphi$)
        Y_CT = self.areaProperties['Y_CT']
        Z_CT = self.areaProperties['Z_CT']
        c_CT = self.areaProperties['c_CT']
        
        # # Get valid nodes coordinates
        vn = self._validNodes
        # Ys = self._nodeCoords[vn,1]
        # Zs = self._nodeCoords[vn,2]
        
        # ys   = Ys.reshape(-1,1)
        # zs   = Zs.reshape(-1,1)
        Ys   = np.array(self._nodeCoords[vn,1]).reshape(-1,1)
        Zs   = np.array(self._nodeCoords[vn,2]).reshape(-1,1)


        # The corrected warping function phi: $\varphi = \varphi^* + Y_{CT}z - Z_{CT}y - c_{CT}$
        self.displacements = self.phiStar + Y_CT*Zs - Z_CT*Ys - c_CT 

        ### Calculate torsional constant
        tic      = time.perf_counter()    
        self.__calcTorsionalConstant(degree)
        toc      = time.perf_counter()    
        self.times.append(f'Torsional constant integration time: {toc-tic:.3f} seconds')
        
        return
    # FIND SHEAR CENTER
    def __solveShearCenter(self, degree):
        '''
        Finds the position of the shear center (CT) and the correction constant (c_CT).

        Parameters
        ----------
        degree : int
            Quadrature degree for numerical integration.
        '''
        points, weights = self._elements[0].getQuadrature(degree)                     ############################
        Iphi, Iphiy, Iphiz = 0,0,0
        Iphi_weighted , Iphiy_weighted, Iphiz_weighted = 0,0,0
        E0 = self.E[0]
        
        # Reads element type
        elemType = self._elements[0].elemType
        
        # # Get valid nodes coordinates
        # vn = self._validNodes
        # Ys = self._nodeCoords[vn,1]
        # Zs = self._nodeCoords[vn,2]

        # Integrals of phi*
        for e,elem in enumerate(self._elements):

            # Degrees of freedom for the element
            elementDof = self.getElemDof(elem)
            # elementDof = np.array(self._elementsNodes)[e, :]
            # print(f'elementDof = {elementDof}')

            # Finding the node coordinates
            ys = self._nodeCoords[elementDof,1]
            zs = self._nodeCoords[elementDof,2]

            # Finding the nodal values of phi* for the element
            phiStar = self.phiStar[elementDof]

            # Integrals of phi*
            if elemType == 'Tri':
                for i, coords in enumerate(points):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)

                    # Shape functions and Jacobian
                    N = elem.getShapeFun(coords)
                    J, _, _ = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                    
                    # phi* _______________________________
                    Fphi  = N@phiStar.reshape(-1,1)
                    Iphi += weights[i]*J*Fphi

                    # Y*phi* _______________________________
                    Fy     = N@ys.reshape(-1,1)
                    Iphiy += weights[i]*J*Fy@Fphi

                    # Z*phi* _______________________________
                    Fz     = N@zs.reshape(-1,1)
                    Iphiz += weights[i]*J*Fz@Fphi

            elif elemType == 'Quad':
                for j, coord_j in enumerate(points):
                    for i, coord_i in enumerate(points):
                        #      i: index
                        # coord: natural coordinates of the point (zeta1,2,3)
                        coords = [coord_i, coord_j]
                        N = elem.getShapeFun(coords)
                        J, _, _ = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                        
                        # phi* _______________________________
                        Fphi  = N@phiStar.reshape(-1,1)
                        Iphi += weights[i]*weights[j]*J*Fphi

                        # Y*phi* _______________________________
                        Fy     = N@ys.reshape(-1,1)
                        Iphiy += weights[i]*weights[j]*J*Fy@Fphi

                        # Z*phi* _______________________________
                        Fz     = N@zs.reshape(-1,1)
                        Iphiz += weights[i]*weights[j]*J*Fz@Fphi

                        ## WEIGHTING BY E
                        E_weight = self.E[e] / E0

                        ## WEIGHTING (BY E)
                        Iphi_weighted    += weights[i]*weights[j]*J*Fphi    * E_weight
                        Iphiy_weighted   += weights[i]*weights[j]*J*Fy@Fphi * E_weight
                        Iphiz_weighted   += weights[i]*weights[j]*J*Fz@Fphi * E_weight
            
        # Assembly of the system of equations
        A   = self.areaProperties['A_w']
        Qy  = self.areaProperties['Qy_w']
        Qz  = self.areaProperties['Qz_w']
        Izz = self.areaProperties['Izz_w']
        Iyy = self.areaProperties['Iyy_w']
        Iyz = self.areaProperties['Iyz_w']

        # Coefficient matrix
        M = np.array([[ -Qz, Qy, A],
                      [-Iyz,Iyy,Qy],
                      [-Izz,Iyz,Qz]])
        # Independent terms
        b = np.array([Iphi_weighted,Iphiy_weighted,Iphiz_weighted])

        # Solution
        x = np.linalg.inv(M)@b.reshape(-1,1)
        Y_CT, Z_CT, c_CT = x

        self.areaProperties.update(Y_CT=float(Y_CT), Z_CT=float(Z_CT), c_CT=float(c_CT))
    
        return
    
    
    
    ## =================================================== ## 
    ##                   Public Methods                    ## 
    ## =================================================== ## 
    
    ## --------------------- Getters --------------------- ##
    ## GET AREA PROPERTIES
    def getAreaProperties(self) -> dict[str,float]:
        '''
        Gets the area properties of the section. 
        
        Note: The reference moduli are always taken to be the first entry on the material properties lists (E0 = E[0], G0 = G[0])

        Returns
        -------
        areaProps : dict
            Dictionary of area properties.
            
        | Key        	| Description                                                                                                                           	| Dimension  	    
        |------------	|---------------------------------------------------------------------------------------------------------------------------------------	|-------------	
        | A_w        	| Weighted Area (∫(E/E0)dA). Area weighted by the ratio of the elastic modulus to a reference modulus (E0).                             	| L^2          	
        | Qy_w       	| Weighted First Moment of Area about Z-axis (∫y(E/E0)dA). Used to locate the centroid.                                                 	| L^3          	
        | Qz_w       	| Weighted First Moment of Area about Y-axis (∫z(E/E0)dA). Used to locate the centroid.                                                 	| L^3          	
        | Y_CG_w     	| Y-coordinate of the Weighted Centroid.                                                                                  	                | L           	
        | Z_CG_w     	| Z-coordinate of the Weighted Centroid.                                                                                  	                | L           	
        | Iyy_w      	| Weighted Second Moment of Area about Z-axis (∫y^2(E/E0)dA), w.r.t the global origin.                                                   	| L^4          	
        | Izz_w      	| Weighted Second Moment of Area about Y-axis (∫z^2(E/E0)dA), w.r.t the global origin.                                                   	| L^4          	
        | Iyz_w      	| Weighted Product of Inertia (∫yz(E/E0)dA), w.r.t the global origin.                                                                   	| L^4          	
        | Iyy_w_cent 	| Weighted Second Moment of Area about the Centroid z-axis.                                                                 	            | L^4          	
        | Izz_w_cent 	| Weighted Second Moment of Area about the Centroid y-axis.                                                                 	            | L^4          	
        | Iyz_w_cent 	| Weighted Product of Inertia about the Centroid.                                                                                 	        | L^4          	
        | J_phi      	| Torsional Constant. Calculated from the warping function φ.                                            	                                | L^4          	
        | J_phi_w_G  	| Torsional Constant Weighted by a reference Shear Modulus G0 (∫(G/G0)(…)dA).                                                               | L^4          	
        | Y_CT       	| Y-coordinate of the Shear Center (Y_CT). Calculated from the BVP correction.                                                           	| L           	
        | Z_CT       	| Z-coordinate of the Shear Center (Z_CT). Calculated from the BVP correction.                                                           	| L           	
        | c_CT       	| Warping Function Correction Constant (c_CT). Part of the final warping solution φ=φ∗+Y_CTz−Z_CTy−c_CT.                                    | L^2          	
        | Qy_cent_w  	| First Moment of Area about Z-axis at the Centroid (∫ycent(E/E0)dA). This value should theoretically be zero but is stored for checks. 	| L^3          	
        | Qz_cent_w  	| First Moment of Area about Y-axis at the Centroid (∫zcent(E/E0)dA). This value should theoretically be zero but is stored for checks. 	| L^3          	
            
        '''
        return self.areaProperties

        ## ----------------- Public Methods ----------------- ## 
    ## GET MASS PROPERTIES
    def getMassProperties(self) -> dict[str,float]:
        '''
        Gets the area mass properties of the section.
        
        Note: The reference moduli are always taken to be the first entry on the material properties lists (E0 = E[0], G0 = G[0])


        Returns
        -------
        massProps : dict
            Dictionary of area properties.
            
        | Key           	| Description                                                                       	| Dimension 	|
        |---------------	|-----------------------------------------------------------------------------------	|-----------	|
        | rho_l         	| Total Mass (Density-Weighted Area) (∫ ρ dA).                                      	| M         	|
        | Qy_rho        	| First Mass Moment about Z-axis (∫ y*ρ dA).                                        	| M *L      	|
        | Qz_rho        	| First Mass Moment about Y-axis (∫ z*ρ dA).                                        	| M *L      	|
        | Iyy_rho       	| Mass Moment of Inertia about Z-axis (∫ y^2*ρ dA).                                 	| M *L^2    	|
        | Izz_rho       	| Mass Moment of Inertia about Y-axis (∫ z^2*ρ dA).                                 	| M *L^2    	|
        | Iyz_rho       	| Mass Product of Inertia (∫ y*z*ρ dA).                                             	| M *L^2    	|
        | I_phi_rho     	| Warping Mass Moment (∫ φ*ρ dA).                                                   	| M *L^3    	|
        | I_phi_phi_rho 	| Second Warping Mass Moment of Inertia (∫ φ^2*ρ dA).                               	| M *L^4    	|
        | I_phi_y_rho   	| Warping-Y Mass Moment of Inertia (∫ φ*y*ρ dA), relative to the Shear Center (CT). 	| M *L^3    	|
        | I_phi_z_rho   	| Warping-Z Mass Moment of Inertia (∫ φ*z*ρ dA), relative to the Shear Center (CT). 	| M *L^3    	|
        | I_phi_yy_rho  	| Warping-Y^2 Mass Moment of Inertia (∫ φ*y^2*ρ dA).                                	| M *L^4    	|
        | I_phi_zz_rho  	| Warping-Z^2 Mass Moment of Inertia (∫ φ*z^2*ρ dA).                                	| M *L^4    	|
        | I_phi_yz_rho  	| Warping-YZ Mass Product of Inertia (∫ φ*y*z*ρ dA).                                	| M *L^4    	|
        '''
        return self.massProperties
    
    ## ------------------ Other methods ------------------ ##    
    ## CALCULATE SHEAR STRESSES AT QUADRATURE POINTS
    def calcShearStresses(self, degree:int=4, twistRate:float=None, mises:bool=False) -> tuple[np.ndarray, np.ndarray]:
        '''Calculate shear stresses at quadrature points of a given degree.
        
        Parameters
        ----------
        degree : int, optional
            Quadrature degree for the sampling points. The default is 4.
        twistRate : float, optional
            Twist rate, if not given stresses are calculated per units of twist rate.
        mises : bool, optional
            Caculate stress magnitude as mises equivalent stress. Default: True
        
        Returns
        -------
        stresses : numpy.ndarray
            XY and XZ shear stress components at gauss points.
        gauss_points : numpy.ndarray
            Position of the gauss points.
        '''
        points, _   = self._elements[0].getQuadrature(degree)   # Calculating at Gaussian points
        # points, _   = [-1, 0, 1]                           # Calculating at nodes
        # points      = [-1, 0, 1] if POS == 'nodes' else self._elements[0].quadrature(degree)
        
        # Initializing stress and gaussian points lists 
        tau_xy      = np.array([])
        tau_xz      = np.array([])
        Y_gauss     = np.array([])
        Z_gauss     = np.array([])
        
        # Get shear center position
        Y_CT = self.areaProperties['Y_CT']
        Z_CT = self.areaProperties['Z_CT']
        
        # Reads element type
        elemType = self._elements[0].elemType
        
        # Iterate through elements
        for e,elem in enumerate(self._elements):

            # Element degrees of freedom 
            elementDof = self.getElemDof(elem)
            # elementDof = np.array(self._elementsNodes)[e, :]

            # Node coordinates w.r.t the shear center
            ys = self._nodeCoords[elementDof,1] - Y_CT
            zs = self._nodeCoords[elementDof,2] - Z_CT
            
            # Nodal warping displacements for the element
            phi = self.displacements[elementDof]
            
            # Element shear modulus
            G   = self.G[e]

            # Quadrature for the element; Iterates through Gaussian points
            if elemType == 'Tri':
                for i, coords in enumerate(points):
                    #      i: index
                    # coord: natural coordinates of the point (zeta1,2,3)

                    # Shape functions and Jacobian
                    N = elem.getShapeFun(coords)
                    J, By, Bz = elem.getJacob(coords)
                    
                    # Normalized tau_xy _______________________________
                    z_gp    = N@zs                       # Gaussian point position w.r.t CT
                    txy     = By@phi - z_gp
                    tau_xy  = np.append(tau_xy, txy)
                    Z_gauss = np.append(Z_gauss, z_gp)

                    # Normalized tau_xz _______________________________
                    y_gp    = N@ys                       # Gaussian point position w.r.t CT
                    txz     = Bz@phi + y_gp
                    tau_xz  = np.append(tau_xz, txz)
                    Y_gauss = np.append(Y_gauss, y_gp)

            elif elemType == 'Quad':
                for j, coord_j in enumerate(points):
                    for i, coord_i in enumerate(points):
                        #      i: index
                        # coord: natural coordinates of the point (zeta1,2,3)
                        coords = [coord_i, coord_j]
                        
                        # Shape functions and Jacobian
                        N = elem.getShapeFun(coords)
                        J, By, Bz = elem.getJacob(coords) # (ys,zs): element nodes, coords: zeta coordinates of the Gaussian point
                            
                        # tau_xy _______________________________
                        Z       = N@zs
                        z_gp    = Z                          # Gaussian point position w.r.t CT
                        txy     = (By@phi - z_gp)*G
                        tau_xy  = np.append(tau_xy, txy)
                        Z_gauss = np.append(Z_gauss, Z)

                        # tau_xz _______________________________
                        Y       = N@ys
                        y_gp    = Y                          # Gaussian point position w.r.t CT
                        txz     = (Bz@phi + y_gp)*G
                        tau_xz  = np.append(tau_xz, txz)
                        Y_gauss = np.append(Y_gauss, Y)

                #      print('End i')
                # print('End j')
        
        
        if mises:
            tau_xy = np.sqrt(3)*tau_xy
            tau_xz = np.sqrt(3)*tau_xz        
        if twistRate is not None:
            tau_xy = twistRate*tau_xy
            tau_xz = twistRate*tau_xz
        
        # Pack results into tuples
        stresses    = tau_xy, tau_xz
        gaussPoints = Y_gauss, Z_gauss
                
        return stresses, gaussPoints
    ## PLOT WARPING DISPLACEMENT
    def plotWarping(self, 
                    levels      :int    = 20, 
                    twistRate   :float  = None, 
                    fignum      :int    = None, 
                    figsize     :tuple  = None,
                    lengthUnits :str    = 'l.u.', 
                    showMesh    :bool   = False, 
                    showCT      :bool   = True, 
                    showCG      :bool   = True, 
                    cbarKwargs  :dict   = {},
                    fontname    :str    = 'serif',
                    **kwargs):
        '''
        Plots the warping function ($\varphi$) distribution or the warping displacement ($\alpha\varphi$) 
        on the cross-section as a contour plot. 
        Results are interpolated from nodal values by the matplotlib.pyplot.tricontourf function.

        Parameters
        ----------
        levels : int, optional
            Number of contour curve divisions. The default is 20.
        twistRate : float, optional
            Torsion rate ($\alpha$). If provided, the warping displacement ($\alpha\varphi$) 
            is plotted instead of the warping function ($\varphi$). The default is None.
        fignum : int, optional
            Figure number. The default is None.
        figsize : tuple, optional
            Figure size (width, height). The default is None.
        lengthUnits : str, optional
            Length unit used in the section's coordinates. The default is 'l.u.' (length units).
        showMesh : bool, optional
            Draws the underlying mesh elements. The default is False.
        showCT : bool, optional
            Marks the Shear Center (CT) on the plot. The default is True.
        showCG : bool, optional
            Marks the Centroid (CG) on the plot. The default is True.
        cbarKwargs : dict, optional
            Keyword arguments to pass to the colorbar function.
        **kwargs : dict
            Additional keyword arguments passed to the primary plotting function (matplotlib.pyplot.tricontourf).

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated figure object.
        ax : matplotlib.axes.Axes
            The axes object with the plotted warping results.
        '''

        phi = np.array(self.displacements).reshape(-1)
        if twistRate != None:
            phi = twistRate*phi
            # label = fr'$\mathrm{{Warping \ Displacement }} \ \alpha\varphi^h$ [\mathrm{{{unit}}}]'+'\n'+ fr'$\alpha = {twistRate:.3e} \ [\mathrm{{rad/{{{unit}}}}}]$'
            label = f'Warping Displacement $\\alpha\\varphi^h$ [{lengthUnits}]'+'\n'+ f'$\\alpha$ = {twistRate:.3e} [rad/{lengthUnits}]'
        else:
            # label = fr'$\mathrm{{Warping \ Function }} \ \varphi^h \ \left[ \mathrm{{{unit}}}^2 \right]$'
            label = f'Warping Function $\\varphi^h$ [{lengthUnits}²/rad]'
        
        # Get all coordinates
        Y   = np.array(self._nodeCoords)[:,1].reshape(-1)
        Z   = np.array(self._nodeCoords)[:,2].reshape(-1)
        
        # Create figure and subplot
        fig, ax = plt.subplots(num=fignum,figsize=figsize)
        
        # Font options
        fontkwargs = dict(fontname=fontname, usetex=False, fontsize=18)
        # fontkwargs = dict(fontname='Times New Roman', usetex=True, fontsize=18)
        
        # Deviding the element into triangular regions, suitable for use in tricontouf
        element_connectivity = np.array(self._elementsNodes)
        match self.nNodes:
            case 3:
                t0 = element_connectivity[:,[0,1,2]]
                connect = t0
            case 6:
                t0 = element_connectivity[:,[5,4,2]]
                t1 = element_connectivity[:,[0,4,5]]
                t2 = element_connectivity[:,[0,3,4]]
                t3 = element_connectivity[:,[3,1,4]]
                connect = np.vstack([t0,t1,t2,t3])
            case 4:
                t0 = element_connectivity[:,[0,1,2]]
                t1 = element_connectivity[:,[2,3,0]]
                connect = np.vstack([t0,t1])
            case 9:
                t0 = element_connectivity[:,[7,8,3]]
                t1 = element_connectivity[:,[0,8,7]]
                t2 = element_connectivity[:,[0,4,8]]
                t3 = element_connectivity[:,[4,1,8]]
                t4 = element_connectivity[:,[1,5,8]]
                t5 = element_connectivity[:,[8,5,2]]
                t6 = element_connectivity[:,[8,2,6]]
                t7 = element_connectivity[:,[8,6,3]]
                connect = np.vstack([t0,t1,t2,t3,t4,t5,t6,t7])
        
        # Padding displacements to account for non used nodes
        phiNew = np.nan*np.ones(self._nodeCoords.shape[0])
        validDofs = np.where(np.logical_not(np.isnan(self._dofs)))
        phiNew[validDofs] = phi
        # phiNew = phiNew.reshape(-1)
        
        # Create triangulation for that type of element
        triang  = tri.Triangulation(Y,Z,connect)

        # Set colobar ticks to 10 equaly spaced values between min and max magnitudes
        minMag = np.min(phiNew)
        maxMag = np.max(phiNew)
        if 'ticks' not in cbarKwargs:
            cbarKwargs['ticks'] = np.linspace(minMag,maxMag,11)
        
        # Create contour plot
        contour = ax.tricontourf(triang, phiNew, cmap='jet',levels=levels)
        cbar    = plt.colorbar(contour, **cbarKwargs)
        cbar.set_label(label, **fontkwargs)

        # Show mesh
        if showMesh:
            for elem in self._elements:
                # Node coordinates
                nodesCoords = elem.getNodeCoords()
                
                # Determining the vertex order
                if self.nNodes == 9:
                    # Case Q9
                    inds = [0,4,1,5,2,6,3,7]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                elif self.nNodes == 6:
                    # Case T6
                    inds = [0,3,1,4,2,5]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                else:
                    # Case linear elements
                    Ys = np.array(nodesCoords)[:,1]
                    Zs = np.array(nodesCoords)[:,2]

                # Vertices
                verts = np.vstack([Ys,Zs]).T

                # poly        = Polygon(verts,fc=fill_color,alpha=.2,ec='k')
                poly        = Polygon(verts,fill=False,alpha=1,ec='k',lw=.3)
                ax.add_patch(poly)
        # Show the Shear Center (CT)
        if showCT:
            ax.scatter(self.areaProperties['Y_CT'], self.areaProperties['Z_CT'], marker='^', edgecolors='r',facecolors='w', s=100, label='CT')
        # Show the Centroid (CG)
        if showCG:
            ax.scatter(self.areaProperties['Y_CG_w'], self.areaProperties['Z_CG_w'], marker='o', edgecolors='b',facecolors='w', s=100, label='CG')

        # ax.set_title(f'Cross Section',fontsize=14,fontweight='bold')
        # ax.set_xlim([np.min(Y),np.max(Y)])
        # ax.set_ylim([np.min(Z),np.max(Z)])
        ax.set_xlabel(f'$Y$ [{lengthUnits}]',**fontkwargs)
        ax.set_ylabel(f'$Z$ [{lengthUnits}]',**fontkwargs)
        ax.tick_params(axis='y',labelsize=12)
        ax.tick_params(axis='x',labelsize=12)
        plt.legend(loc='center right')
        # ax.axis('tight')
        ax.set_aspect('equal', adjustable='box')
        fig.tight_layout()
        
        return fig, ax
    ## PLOT CROSS-SECTION GEOMETRY
    def plotSection(self, 
                showElemLabel   =   False, 
                showNodeLabel   =   False, 
                showGaussPoints =   False,
                fillColor:str   =   '#59c1f9', 
                figsize:tuple   =   None, 
                fontsize        =   12) -> tuple[plt.Figure, plt.Axes]:
        '''
        Plots the cross-section geometry.

        Parameters
        ----------
        showElemLabel : bool, optional
            Show element labels, default False
        showNodeLabel : bool, optional
            Show node labels, default False
        showGaussPoints : bool, optional
            Show gauss point positions, default False
        fillColor : str, optional
            The color to fill the element polygons. The default is '#59c1f9'.
        figsize : tuple, optional
            Figure size, default: None
        fontsize: int, optional
            Reference font size. Element label and node label sizes are smaller. Default: 12.
        '''
        fig, ax      = plt.subplots(figsize=figsize)      
                   
        # Colors
        colors = [fillColor, '#FFA55C', '#BDFA96', '#FADD8E', '#EA8EFA',
                  '#A3E9FA', '#FAA2C1', '#FAF2A2', '#879FA5', '#7A7750']
        
        # Mesh element type
        elemType = self.elemType 
        # Gauss points notural coordinates
        points, _ = self._elements[0].getQuadrature(self.intDegree)
        # Lists to store gauss points
        Z_gauss = []
        Y_gauss = []
        
        # Track material change to update fill color
        uniqueMaterials = np.array([])
        newMaterial = False
        # colorIndex = -1
        
        # Loop through elements
        for e,elem in enumerate(self._elements):
            # Node coordinates
            nodesCoords = elem.getNodeCoords()
            
            # Read element material 
            material = self.G[e] 
                            
            # Check if material is new and update color index if it is
            newMaterial = False
            if material not in uniqueMaterials:
                uniqueMaterials = np.append(uniqueMaterials, material)
                newMaterial = True
                # colorIndex = int(colorIndex +1)
             
            # Get material color index
            try:
                colorIndex = int(np.where(uniqueMaterials==material)[0])
            except IndexError:
                # If there are more than 10 materials, set new ones to the first color
                colorIndex = 0
            # print(f'{colorIndex = }, {material = }')
            
            # Get node number in anti-clockwise order
            match len(nodesCoords):
                case 9:
                    connect = [0,4,1,5,2,6,3,7]
                case 6:
                    connect = [0,3,1,4,2,5]
                case 4:
                    connect = [0,1,2,3]
                case 3:
                    connect = [0,1,2]
                case _:
                    raise Exception('Invalid number of nodes')
            
            vertYs = np.array(nodesCoords)[connect,1]
            vertZs = np.array(nodesCoords)[connect,2]
            
            
            verts = np.vstack([vertYs,vertZs]).T

            # Vertices
            if newMaterial:
                # Create a new label for the legend for every new material
                poly        = Polygon(verts,fc=colors[colorIndex],alpha=.2,ec='k',zorder=0, label=f'$\mu$ = {material:.2e}')
            else:
                poly        = Polygon(verts,fc=colors[colorIndex],alpha=.2,ec='k',zorder=0)
            ax.add_patch(poly)
            
            # Show element number at the centroid of the first 3 node coordinates
            if showElemLabel:
                coordsTriangle = np.array(nodesCoords)[:3]
                triCentroid = 0, np.sum(coordsTriangle[:,1])/3, np.sum(coordsTriangle[:,2])/3 
                ax.text(triCentroid[1],triCentroid[2], elem.getIndex(), 
                        fontsize=fontsize-2,ha='center',va='center')
            
            # Show node numbers
            if showNodeLabel:
                nodesIndices = elem.getNodeIndices()
                for k,coord in enumerate(nodesCoords):
                    ax.text(coord[1],coord[2], nodesIndices[k], 
                            fontsize=fontsize-4, c='m', ha='left',va='top', style='italic')
            
            # Show gauss points
            if showGaussPoints:
                # Element degrees of freedom 
                elementDof = self.getElemDof(elem)

                # Node coordinates 
                ys = self._nodeCoords[elementDof,1]
                zs = self._nodeCoords[elementDof,2]


                # Quadrature for the element; Iterates through Gaussian points
                if elemType == 'Tri':
                    for coords in points:
                        #      i: index
                        # coord: natural coordinates of the point (zeta1,zeta2,zeta3)

                        # Shape functions and Jacobian
                        N = elem.getShapeFun(coords)
                        
                        # Z coordinate _______________________________
                        z_gp    = N@zs                       # Gaussian point position 
                        Z_gauss = np.append(Z_gauss, z_gp)

                        # Y coordinate _______________________________
                        y_gp    = N@ys                       # Gaussian point position 
                        Y_gauss = np.append(Y_gauss, y_gp)

                elif elemType == 'Quad':
                    for coord_j in points:
                        for coord_i in points:
                            #      i: index
                            # coord: natural coordinates of the point (zeta1,zeta2)
                            coords = [coord_i, coord_j]
                            
                            # Shape functions and Jacobian
                            N = elem.getShapeFun(coords)
                            
                            # Z coordinate _______________________________
                            z_gp    = N@zs
                            Z_gauss = np.append(Z_gauss, z_gp)

                            # Y coordinate _______________________________
                            y_gp    = N@ys
                            Y_gauss = np.append(Y_gauss, y_gp)

        
        # Plot gauss points
        if showGaussPoints:
            ax.scatter(Y_gauss,Z_gauss,c='k',alpha=.25,marker='x',s=2)
        
        Ys = self._nodeCoords[:,1]
        Zs = self._nodeCoords[:,2]

        ax.set_aspect('equal')
        ax.set_title(f'Cross Section',fontsize=14,fontweight='bold')
        ax.set_xlim([np.min(Ys),np.max(Ys)])
        ax.set_ylim([np.min(Zs),np.max(Zs)])
        ax.legend(title='Material',bbox_to_anchor=(1, 1))
        # plt.tight_layout()
        return fig, ax
    ## PLOT SHEAR STRESS DISTRIBUITION
    def plotShearStresses(self, 
                          degree        :int    = 2, 
                          twistRate     :float  = None, 
                          mises         :bool   = False,
                          mode          :str    = 'vector', 
                          stressUnits   :str    = 's.u.',
                          lengthUnits   :str    = 'l.u.',
                          showMesh      :bool   = True, 
                          vectorUnits   :bool   = True, 
                          vectorStep    :int    = 1, 
                          figsize       :tuple  = (6,4), 
                          cbarKwargs    :dict   = dict(),
                          cmap          :str    = 'jet',
                          fontname      :str    = 'serif',
                          **kwargs) -> tuple[plt.Figure, plt.Axes]:
        '''
        Plots shear stresses for visualization of the shear flow.

        Parameters
        ----------
        degree : int, optional
            Quadrature degree for calculation points. The default is 2.
        mode : str, optional
            Plotting mode:
            'vector' : shows the vector field of stresses.
            'scalar' : shows the scalar field of each component.
            The default is 'vector'.
        showMesh : bool, optional
            Show the mesh on the plot. The default is True.
        twistRate : float, optional
            Torsion rate ($\alpha$). If not provided, stresses are given per unit of alpha. The default is None.
        lengthUnits : str, optional
            Length unit of measure. The default is 'mm'.
        stressUnits : str, optional
            Stress unit of measure. The default is 'Mpa'.
        vectorUnits : bool, optional
            Use normalized vectors in vector mode. The default is True.
        vectorStep : int, optional
            Step size for sampling vectors in vector mode. For example, if vectorStep=2, every second vector is displayed. The default is 1.
        figsize : tuple, optional
            Figure size (width, height). The default is (6, 4).
        mises : bool, optional
            Plot the equivalent von Mises shear stress. The default is False.
        cbarKwargs : dict, optional
            Keyword arguments to pass to the plt.colorbar function. The default is dict().
        cmap : str, optional
            Color map for filled contour plot. Check possible options on ``matplotlib.colormaps``. Default: 'jet'. 
        fontname: str, optional
            Font for the plots.
        **kwargs : dict
            Additional keyword arguments for the plotting function (e.g., quiver or tricontourf).

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated figure object.
        ax : matplotlib.axes.Axes or list of Axes
            The axes object(s) used for plotting.
        '''
        # Calculate stresses and get gauss points
        stresses, gaussPoints   = self.calcShearStresses(degree, twistRate=twistRate, mises=mises)
        tau_xy, tau_xz          = stresses
        Y_gauss, Z_gauss        = gaussPoints
        
        # Stress magnitudes
        mags     = np.sqrt(tau_xy**2 + tau_xz**2)
        maxMag = np.max(mags)
        minMag = np.min(mags)
        
        # Shear center
        Y_CT = self.areaProperties['Y_CT']
        Z_CT = self.areaProperties['Z_CT']
        

        ## ------------------------- PLOTTING ------------------------- ##
        # Font options
        fontkwargs = dict(fontname=fontname, usetex=False, fontsize=18)
        
        # Determine colorbar label 
        if twistRate is not None:
            if mises:
                labelBar = r'$\sqrt{3}| \tau | $' + f' [{stressUnits}]'
            else:
                labelBar = r'$| \tau | $' + f' [{stressUnits}]'            
        else:
            if mises:
                labelBar = r'$\sqrt{3}| \tau | / \alpha $' + fr' $[\mathrm{{{stressUnits}}}\cdot\mathrm{{{lengthUnits}/rad}}]$'
            else:    
                labelBar = r'$| \tau | / \alpha $' + fr' $[\mathrm{{{stressUnits}}}\cdot\mathrm{{{lengthUnits}/rad}}]$'
            
        # Select the plot type: vector or scalar
        match mode:
            ## SHOW STRESSES AS VECTOR FIELD
            case 'vector':
                # Get gauss point positions
                Y       = Y_gauss
                Z       = Z_gauss
                
                # Create figure and axes
                fig, ax = plt.subplots(figsize=figsize, num='Shear Stresses Vector')
                
                # Set colobar ticks to 10 equaly spaced values between min and max magnitudes
                if 'ticks' not in cbarKwargs:
                    cbarKwargs['ticks'] = np.linspace(minMag,maxMag,11)                    
                
                # Calculate normilized vectors
                (U, V) = (tau_xy/mags, tau_xz/mags) if vectorUnits else (tau_xy, tau_xz)
                    
                # Plot vector field
                plot    = ax.quiver(Y[::vectorStep]+Y_CT, Z[::vectorStep]+Z_CT,
                                    U[::vectorStep], V[::vectorStep], 
                                    mags[::vectorStep], cmap=cmap, pivot='mid', 
                                    units='xy', **kwargs)
                
                # Create color bar
                cbar    = plt.colorbar(plot, **cbarKwargs)
                cbar.set_label(labelBar, **fontkwargs)

                # Add labels
                ax.set_xlabel(r'$Y$', **fontkwargs)
                ax.set_ylabel(r'$Z$', **fontkwargs)
                # ax.set_title('Shear Stress Distribution ($yz$-plane)'  , **fontkwargs)
                ax.set_aspect('equal', adjustable='box')
                # ax.axis('tight')

            ## COMPONENTS AS SCALAR FIELDS
            case 'scalar':

                # Putting vectors in n_elem x n_gp format
                Y_gauss     = np.reshape(Y_gauss,[self.totalElements,-1]) + Y_CT
                Z_gauss     = np.reshape(Z_gauss,[self.totalElements,-1]) + Z_CT

                # Create figure and axes
                fig, ax = plt.subplots(figsize=(figsize[0],2*figsize[1]), nrows=2, sharex=True, sharey=True, num='Shear Stresses Scalar')         
                
                
                # Creating a triangulation connectivity for gauss points, suitable for use in tricontouf
                nPoints = Y_gauss.shape[0]*Y_gauss.shape[1]                 # Total of gauss points
                gauss_connec = np.arange(nPoints).reshape(Y_gauss.shape)    # Connectivity for gausspoints
                
                # Triangulation for quads is straightforward
                #
                # (k + Ngp)  --------------- (k + Ngp + 1) 
                #    |                       .     |
                #    |      t2          .          |
                #    |             .               |
                #    |        .          t1        |
                #    |   .                         |
                #   (k) ------------------------ (k+1) 
                # k = n + Ngp*m
                # t1 = (k, k+1, k+Ngp+1)
                # t2 = (k, (k+Ngp+1, k+Ngp)
                if self.elemType == 'Quad':
                    # Number of gauss points per axis
                    Ngp = int(np.sqrt(Y_gauss.shape[1]))
                    # Initialization of the connectivy for the triangles
                    connectTrian = np.array([0,0,0])
                    for m,n in product(range(Ngp-1),range(Ngp-1)):
                        k = n + Ngp*m
                        t1 = gauss_connec[:,[k, k+1     ,k+Ngp+1]]
                        t2 = gauss_connec[:,[k, k+Ngp+1 ,k+Ngp]]
                        connectTrian = np.vstack([connectTrian, t1, t2])
                    # Delete first row
                    connectTrian = np.delete(connectTrian,0,axis=0)
                    
                else: 
                    # Triangulation for tris is complicated, has to be dealt in case by case
                    # because of gauss point enumeration
                    match Y_gauss.shape[1]:
                        case 3:
                            t0 = gauss_connec[:,[0,1,2]]
                            connectTrian = t0
                        case 6:
                            t0 = gauss_connec[:,[4,0,2]]
                            t1 = gauss_connec[:,[0,5,2]]
                            t2 = gauss_connec[:,[5,1,2]]
                            t3 = gauss_connec[:,[1,3,2]]
                            connectTrian = np.vstack([t0,t1,t2,t3])
                        case 7:
                            t0 = gauss_connec[:,[2,4,0]]
                            t1 = gauss_connec[:,[4,3,0]]
                            t2 = gauss_connec[:,[3,5,0]]
                            t3 = gauss_connec[:,[5,1,0]]
                            t4 = gauss_connec[:,[1,6,0]]
                            t5 = gauss_connec[:,[6,2,0]]
                            connectTrian = np.vstack([t0,t1,t2,t3,t4,t5])
                            
                        case _:
                            raise NotImplementedError(f'Triagulation not implemented for degree {degree} quadrature points.')

                # Create triangulation for tricontour plot
                triang  = tri.Triangulation(Y_gauss.reshape(-1),Z_gauss.reshape(-1),connectTrian)


                ## ------------------------ Tau_xy ------------------------ ##
                # Max and min stresses
                maxPlot     = np.max(tau_xy) 
                minPlot     = np.min(tau_xy) 
                # Set colobar ticks to 10 equaly spaced values between min and max magnitudes
                if 'ticks' not in cbarKwargs:
                    cbarKwargs['ticks'] = np.linspace(minPlot,maxPlot,11)
                
                plot    = ax[0].tricontourf(triang, tau_xy, cmap=cmap, **kwargs)

                # Add labels
                ax[0].set_title(r'$\tau_{xy}$'  , usetex=True, fontsize=24)
                ax[0].set_xlabel(r'$Y$'         , **fontkwargs)
                ax[0].set_ylabel(r'$Z$'         , **fontkwargs)
                ax[0].set_aspect('equal', adjustable='box')
                # ax[0].axis('tight')
                
                # Colorbar
                cbar    = plt.colorbar(plot, ax=ax[0], fraction=0.035, pad=0.04, use_gridspec=True, **cbarKwargs)
                cbar.set_label(labelBar , **fontkwargs)


                ## ------------------------ Tau_xz ------------------------ ##
                # Max and min stresses
                maxPlot     = np.max(tau_xz) 
                minPlot     = np.min(tau_xz) 
                
                plot    = ax[1].tricontourf(triang, tau_xz, cmap=cmap, 
                                            vmax=maxPlot, vmin=minPlot, 
                                            **kwargs)

                # Add labels
                ax[1].set_title(r'$\tau_{xz}$'  , usetex=True, fontsize=24)
                ax[1].set_xlabel(r'$Y$'         , **fontkwargs)
                ax[0].set_ylabel(r'$Z$'         , **fontkwargs)
                ax[1].set_aspect('equal', adjustable='box')
                # ax[1].axis('tight')
                
                # Colorbar
                cbar    = plt.colorbar(plot, ax=ax[1], 
                                       fraction=0.035, pad=0.04, 
                                       use_gridspec=True, 
                                       **cbarKwargs)
                cbar.set_label(labelBar , **fontkwargs)

            case _:
                raise Exception(f"Invalid mode: {mode}. Should be 'vector' or 'scalar'.")     # Invalid mode

        # If showMesh is requested
        if showMesh:
            for elem in self._elements:
                # Node coordinates
                nodesCoords = elem.getNodeCoords()
                
                # Determining the vertex order
                if self.nNodes== 9:
                    # Case Q9
                    inds = [0,4,1,5,2,6,3,7]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                elif self.nNodes== 6:
                    # Case T6
                    inds = [0,3,1,4,2,5]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                else:
                    # Case linear elements
                    Ys = np.array(nodesCoords)[:,1]
                    Zs = np.array(nodesCoords)[:,2]

                # Vertices
                verts = np.vstack([Ys,Zs]).T

                # Polygons
                poly    = Polygon(verts,fill=False,alpha=1,ec='k',lw=.3)
                if type(ax) == Axes:
                    # Case where the plot has only one axes (vector case)
                    ax.add_patch(poly)
                elif isinstance(ax, np.ndarray) and ax.ndim == 1 and isinstance(ax[0], Axes):
                    # Case where the plot has more than one axes (scalar case)
                    poly2    = Polygon(verts,fill=False,alpha=1,ec='k',lw=.3)
                    
                    ax[0].add_patch(poly)
                    ax[1].add_patch(poly2)
                    # {axes.add_patch(poly) for axes in ax}
                else:
                    raise Exception(f'Plot object type {type(ax)} not valid.')
        
        fig.tight_layout()

        return fig, ax
    ## CREATE GID POST-PROCESSING FILES
    def toGid(self, filename:str):
        '''
        Generates post-processing files for reading in GiD software.

        Parameters
        ----------
        filename : str
            The base filename for the .flavia.msh and .flavia.res files.
        '''

        nelem   = self.totalElements        # number of elements
        nnode   = self.nNodes           # number of nodes per element
        npnod   = self.totalDofs       # total number of nodes

        # Check mesh element type
        types = [elem.elemType for elem in self._elements]
        types = np.array(types)
        if np.all(types == 'Tri'):
            eletyp = 'Triangle'
        elif np.all(types == 'Quad'):
            eletyp = 'Quadrilateral'
        else:
            raise Exception('The mesh cotains more than one type of element, which is not supported by GiD')

        # msh and res filenames
        msh_file = filename + '.flavia.msh'
        res_file = filename + '.flavia.res'
               
        # Mesh file
        with open(msh_file,'w') as file:
            file.write('# \n')
            file.write(f'MESH dimension 2   Elemtype {eletyp}   Nnode {nnode:.0f} \n')
            file.write('coordinates \n')
            for i in range(npnod): #i = 1 : npnod
                file.write('{:6d} {:12.5f} {:12.5f} \n'.format(i+1,self._nodeCoords[i,1],self._nodeCoords[i,2]))
            file.write('end coordinates \n \n')

            file.write('elements \n')
            for i in range(nelem): # i=1:nelem
                # file.write('{:6d} {} \n'.format(i+1,self._elementsNodes[i].__str__()[1:-1]))
                print(f'{i+1:6d}',*np.array(self._elementsNodes[i][:])+1, file=file)
            file.write('end elements \n \n')
        # Console message
        print(f'File Saved: {msh_file}')

        # "result name" "analysis name" step_value result_type result_location "location name"
        # Results File (Warping)
        phi = self.displacements.reshape(-1)
        with open(res_file,'w') as file:
            file.write('Gid Post Results File 1.0 \n')
            file.write('# \n')
            file.write('Result "Warping" "analysis name" 1 Scalar OnNodes \n')
            file.write('ComponentNames "Warping Function" \n')
            file.write('Values \n')
            for i in range(self.totalDofs):# i=1:nelem
                file.write(f'{i+1:6d} {phi[i]:12.5f} \n')
            file.write('End Values \n')
        
        # Console message
        print(f'File Saved: {res_file}\n')
    ## START CROSS-SECTION CALCULATIONS
    def simulate(self, displayTimes=True):
        '''Start simulation: calculate area properties, solve warping, solve shear center'''
        ## ------------------------------------------------------------ ##
        ##                       SIMULATION START                       ##
        ## ------------------------------------------------------------ ##
        # Determine area properties
        tic = time.perf_counter()    
        self.__calcAreaProperties()
        toc = time.perf_counter()    
        self.times.append(f'Area properties integration time:    {toc-tic:.3f} seconds')     

        # Assemble global stiffness
        tic = time.perf_counter()    
        self.__globalStiffness(degree=self.intDegree)
        toc = time.perf_counter()    
        self.times.append(f'Global stiffness assembly time:      {toc-tic:.3f} seconds')
        
        
        # Node where phi* = 0
        prescribedDof = np.array([0])
        self.__boundaryConditions(prescribedDof)

        # Warping solution
        self.__solveWarping(degree=self.intDegree)

        # Determine mass properties
        self.__calcMassProperties()
        
        # Display times
        if displayTimes:
            [print(t) for t in self.times]
   