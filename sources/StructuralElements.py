##  ================================================================  ##
##                   CLASSES DE ELEMENTOS LINEARES                    ##
##  ================================================================  ##

from FEM_Secao2D import Elem1D

class Barra1D(Elem1D):
    '''Classe para elementos de barra 1D'''
    def __init__(self, nodes:tuple[Node,Node], secao:Secao|Secao2D, index:int, rot_x:float=0):
        # Herança da classe de elementos 1D
        Elem1D.__init__(self, nodes, index, rot_x)
        # Salvando dados da secao transversal
        self.secao = secao 

        A  = secao.propsArea['A_pond']  # Área ponderada
        E0 = secao.E[0]                 # Módulo de elasticidade de referência
        L  = self.getLength()           # Comprimento do elemento


        # Matriz de rigidez local
        self.stiffness = E0*A/L * np.array([[ 1,-1],
                                            [-1, 1]])

    def __str__(self):
        return 'Barra 1D'
    
    
    def getStiffness(self):
        '''## Retorna:
            stiff: matriz de rigidez local do elemento
            [
                u1,u2,              : deslocamentos de barra
            ]'''
        return self.stiffness
    
class Eixo1D(Elem1D):
    '''Classe para elementos de eixo 1D'''
    def __init__(self, nodes:tuple[Node,Node], secao:Secao|Secao2D, index:int, rot_x:float=0):
        # Herança da classe de elementos 1D
        Elem1D.__init__(self, nodes, index, rot_x)
        # Salvando dados da secao transversal
        self.secao = secao 

        Jphi = secao.propsArea['Jphi']      # Constante de torção
        L  = self.getLength()               # Comprimento do elemento
        G0 = secao.G[0]                     # Módulo de cisalhamento de referência

        # Matriz de rigidez local
        self.stiffness = G0*Jphi/L * np.array([[ 1,-1],
                                               [-1, 1]])
        
    def __str__(self):
        return 'Eixo 1D'
    
    def getStiffness(self):
        '''## Retorna:
            stiff: matriz de rigidez local do elemento
            [
                tx1, tx2            : deslocamentos de eixo
            ]'''
        return self.stiffness

class VigaAssim1D(Elem1D):
    '''Classe para elementos de viga de seção assimétrica 1D'''
    def __init__(self, nodes:tuple[Node,Node], secao:Secao|Secao2D, index:int, rot_x:float=0):
        # Herança da classe de elementos 1D
        Elem1D.__init__(self, nodes, index, rot_x)
        # Salvando dados da secao transversal
        self.secao = secao 

        Izz  = secao.propsArea['Izz_cent_pond']  # Momento de inercia ponderado e.r. ao centroide
        Iyy  = secao.propsArea['Iyy_cent_pond']  # Momento de inercia ponderado e.r. ao centroide
        Iyz  = secao.propsArea['Iyz_cent_pond']  # Momento de inercia ponderado e.r. ao centroide
        E0   = secao.E[0]                   # Módulo de elasticidade de ref.
        L    = self.getLength()             # Comprimento do elemento

        stif      = np.zeros([8,8])

        # indViga   = np.arange(4)
        # indViga2  = np.arange(4)+4


        # L1 , C1  = np.meshgrid(indViga  ,indViga  )
        # L2 , C2  = np.meshgrid(indViga2 ,indViga2 )


        stif[0:4,0:4]  = 2*E0*Izz/(L**2) *np.array([[  6/L,    3, -6/L,    3],
                                                    [    3,  2*L,   -3,    L],
                                                    [ -6/L,   -3,  6/L,   -3],
                                                    [    3,    L,   -3,  2*L]])
        
        stif[4:8,4:8]  = 2*E0*Iyy/(L**2) *np.array([[  6/L,   -3, -6/L,   -3],
                                                    [   -3,  2*L,    3,    L],
                                                    [ -6/L,    3,  6/L,    3],
                                                    [   -3,    L,    3,  2*L]])
        
        stif[0:4,4:8]  = 2*E0*Iyz/(L**2) *np.array([[  6/L,   -3, -6/L,   -3],
                                                    [   3,  -2*L,   -3,   -L],
                                                    [ -6/L,    3,  6/L,    3],
                                                    [    3,   -L,   -3, -2*L]])
        stif[4:8,0:4]  = stif[0:4,4:8].T

        self.stiffness   = stif

    def __str__(self):
        return 'Viga Assim. 1D'
    
    def getStiffness(self):
        '''## Retorna:
            stiff: matriz de rigidez local do elemento
            [
                v1, tz1, v2, tz2,   : deslocamentos de viga no plano xy local
                w1, ty1, w2, ty2,   : deslocamentos de viga no plano xz local
            ]'''
        return self.stiffness

class VigaAssim3D(Elem1D):
    '''Classe para elementos de viga de seção assimétrica 3D
    A rigidez local é dada no formato separado:
    
    u = {u1,u2, }
    
    '''
    def __init__(self, nodes:tuple[Node,Node], secao:Secao|Secao2D, index:int, rot_x:float=0):
        # Herança da classe de elementos 1D
        Elem1D.__init__(self, nodes, index, rot_x)
        # Salvando dados da secao transversal
        self.secao = secao 

        A    = secao.propsArea['A_pond']         # Área ponderada
        Izz  = secao.propsArea['Izz_pond_cent']  # Momento de inercia ponderado e.r. ao centroide
        Iyy  = secao.propsArea['Iyy_pond_cent']  # Momento de inercia ponderado e.r. ao centroide
        Iyz  = secao.propsArea['Iyz_pond_cent']  # Momento de inercia ponderado e.r. ao centroide
        Jphi = secao.propsArea['Jphi_pond_G']    # Constante de torção

        E0   = secao.E[0]                   # Módulo de elasticidade de ref.
        G0   = secao.G[0]                   # Módulo de cisalhamento de referência
        L    = self.getLength()             # Comprimento do elemento

        stif      = np.zeros([12,12])

        # indBarra  = np.arange(2)
        # indViga   = np.arange(4)+2
        # indViga2  = np.arange(4)+6
        # indEixo   = np.arange(2)+10


        # L0 , C0  = np.meshgrid(indBarra,indBarra)
        # L1 , C1  = np.meshgrid(indViga ,indViga )
        # L2 , C2  = np.meshgrid(indViga2 ,indViga2 )
        # L3 , C3  = np.meshgrid(indEixo ,indEixo )


        stif[0:2,0:2]      = E0*A/L*np.array([[ 1,-1],[-1, 1]])
        stif[2:6,2:6]      = 2*E0*Izz/(L**2) *np.array([[  6/L,    3, -6/L,    3],
                                                        [    3,  2*L,   -3,    L],
                                                        [ -6/L,   -3,  6/L,   -3],
                                                        [    3,    L,   -3,  2*L]])
        
        stif[2:6,6:10]     = 2*E0*Iyz/(L**2) *np.array([[  6/L,   -3, -6/L,   -3],
                                                        [   3,  -2*L,   -3,   -L],
                                                        [ -6/L,    3,  6/L,    3],
                                                        [    3,   -L,   -3, -2*L]])
        
        stif[6:10,6:10]    = 2*E0*Iyy/(L**2) *np.array([[  6/L,   -3, -6/L,   -3],
                                                        [   -3,  2*L,    3,    L],
                                                        [ -6/L,    3,  6/L,    3],
                                                        [   -3,    L,    3,  2*L]])
        stif[6:10,2:6]     = stif[2:6,6:10].T
        stif[10:12,10:12]  = G0*Jphi/L*np.array([[ 1,-1],[-1, 1]])

        self.stiffness   = stif
        
        ## --------------------------- Matriz de Massa me_ff --------------------------- ##
        yCT     = secao.propsArea['Y_CT']
        zCT     = secao.propsArea['Z_CT']
        rho_l   = secao.propsInercia['rho_l'        ]
        Qy      = secao.propsInercia['Qy_rho'       ]
        Qz      = secao.propsInercia['Qz_rho'       ]
        Iyy     = secao.propsInercia['Iyy_rho'      ]
        Izz     = secao.propsInercia['Izz_rho'      ]
        Iyz     = secao.propsInercia['Iyz_rho'      ]
        I_phi   = secao.propsInercia['I_phi_rho'    ]
        I_yphi  = secao.propsInercia['I_phi_y_rho'  ]
        I_zphi  = secao.propsInercia['I_phi_z_rho'  ]
        I_phi2  = secao.propsInercia['I_phi_phi_rho']
        # Argumentos:
        s = (    L,    yCT,    zCT,
             rho_l,     Qy,     Qz,
               Iyy,    Izz,    Iyz,
             I_phi, I_yphi, I_zphi, I_phi2)
        
        self.me_ff      = me_ff_func(*s)
        
        

    def __str__(self):
        return 'Viga Assim. 3D'
    def __repr__(self):
        stif    = self.getStiffness()
        # stif    = np.round(self.getStiffness(),3)
        k_uu    = stif[0:2,0:2]     # rigidez de barra
        k_viga  = stif[2:10,2:10] # rigidez de viga assim.
        k_tt    = stif[10:12,10:12]    # rigidez de eixo
        # k_vv = stif[2:6,2:6];   k_vw = stif[2:6,6:10]
        # k_ww = stif[6:10,6:10]; k_wv = stif[6:10,2:6]     
        return f'''\n({self._index}) Viga Assim. 3D:\n\t Rig. Barra: \n{k_uu},\n\t Rig. Viga:  \n{k_viga},\n\t Rig. Eixo:  \n{k_tt}'''

    def getStiffness(self):
        '''## Retorna:
            stiff: matriz de rigidez local do elemento
            [
                u1,u2,              : deslocamentos de barra
                v1, tz1, v2, tz2,   : deslocamentos de viga no plano xy local
                w1, ty1, w2, ty2,   : deslocamentos de viga no plano xz local
                tx1, tx2            : deslocamentos de eixo
            ]'''
        return self.stiffness