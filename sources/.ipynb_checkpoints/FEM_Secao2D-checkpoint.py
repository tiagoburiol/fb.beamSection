## ------------------------------ ATENÇÃO ------------------------------ ##
## É necessário instalar o pacote nptyping (>>> pip install nptyping)
## para os indicadores dos tipos do numpy
## ===================================================================== ##
import numpy as np
from itertools import product                       # Para loops aninhandos
from nptyping import NDArray, Int, Shape, Float     # Para anotações de tipos do numpy
import matplotlib.pyplot as plt

np.set_printoptions(precision=2,threshold=np.inf,linewidth=np.inf)


plt.rcParams.update({    
    # 'text.usetex'       :   False,
    # 'mathtext.fontset'  :   'custom',
    # "mathtext.bf"       :   'serif:bold',
    # 'mathtext.it'       :   'serif:italic',
    # 'mathtext.rm'       :   'serif',
    
    'mathtext.fontset'  :   'dejavuserif',
    
    
    # 'font.family'       :   'serif',
    # 'font.serif'        :   'dejavuserif',
    # 'font.serif'        :   'Times New Roman',
    'font.size'         :   10.0
})

##  ================================================================  ##
##                             Classe Nós                             ##
##  ================================================================  ##
class Node:
    '''Classe para nós'''
    ## Atributos da classe
    __nodeCount = 0   # Contador para o número atual de instâncias da classe
    
    ##  ---------------------  Construtor  ---------------------  ##
    def __init__(self, coords:np.ndarray, index:int = None):
        # Converte entrada para uma lista
        if type(coords) is not np.ndarray:
            self._coords = np.array(coords)
        else:
            self._coords = coords

        # Associa a contagem de nós ao índice do nó
        if index is None:
            self._index = Node.__nodeCount
        else:
            # Caso seja dado um índice
            self._index = index

        # Incrementa o contador de nós
        Node.__nodeCount += 1
    def __repr__(self):
        return str(np.round(self.getCoords(),3))

    ##  --------------------- Getters ---------------------  ##
    def getCoords(self)                     -> np.ndarray:
        return self._coords
    def getIndex(self)                      -> int:
        return self._index


##  ================================================================  ##
##                        CLASSES DE ELEMENTOS                        ##
##  ----------------------------------------------------------------  ##
class Elem:
    '''Classe genérica para elementos'''
    ##  ---------------------  Construtor  ---------------------  ##
    def __init__(self, nodes:tuple[Node,Node], index:int):
        self._nodes = nodes  # Nós do elemento (conectividade)
        self._index = index  # Índice do elemento 
    def __repr__(self):
        return 'Nós: ' + str(self.getNodes())

    ##  ---------------------  Getters  ---------------------  ##
    def getNodes(self)                      -> tuple[Node, Node]:
        '''Retorna uma tupla contendo os nós do elemento'''
        return self._nodes
    def getNodesIndexes(self)                      -> tuple[Node, Node]:
        '''Retorna uma tupla contendo os indices dos nos do elemento'''
        return [node._index for node in self._nodes]
    def getNodesCoords(self)                -> list:
        '''Retorna uma lista contendo as coordenadas dos nós do elemento'''
        return list(node.getCoords() for node in self.getNodes())
    def getNCs(self)                        -> list:
        '''Retorna uma lista contendo as coordenadas dos nós do elemento'''
        return list(node.getCoords() for node in self.getNodes())


##  ----------------------------------------------------------------  ##
##                            Elementos 1D                            ##
##  ----------------------------------------------------------------  ##
class Elem1D(Elem):
    '''Classe de elementos 1D'''
    ##  ---------------------  Construtor  ---------------------  ##
    def __init__(self, nodes:tuple[Node,Node], index:int, rot_x:float=0):
        ## Herança da classe elementos
        Elem.__init__(self, nodes, index)

        
        ## ___ Inicializações
        self._rot_x     = rot_x                     # [rad] Rotação em torno de x
        self._length    = np.nan
        self._localAxes = np.nan * np.ones([3,3])   # vetor de 3 vetores: [ex, ey, ez]
        self._midPoint  = np.nan


        ## ___ Var Locais
        p1, p2 = self.getNCs()  # obtendo as coordenadas dos nós
        #p1, p2 = nos[0], nos[1]

        ## ___ Cálculos
        self._length    = self.__calcLength__(p1,p2)     # calculando o comprimento 
        self._localAxes = self.__calcLocalAxes__(p1,p2)  # calculando os eixos locais
        self._midPoint  = self.__calcMidPoint__(p1,p2)   # calculando o ponto médio

        # Caso seja prescrita rotação em torno de x
        if rot_x != 0:
            self.setRotation(rot_x) 
    def __repr__(self):
        nodes  = self.getNodes()
        length = self.getLength()
        rot    = self.getRotation()
        return 'Nós: {}, Compr: {:.4e}, Rotacao x: {:.4e}'.format(nodes,length,rot)

    ##  ---------------------  Cálculos  ---------------------  ##
    def __calcLength__(self,p1,p2)          -> float:
        # Comprimento 
        L  = np.linalg.norm(p2-p1)
        return L
    def __calcLocalAxes__(self,p1,p2)       -> np.ndarray:
        # I = np.array([1,0,0])
        J = np.array([0,1,0])
        K = np.array([0,0,1])
        
        ## Direção x local
        i = (p2 - p1) / self._length
        # Cossenos diretores do eixo x local
        # cos_11 = np.sum(i*I)
        cos_12 = np.sum(i*J)
        cos_13 = np.sum(i*K)
        # Angulos diretores do eixo x local
        # lamb_11 = np.arccos(cos_11)
        lamb_12 = np.arccos(cos_12)
        lamb_13 = np.arccos(cos_13)

        ## Direção y local
        prod  = np.cross(K,i)
        norma = np.linalg.norm(prod)
        if norma < 1e-6 or np.isnan(norma):
            # Caso ey seja paralelo a (0,0,1)
            j = np.cross(J,i)/np.sin(lamb_12)
        else:
            j = prod/np.sin(lamb_13)

        ## Direção z local
        k = np.cross(i,j)
        return np.array([i,j,k])  
    def __calcMidPoint__(self,p1,p2)        -> np.ndarray:
        # Ponto Médio
        dp = p2-p1
        return p1 + np.array(dp)/2
        

    ##  ---------------------  Getters  ---------------------  ##
    def getLength(self)                     -> float:
        '''Retorna:
            L: comprimento do elemento'''
        return self._length
    def getLocalAxes(self)                  -> list[np.ndarray,np.ndarray,np.ndarray]:
        '''Retorna:
            [ex,ey,ez]: lista de vetores diretores dos eixos locais do elemento (ndarrays)'''
        return self._localAxes
    def getMidPoint(self)                   -> np.ndarray:
        '''Retorna:
            [mx,my,mx]: coordenadas do ponto médio do elemento (ndarray[3])'''
        return self._midPoint
    def getJacob(self)                      -> list[float,np.ndarray]:
        '''Calcula o jacobiano do elemento
        Retorna: 
            J: jacobiano (float)
            Jacob: matriz jacobiana (ndarray)'''
        J = self._length / 2
        Jacob = np.array(J)
        return J, Jacob
    def getRotation(self)                   -> float:
        '''Obtem o angulo de rotação do elemento em rad'''
        return self._rot_x

    ##  ---------------------  Setters  ---------------------  ##
    def setRotation(self, angulo):
        '''Rotaciona o elemento en torno do seu eixo x local'''
        C = np.cos(angulo)
        S = np.sin(angulo)
        # Rotação em torno de 'x'
        Rx = np.array([[1,0,0],[0,C,S],[0,-S,C]])
        # Atualiza os eixos locais
        self._localAxes = Rx @ self._localAxes
        self._rot_x     = angulo


##  ----------------------------------------------------------------  ##
##                  Elementos 2D (Seções Transversais)                ##
##  ----------------------------------------------------------------  ##
class Elem2D(Elem):
    '''Classe de elementos 2D para seções transversais'''
    
    ##  ---------------------  Construtor  ---------------------  ##
    def __init__(self, nodes:tuple[Node,...]):
        ## Herança da classe elementos
        Elem.__init__(self, nodes, None)

        self.nNos = len(nodes)

        # Tipo de elemento
        self.tipoElem       = ''
        if self.nNos in [3,6]:
            self.tipoElem   = 'Tri'
        elif self.nNos in [4,8,9]:
            self.tipoElem   = 'Quad' 
        else:
            raise Exception('Tipo de elemento não implementado')
    

    ##  JACOBIANO E DERIVADAS FÍSICAS
    def getJacob(self, zetas:np.ndarray)     -> list[float, np.ndarray, np.ndarray]:
        '''Calcula o jacobiano e derivadas físicas para elementos 2D
        Parâmetros:
            zetas: posição natural do ponto no interior do elemento 3d-array
        Retorna:
            J: jacobiano
            dNdy: derivada das func de interp e.r. a y
            dNdz: derivada das func de interp e.r. a z
            '''
        # Zetas (Necessário para elementos de ordem maior)
        ncs = np.array(self.getNCs())
        # print(f'{ncs = }')
        ys  = ncs[:,1]
        zs  = ncs[:,2]

        match self.nNos:
            ## T3  ___________________________________________
            case 3:
                # Coordenadas naturais do ponto
                zeta1, zeta2, zeta3 = zetas

                # Coordenadas dos nós do elemento
                y1,y2,y3 = ys
                z1,z2,z3 = zs

                Jy1 = y1
                Jy2 = y2
                Jy3 = y3

                Jz1 = z1
                Jz2 = z2
                Jz3 = z3

                # Matriz jacobiana
                matrizJacob = np.array([[   1,  1,  1],
                                         [Jy1,Jy2,Jy3],
                                         [Jz1,Jz2,Jz3]])

                # Jacobiano
                J = .5*np.linalg.det(matrizJacob) # Fellipa p. 364

                # Vetor de derivadas parciais
                derivs = 1/(2*J)*np.array([[z2-z3, z3-z1, z1-z2], 
                                           [y3-y2, y1-y3, y2-y1]])
                
                dNdy = derivs[0,:]
                dNdz = derivs[1,:]

            ## T6  __________________________________________
            case 6:
                # Coordenadas naturais do ponto
                zeta1, zeta2, zeta3 = zetas
                
                # Coordenadas dos nós do elemento
                y1,y2,y3,y4,y5,y6 = ys
                z1,z2,z3,z4,z5,z6 = zs

                # N1 = zeta1*(2*zeta1-1)
                # N2 = zeta2*(2*zeta2-1)
                # N3 = zeta3*(2*zeta3-1)
                # N4 = 4*zeta1*zeta2
                # N5 = 4*zeta2*zeta3
                # N6 = 4*zeta3*zeta1
                
                # Derivadas
                dNdzeta1 = np.array([4*zeta1-1,         0,         0, 4*zeta2,       0, 4*zeta3])
                dNdzeta2 = np.array([        0, 4*zeta2-1,         0, 4*zeta1, 4*zeta3,       0])
                dNdzeta3 = np.array([        0,         0, 4*zeta3-1,       0, 4*zeta2, 4*zeta1])

                # Componentes da matriz jacobiana
                Jy1 = np.dot(ys,dNdzeta1)
                Jy2 = np.dot(ys,dNdzeta2)
                Jy3 = np.dot(ys,dNdzeta3)

                Jz1 = np.dot(zs,dNdzeta1)
                Jz2 = np.dot(zs,dNdzeta2)
                Jz3 = np.dot(zs,dNdzeta3)

                # Matriz jacobiana
                matrizJacob = np.array([[  1,  1,  1],
                                        [Jy1,Jy2,Jy3],
                                        [Jz1,Jz2,Jz3]])
                
                # Jacobiano
                J = .5*np.linalg.det(matrizJacob)

                # Vetor de derivadas parciais
                P = 1/(2*J)*np.array([[Jz2-Jz3, Jy3-Jy2],
                                      [Jz3-Jz1, Jy1-Jy3],
                                      [Jz1-Jz2, Jy2-Jy1]])
                derivs = P.T@np.array([dNdzeta1, 
                                       dNdzeta2,
                                       dNdzeta3])
                
                dNdy = derivs[0,:]
                dNdz = derivs[1,:]

            ## Q4  __________________________________________
            case 4:
                # Coordenadas naturais do ponto
                eta, zeta = zetas

                # # Coordenadas dos nós do elemento
                # y1,y2,y3,y4, = ys

                # N1 = 1/4*(1-eta)*(1-zeta)
                # N2 = 1/4*(1+eta)*(1-zeta)
                # N3 = 1/4*(1+eta)*(1+zeta)
                # N4 = 1/4*(1-eta)*(1+zeta)
                
                # Derivadas
                dNdeta  = 1/4*np.array([zeta-1,   1-zeta, 1+zeta, -(1+zeta)])
                dNdzeta = 1/4*np.array([ eta-1, -(1+eta),  1+eta,     1-eta])

                # Componentes da matriz jacobiana
                Jy_eta  = np.dot(ys,dNdeta)
                Jy_zeta = np.dot(ys,dNdzeta)

                Jz_eta  = np.dot(zs,dNdeta)
                Jz_zeta = np.dot(zs,dNdzeta)

                # Matriz jacobiana
                matrizJacob = np.array([[ Jy_eta, Jz_eta],  # Componente Jyi multiplica a derivada dy
                                        [Jy_zeta,Jz_zeta]]) # Componente Jzi multiplica a derivada dz

                # Jacobiano
                J = np.linalg.det(matrizJacob)

                # Jacobiana inversa
                jacobInv = np.linalg.inv(matrizJacob)

                # Derivadas Físicas
                dNdksi = np.vstack((dNdeta, dNdzeta))#.reshape(2,-1)
                dNdy = jacobInv[0,:]@dNdksi
                dNdz = jacobInv[1,:]@dNdksi

            ## Q8  __________________________________________
            case 8:
                # Coordenadas naturais do ponto
                eta, zeta = zetas

                # # Coordenadas dos nós do elemento
                # y1,y2,y3,y4,y5,y6,y7,y8, = ys
                # z1,z2,z3,z4,z5,z6,z7,z8, = zs

                # N1 = 1/4*eta*zeta*(eta-1)*(zeta-1)
                # N2 = 1/4*eta*zeta*(eta+1)*(zeta-1)
                # N3 = 1/4*eta*zeta*(eta+1)*(zeta+1)
                # N4 = 1/4*eta*zeta*(eta-1)*(zeta+1)

                # N5 = 1/2*zeta*(1-eta**2)*(   zeta-1)
                # N6 = 1/2* eta*(   eta+1)*(1-zeta**2)
                # N7 = 1/2*zeta*(1-eta**2)*(   zeta+1)
                # N8 = 1/2* eta*(   eta-1)*(1-zeta**2)
                
                # Derivadas
                dNdeta  = 1/4*np.array([ zeta*(2*eta-1)*(zeta-1),   zeta*(2*eta+1)*(zeta-1),   zeta*(2*eta+1)*(zeta+1),  zeta*(2*eta-1)*(zeta+1),
                                             4*eta*zeta*(1-zeta),  -2*(2*eta+1)*(zeta**2-1),      -4*eta*zeta*(zeta+1),  2*(1-2*eta)*(zeta**2-1)])
                dNdzeta = 1/4*np.array([  eta*(eta-1)*(2*zeta-1),    eta*(eta+1)*(2*zeta-1),    eta*(eta+1)*(2*zeta+1),   eta*(eta-1)*(2*zeta+1),
                                         2*(1-2*zeta)*(eta**2-1),       -4*eta*zeta*(eta+1),  -2*(eta**2-1)*(2*zeta+1),       4*eta*zeta*(1-eta)])

                # Componentes da matriz jacobiana
                Jy_eta  = np.dot(ys,dNdeta)
                Jy_zeta = np.dot(ys,dNdzeta)

                Jz_eta  = np.dot(zs,dNdeta)
                Jz_zeta = np.dot(zs,dNdzeta)

                # Matriz Jacobiana
                matrizJacob = np.array([[ Jy_eta, Jz_eta],  # Componente Jyi multiplica a derivada d_eta
                                        [Jy_zeta,Jz_zeta]]) # Componente Jzi multiplica a derivada d_zeta
                
                # Jacobiano
                J = np.linalg.det(matrizJacob)
                
                # Jacobiana Inversa
                jacobInv = np.linalg.inv(matrizJacob)

                # Derivadas Físicas
                dNdksi = np.vstack((dNdeta, dNdzeta))#.reshape(2,-1)
                dNdy = jacobInv[0,:]@dNdksi
                dNdz = jacobInv[1,:]@dNdksi
            
            ## Q9  __________________________________________
            case 9:
                # Coordenadas naturais do ponto
                eta, zeta = zetas

                # # Coordenadas dos nós do elemento
                # y1,y2,y3,y4,y5,y6,y7,y8, = ys
                # z1,z2,z3,z4,z5,z6,z7,z8, = zs

                # N1 = 1/4*eta*zeta*(eta-1)*(zeta-1)
                # N2 = 1/4*eta*zeta*(eta+1)*(zeta-1)
                # N3 = 1/4*eta*zeta*(eta+1)*(zeta+1)
                # N4 = 1/4*eta*zeta*(eta-1)*(zeta+1)

                # N5 = 1/2*zeta*(1-eta**2)*(   zeta-1)
                # N6 = 1/2* eta*(   eta+1)*(1-zeta**2)
                # N7 = 1/2*zeta*(1-eta**2)*(   zeta+1)
                # N8 = 1/2* eta*(   eta-1)*(1-zeta**2)

                # N9 = (1-eta**2)*(1-zeta**2)
                
                # Derivadas
                dNdeta  = 1/4*np.array([ zeta*(2*eta-1)*(zeta-1),   zeta*(2*eta+1)*(zeta-1),   zeta*(2*eta+1)*(zeta+1),  zeta*(2*eta-1)*(zeta+1),
                                             4*eta*zeta*(1-zeta),  -2*(2*eta+1)*(zeta**2-1),      -4*eta*zeta*(zeta+1),  2*(1-2*eta)*(zeta**2-1),
                                             8*eta*(zeta**2-1)])
                dNdzeta = 1/4*np.array([  eta*(eta-1)*(2*zeta-1),    eta*(eta+1)*(2*zeta-1),    eta*(eta+1)*(2*zeta+1),   eta*(eta-1)*(2*zeta+1),
                                         2*(1-2*zeta)*(eta**2-1),       -4*eta*zeta*(eta+1),  -2*(eta**2-1)*(2*zeta+1),       4*eta*zeta*(1-eta),
                                         8*zeta*(eta**2-1)])

                # Componentes da matriz jacobiana
                Jy_eta  = np.dot(ys,dNdeta)
                Jy_zeta = np.dot(ys,dNdzeta)

                Jz_eta  = np.dot(zs,dNdeta)
                Jz_zeta = np.dot(zs,dNdzeta)

                # Matriz Jacobiana
                matrizJacob = np.array([[ Jy_eta, Jz_eta],  # Componente Jyi multiplica a derivada d_eta
                                        [Jy_zeta,Jz_zeta]]) # Componente Jzi multiplica a derivada d_zeta
                
                # Jacobiano
                J = np.linalg.det(matrizJacob)
                
                # Jacobiana Inversa
                jacobInv = np.linalg.inv(matrizJacob)

                # Derivadas Físicas
                dNdksi = np.vstack((dNdeta, dNdzeta))#.reshape(2,-1)
                dNdy = jacobInv[0,:]@dNdksi
                dNdz = jacobInv[1,:]@dNdksi

            case _:
                raise Exception('Número de nós inválido')
            
        return J, dNdy.reshape(1,-1), dNdz.reshape(1,-1)

 
    ##  FUNÇÕES DE INTERPOLAÇÃO
    def getShapeFun(self, zetas:np.ndarray)     -> np.ndarray:
        '''Determina o valor das funções de interpolação  em um ponto de um elemento 2D
        Parâmetros:
            zetas: posição natural do ponto no interior do elemento (3d-array)
        Retorna:
            N: vetor com o valor funções de interpolação avaliado em 'zetas'
            '''
        match self.nNos:
            ## T3 - LINEAR ------------------------------------------------------ ##
            case 3:
                # Coordenadas naturais do ponto
                zeta1,zeta2,zeta3 = zetas
                N1 = zeta1
                N2 = zeta2
                N3 = zeta3
                N  = np.array([N1,N2,N3])

            ## T6 - BILINEAR ---------------------------------------------------- ##
            case 6:
                # Coordenadas naturais do ponto
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
                # Coordenadas naturais do ponto
                eta,zeta = zetas
                N1 = 1/4*(1-eta)*(1-zeta)
                N2 = 1/4*(1+eta)*(1-zeta)
                N3 = 1/4*(1+eta)*(1+zeta)
                N4 = 1/4*(1-eta)*(1+zeta)
                N  = np.array([N1,N2,N3,N4])

            ## Q8 - QUADRÁTICO -------------------------------------------------- ##
            # case 8:
            #     # Coordenadas naturais do ponto
            #     eta,zeta = zetas
            #     N1 = 1/4*eta*zeta*(eta-1)*(zeta-1)
            #     N2 = 1/4*eta*zeta*(eta+1)*(zeta-1)
            #     N3 = 1/4*eta*zeta*(eta+1)*(zeta+1)
            #     N4 = 1/4*eta*zeta*(eta-1)*(zeta+1)
            #     N5 = 1/2*zeta*(1-eta**2)*(   zeta-1)
            #     N6 = 1/2* eta*(   eta+1)*(1-zeta**2)
            #     N7 = 1/2*zeta*(1-eta**2)*(   zeta+1)
            #     N8 = 1/2* eta*(   eta-1)*(1-zeta**2)
            #     N  = np.array([N1,N2,N3,N4,N5,N6,N7,N8])

            ## Q9 - BIQUADRÁTICO ------------------------------------------------ ##
            case 9:
                # Coordenadas naturais do ponto
                eta,zeta = zetas
                N1 = 1/4*eta*zeta*(eta-1)*(zeta-1)
                N2 = 1/4*eta*zeta*(eta+1)*(zeta-1)
                N3 = 1/4*eta*zeta*(eta+1)*(zeta+1)
                N4 = 1/4*eta*zeta*(eta-1)*(zeta+1)
                N5 = 1/2*zeta*(1-eta**2)*(   zeta-1)
                N6 = 1/2* eta*(   eta+1)*(1-zeta**2)
                N7 = 1/2*zeta*(1-eta**2)*(   zeta+1)
                N8 = 1/2* eta*(   eta-1)*(1-zeta**2)
                N9 = (1-eta**2)*(1-zeta**2)
                N  = np.array([N1,N2,N3,N4,N5,N6,N7,N8,N9])

            case _:
                raise Exception('Número de nós inválido')

        return N.reshape(1,-1)


    
    ##  QUADRATURAS 
    # Escolha da quadratura em função do tipo de elemento
    def quadratura(self, grau):
        # Tipo de elementos
        if self.tipoElem == 'Tri':
            return self.quadTri(grau)
        elif self.tipoElem == 'Quad':
            return self.quadQuad(grau)
        else:
            raise Exception('Tipo de elemento não implementado')
    
    # Triangulares
    def quadTri(self, grau):
        ''' Retorna os pontos (em coordenadas triangulares) e pesos para quadratura em triangulos
        ### SINTAXE:
            pnts, pesos = quadTri(grau)
        ### ENTRADAS:
            ``grau``: escalar
                Grau da quadratura (1 a 5)
        ### SAÍDAS:
            ``pnts``: array (n, 3)
            ``pesos``: array (1, n)
        '''

        # Pontos Internos
        match grau:
            case 1:
                pnts = [[1/3,1/3,1/3]]
                pesos = [1]

            case 2:
                w  = 1/3
                p1 = 2/3; p2 = 1/6
                pnts = [[p1,p2,p2],
                        [p2,p1,p2],
                        [p2,p2,p1]]
                pesos = [w,w,w]

            case 3 | 4:
                # Quadratura de grau 4, não há quadratura estável de grau 3
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
                pesos = [w1,w1,w1,w2,w2,w2]

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
                pesos = [w1,w2,w2,w2,w3,w3,w3]
            case _:
                raise NotImplementedError(f'{grau} não é um grau válido ou implementado')
        return np.array(pnts), np.array(pesos)
    
    # Quadrilaterais
    def quadQuad(self, grau):
        ''' Retorna os pontos (em coordenadas naturais) e pesos para quadratura 1D
        ### SINTAXE:
            pnts, pesos = quadQuad(grau)
        ### ENTRADAS:
            ``grau``: escalar
                Grau da quadratura (1 a )
        ### SAÍDAS:
            ``pnts``: array (1, n)
            ``pesos``: array (1, n)
        '''
        # Número de Pontos Necessários
        n_gp = np.ceil((grau+1)/2)

        # Pontos Internos
        # Dados de Fish (2007)
        match n_gp:
            case 1:
                pnts = [0]
                pesos = [2]

            case 2:
                p  = 1/np.sqrt(3)
                
                w  = 1

                pnts  = [-p, p]
                pesos = [w,w]

            case 3 :
                p1  = 0.7745966692
                p2  = 0

                w1  = 0.5555555556
                w2  = 0.8888888889
                
                pnts  = [-p1, p2, p1]
                pesos = [ w1, w2, w1]

            case 4:
                p1  = 0.8611363116
                p2  = 0.3399810436

                w1  = 0.3478548451
                w2  = 0.6521451549

                pnts  = [-p1, -p2, p2, p1]
                pesos = [ w1,  w2 ,w2, w1]
                
            case 5:
                p1  = 0.9061798459
                p2  = 0.5384693101
                p3  = 0

                w1  = 0.2369268851
                w2  = 0.4786286705
                w3  = 0.5688888889

                pnts  = [-p1, -p2, p3, p2, p1]
                pesos = [ w1,  w2, w3, w2, w1]

            case 6:
                p1  = 0.9324695142
                p2  = 0.6612093865
                p3  = 0.2386191861

                w1  = 0.1713244924
                w2  = 0.3607615730
                w3  = 0.4679139346

                pnts  = [-p1, -p2, -p3, p3, p2, p1]
                pesos = [ w1,  w2,  w3, w3, w2, w1]

            case _:
                raise NotImplementedError(f'{grau} não é um grau válido ou implementado')
        return np.array(pnts), np.array(pesos)
    
    
    ##  MATRIZ DE RIGIDEZ LOCAL  
    def getKlocal(self,grau:int)                -> np.ndarray:
        '''Calcula a matriz de rigidez local do elemento
        ##  Parâmetros:
            grau: grau da quadratura para a integração numérica
        ## Retorna:
            soma: matriz de rigidez local'''

        zetas, pesos = self.quadratura(grau)

        soma = 0
        # Quadratura para o elemento; Varre os pontos gaussianos
        if self.tipoElem == 'Tri':
            for i, coords in enumerate(zetas):
                #     i: indice
                # coord: coordenadas naturais do ponto (zeta1,2,3)

                J, By, Bz = self.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss

                # Termo Integrando
                BB = By.T @ By + Bz.T @ Bz
                #BB = By[...,None] @ By + Bz[...,None] @ Bz


                soma += pesos[i]*J*BB

        elif self.tipoElem == 'Quad':
            for i, zeta_i in enumerate(zetas):
                for j, zeta_j in enumerate(zetas):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,2,3)
                    coords = [zeta_i, zeta_j]
                    J, By, Bz = self.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss

                    # Termo Integrando
                    BB = By.T @ By + Bz.T @ Bz
                    #BB = By[...,None] @ By + Bz[...,None] @ Bz


                    soma += pesos[i]*pesos[j]*J*BB

        return soma


    ##  VETOR DE FORÇAS LOCAL
    def getFlocal(self,grau:int)                -> np.ndarray: 
        '''Calcula o vetor de forças local do elemento
        ##  Parâmetros:
            grau: grau da quadratura para a integração numérica
        ## Retorna:
            soma: vetor de forças local'''

        pontos, pesos = self.quadratura(grau)
        ncs = np.array(self.getNCs())
        ys  = ncs[:,1]
        zs  = ncs[:,2]

        soma = 0
        # Quadratura para o elemento; Varre os pontos gaussianos

        if self.tipoElem == 'Tri':
            for i, coords in enumerate(pontos):
                #     i: indice
                # coord: coordenadas naturais do ponto (zeta1,2,3)

                # Termo Integrando
                N           = self.getShapeFun(coords)
                J, By, Bz   = self.getJacob(coords) 
                # print(f'{ys.shape = }')
                # print(f'{zs.shape = }')
                # print(f'{N.shape = }')
                # print(f'{By.shape = }')
                # print(f'{Bz.shape = }')
                # F           = By.T @ zs.reshape(1,-1) @ N.T - Bz.T @ ys.reshape(1,-1) @ N.T
                F           = By.T @ N @ zs.reshape(-1,1) - Bz.T @ N @ ys.reshape(-1,1)

                soma += pesos[i]*J*F

        elif self.tipoElem == 'Quad':
            for j, zeta_j in enumerate(pontos):
                for i, zeta_i in enumerate(pontos):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,2,3)
                    coords = [zeta_i, zeta_j]
                    J, By, Bz = self.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss

                    # Termo Integrando
                    N = self.getShapeFun(coords)
                    # F = By.T @ zs.reshape(1,-1) @ N.T - Bz.T @ ys.reshape(1,-1) @ N.T
                    F = By.T @ N @ zs.reshape(-1,1) - Bz.T @ N @ ys.reshape(-1,1)

                    soma += pesos[i]*pesos[j]*J*F

        return soma


##  ================================================================  ##
##                         CLASSES DE MALHAS                          ##
##  ----------------------------------------------------------------  ##
class Malha:
    '''Classe genérica para malhas'''
    def __init__(self, tipoElem: Elem1D | Elem2D, elems:list[Elem]=None, coordenadas:list[list[float,float,float]]=None, conectividade:list[list[int,int]]=None):
        # # Inicialização dos Elementos

        if elems is not None:
            #### Caso seja fornecida uma lista de elementos
            self._elems = list(elems)

            elementsNodes = []
            coords = []
            for elem in elems:
                elementsNodes.append([node._index for node in elem._nodes]) # Conectividades
            self._elementsNodes = elementsNodes
            # print(self._elementsNodes)
            
            # Lista ordenada dos nós
            elemNodes = np.array(elementsNodes)
            elemNodes = set(elemNodes.flatten())
            nodeIndexes = list(elemNodes)
            indices  = []
            nodes    = []
            for i in nodeIndexes:
                for elem in elems:
                    for node in elem._nodes:
                        if node._index == i and node._index not in indices:
                            nodes.append(node)
                            coords.append(node._coords)
                            indices.append(node._index)

            self._nodesCoords = np.array(coords)
            self._nodes       = nodes
            # print(self._nodesCoords)


        elif (coordenadas is not None) and (conectividade is not None):
            #### Caso sejam fornecidos vetores de coordenadas e conectividade
            # Criação da lista de instancias dos nós
            nos   = [Node(nodeCoords,indNode) for indNode,nodeCoords in enumerate(coordenadas)]

            # Varrendo o vetor de elementos
            elems = []
            for elemNodes in conectividade:
                nos_elem = []
                nos_malha = []
                for node in elemNodes:
                    # Buscando instancia dos nós do elemento e pondo em uma lista
                    nos_elem.append(nos[node])
                    # Criando uma lista de nós da malha
                    if nos[node] not in nos_malha:
                        nos_malha.append(nos[node])
                # Passando lista de instancia de nós para um elemento
                if tipoElem == Elem1D:
                    elems.append(Elem1D(nos_elem))
                elif tipoElem == Elem2D:
                    elems.append(Elem2D(nos_elem))
                    
            self._elems = elems
            self._nodes = nos_elem
            
            self._elementsNodes = conectividade
            self._nodesCoords   = coordenadas
    def __repr__(self):
        lista = [f'Elemento: {elem._index}, {elem}' for elem in self._elems]
        return str(lista)

    ##  ---------------------  Getters  ---------------------  ##
    def getElems(self):
        return self._elems

    ##  ---------------------  Métodos  ---------------------  ##
    def add(self, elem:Elem):
        '''Adicionar um elemento à malha'''
        self._elems.append(elem)
    def delete(self, indice:int):
        '''remover um elemento da malha'''
        self._elems.remove(self._elems[indice])
    def indice(self, elem:Elem):
        '''Busca o índice de um elemento da malha'''
        self._elems.index(elem)
    def dataframe(self, precisao:int = 4):
        '''Mostra a malha como um dataframe'''
        import pandas as pd

        col0 = 'Elemento'
        # Key para os indides dos nó
        col1 = 'Nós'
        # Key para coordenadas
        col2 = 'Coordenadas'
        col3 = 'Rotação x [rad]'
        col4 = 'i_e'
        col5 = 'j_e'
        col6 = 'k_e'
        dic  = {col0:[], col1:[], col2:[], col3:[], col4:[], col5:[], col6:[]}
        dfIndex = []
        for elem in self._elems: 
            indElem = elem._index
            nos     = elem._nodes
            indNos  = [node._index  for node in nos]
            coords  = [np.round(node._coords, precisao) for node in nos]
            rot     = np.round(elem._rot_x, precisao)
            i,j,k   = np.round(elem.getLocalAxes(),3)
            # classe  = elem.__str__()
            
            # Montando dicionario
            dic[col0].append(indElem)
            dic[col1].append(indNos)
            dic[col2].append(coords)
            dic[col3].append(rot)
            dic[col4].append(i)
            dic[col5].append(j)
            dic[col6].append(k)
            # dic[col4].append(classe)
            dfIndex.append('')              # para inibir a coluna de índices
        df = pd.DataFrame(data=dic,index=dfIndex)
        # df.index.name = ''
        return df
    

##  ----------------------------------------------------------------  ##
##                               Malha 1D                             ##
##  ----------------------------------------------------------------  ##
class Malha1D(Malha):
    '''Classe para malhas reticuladas, ou seja, de elementos retilíneos no espaço.
    Recebe uma lista de elementos'''
    def __init__(self, elems:list[Elem]=None, coordenadas:list[list[float,float,float]]=None, conectividade:list[list[int,int]]=None):
        ## Herança da classe de malhas
        Malha.__init__(self,Elem1D,elems,coordenadas,conectividade)            
        
    ##  ---------------------  Plot da malha  ---------------------  ##
    def plotarMalha(self, mostrarLabelElems:bool=True,mostrarLabelNos:bool=True,mostrarEixosLocais:bool=True,mostrarSecoes:bool=False,fatorEscSec:float=1):
        '''Plotar a malha reticulada'''
        import matplotlib.pyplot as plt

        # Criação dos axes
        elevacao    = 30
        azimute     = 45
        ax          = plt.axes(projection='3d',azim=azimute, elev=elevacao)
        

        # Loop para os elementos
        for indElem, elem in enumerate(self._elems): 
            ## ----------------------------------------------------------------- ##
            ##                       MOSTRAR OS ELEMENTOS                        ##
            ## ----------------------------------------------------------------- ##
            # Obtendo as coordenadas dos nós e convertendo em array do numpy
            nos = np.array( elem.getNCs() ) 
            # nos = list([X1,Y1,Z1],[X2,Y2,Z3])
            X = nos[:,0]
            Y = nos[:,1]
            Z = nos[:,2]
            # Plot dos elementos
            ax.plot3D(X,Y,Z,marker='o',lw=1,c='k',mec='#04a301',mfc='w',markersize=5,label=None)

            ## ----------------------------------------------------------------- ##
            ##                MOSTRAR AS ETIQUETAS DOS ELEMENTOS                 ##
            ## ----------------------------------------------------------------- ##
            # Obtem o ponto médio
            [Xm, Ym, Zm] = elem._midPoint

            # Etiqueta dos Elementos
            if mostrarLabelElems:
                label = '  ('+str(indElem)+')'    
                ax.text(Xm,Ym,Zm,label,zdir=None,c='k',va='top',ha='left',fontweight='bold',fontstyle='italic')
            
            # Etiqueta dos Nós
            if mostrarLabelNos:
                # Buscando nós de cada elemento 
                nos     = elem._nodes
                indNos  = [node._index  for node in nos]    # indices dos nós
                coords  = [node._coords for node in nos]    # coordenadas dos nós
                # Buscando o índice e coordenadas de cada nó
                for i,ind in enumerate(indNos):
                    label = str(ind) + '  '                 # etiqueta para o nó
                    X,Y,Z = coords[i]                       # posição do nó
                    ax.text(X,Y,Z,label,zdir=None,c='#04a301',va='top',ha='right',fontweight='bold')
            
            ## ----------------------------------------------------------------- ##
            ##                      MOSTRAR OS EIXOS LOCAIS                      ##
            ## ----------------------------------------------------------------- ##
            if mostrarEixosLocais:
                ### Loop para cada cosseno diretor
                for i in range(3):
                    # cDir = self.cosDir[nElem][i] #+ Xm
                    # pm   = [Xm,Ym,Zm] # Ponto médio do elemento
                    cDir = elem._localAxes[i]
                    pm   = elem._midPoint
                    p2   = pm + .5*cDir
                    
                    # Eixos locais
                    x = [pm[0], p2[0]] 
                    y = [pm[1], p2[1]]
                    z = [pm[2], p2[2]]
                    
                    # Determinação da cor da linha
                    if   i == 0:
                        cor = 'r' # x-local
                    elif i == 1:
                        cor = '#ffd20a' # y-local
                    else:
                        cor = 'b' # z-local

                        
                    # Plot das linhas indicadoras do eixos locais
                    ax.plot3D(x,y,z,marker=None,lw=2,c=cor)
            
            ## ----------------------------------------------------------------- ##
            ##                        MOSTRAR AS SEÇÕES                          ##
            ## ----------------------------------------------------------------- ##
            if mostrarSecoes:
                y_cg        = elem.secao.propsArea['Y_CG']  # Posição do centroide
                z_cg        = elem.secao.propsArea['Z_CG']
                elemNodes   = elem.secao._elementsNodes     # Conectividade da malha da seção
                nNos        = elem.secao.nNodes
                # Corrigindo a origem do sistema para o cg e ponto médio do elemento
                nos2D            = elem.secao._nodesCoords
                nodesCoords      = np.zeros(nos2D.shape)    # Coordenadas da malha da seção

                # Corrigindo a origem do sistema para o centroide e aplicando o fator de escala
                nodesCoords[:,1] = fatorEscSec*(nos2D[:,1] - y_cg)  
                nodesCoords[:,2] = fatorEscSec*(nos2D[:,2] - z_cg)
                
                # print('elemNodes: \n', elemNodes)
                # Determinando a ordem dos vertices
                match nNos:
                    case 9:
                        inds = [0,4,1,5,2,6,3,7]
                    case 6:
                        inds = [0,3,1,4,2,5]
                    case 4:
                        inds = [0,1,2,3]
                    case 3:
                        inds = [0,1,2]
                # MCD do elemento 1D
                Ae         = elem.getLocalAxes()            
                for elem2D in elemNodes:                        # Percorre os elementos da seção transv
                    # print('elem2D: \n',elem2D)
                    
                    if type(elem2D) != np.ndarray :
                        elem2D = np.array(elem2D)
                    # Reordenando os nós para manter geometria coerente
                    elem2D = elem2D[inds]

                    for k, node in enumerate(elem2D):           # Percorre os nós do elemento
                        # Indices do elemento 2D
                        IND0 = elem2D[k]
                        try:
                            IND1 = elem2D[k+1]
                        except IndexError:
                            # Se a pos0 for a última, a proxima será a inicial
                            IND1 = elem2D[0]

                        # Transformação os eixos do corpo
                        X0,Y0,Z0 = Ae.T@nodesCoords[IND0,:] # nó 0
                        X1,Y1,Z1 = Ae.T@nodesCoords[IND1,:] # nó 1

                        Xs = np.array([X0,X1])
                        Ys = np.array([Y0,Y1])
                        Zs = np.array([Z0,Z1])

                        # Determinando a ordem dos vertices
                        match nNos:
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


        # ax.set_proj_type('ortho')
        ax.set_xlabel(r'$X$',fontsize=18,usetex=True);  # ax.set_xlim([np.min(X),np.max(X)])
        ax.set_ylabel(r'$Y$',fontsize=18,usetex=True);  # ax.set_ylim([np.min(Y),np.max(Y)])
        ax.set_zlabel(r'$Z$',fontsize=18,usetex=True);  # ax.set_zlim([np.min(Z),np.max(Z)])

        # ax.set_aspect('equal')
        # ax.set_box_aspect([1,1,1])
        if mostrarEixosLocais:
            ax.legend(['Malha','x-local','y-local','z-local'])
        else:
            ax.legend(['Malha'])

        return ax
    


##  ================================================================  ##
##                          CLASSES DE SEÇÕES                         ##
##  ----------------------------------------------------------------  ##
class Secao:
    '''Classe para criação das propriedades de área de uma seção genérica.
    Parâmetros:
        propsArea: {key:value} com pares de propriedade e valor
        {
            keys:       'A': Área, 
                       'Qy': Mom. Estático de Área e.t. de y-local, 
                       'Qz': Mom. Estático de Área e.t. de z-local, 
                      'Iyy': Mom. de Inércia de Área e.t. de y-local, 
                      'Izz': Mom. de Inércia de Área e.t. de z-local, 
                      'Iyz': Produto de Inércia de Área no plano yz-local, 
                     'Jphi': Constante de Torção 
            values: valor da propriedade (float)
        }
    Exemplo:
        propsArea = dict(A=10,Izz=30,...)'''
    def __init__(self, propsArea:dict[str,float], E, G):
        self.propsArea = propsArea
        self.E = E
        self.G = G

    ##  ---------------------  Getters  ---------------------  ##
    def getPropsArea(self)                  -> dict[str,float]:
        '''Obtem as propriedades de área da seção'''
        return self.propsArea
    
    ##  ---------------------  Setters  ---------------------  ##
    def setPropsArea(self, propsArea:dict[str,float]):
        '''Atribui as propriedades de área à seção'''
        self.propsArea = propsArea


##  ----------------------------------------------------------------  ##
##                              Seção 2D                              ##
##  ----------------------------------------------------------------  ##
class Secao2D(Malha):
    '''Subclasse para criação de uma seção transversal a partir de 
    uma malha de elementos 2D
    ## Parâmetros
        elems: lista de objetos Elem2D (opcional)
        coords: lista de coordenadas dos nós
        conect: lista de nós dos elementos
        young: lista dos módulos de elasticidade dos elementos
    '''
    def __init__(self, 
                 elems:         list[Elem2D]        = None, 
                 coordenadas:   list[list[float]]   = None, 
                 conectividade: list[list[int]]     = None,
                 young:         list[list[int]]     = None,
                 nu:            list[list[int]]     = None,
                 rho:           list[list[int]]     = None,
                 grauInteg:     int                 = None,
                 ):

        ## Herança da classe de malha 2D
        Malha.__init__(self,Elem2D,elems,coordenadas,conectividade)
        
        # Propriedades dos materiais
        self.E              = young
        self.nu             = nu
        self.G              = self.E / (2*(1 + self.nu))
        self.rho            = rho
        self.nElems         = len(conectividade)
        self.nNodes         = len(conectividade[0])
        self.totalNodes     = len(coordenadas)
        self.grauInteg      = grauInteg
        
        # Propriedades para o PVC do empenamento
        self.activeDof      = None
        self.prescribedDof  = None
        self.stiffness      = np.zeros([self.totalNodes,self.totalNodes])
        self.displacements  = np.zeros([self.totalNodes,1])
        self.forces         = np.zeros([self.totalNodes,1])
        self.phiStar        = np.zeros([self.totalNodes,1])
        self.propsArea      = {}
        self.propsInercia   = {}

        # Tipo de elemento
        self.tipoElem       = ''
        if self.nNodes in [3,6]:
            self.tipoElem   = 'Tri'
        elif self.nNodes in [4,8,9]:
            self.tipoElem   = 'Quad' 
        else:
            raise Exception('Tipo de elemento não implementado')

        # Ordem de integração
        if (self.nNodes in [3,4]) and (self.grauInteg is None):
            self.grauInteg = 2      # Grau 2 para elementos lineares
        else:
            self.grauInteg = 4      # Grau 4 para elementos quadráticos
        
        # Determina as prop area
        self.__calcpropsArea()

        # Montar rigidez global
        self.__rigidezGlobal(grau=self.grauInteg)
        # self.__rigidezGlobal(grau=4)

        # Nó onde phi* = 0
        prescribedDof = np.array([0])
        self.__condContorno(prescribedDof)

        # Solução do empenamento
        self.__solucao(grau=self.grauInteg)
        # self.__solucao(grau=4)

        # Determina as prop de inercia
        self.__calcpropsInercia()
        
        
    ##  ---------------------  Getters  ---------------------  ##
    def getPropsArea(self)                  -> dict[str,float]:
        '''Obtem as propriedades de área da seção'''
        return self.propsArea

        ##  -----------------  Métodos Públicos -----------------  ## 
    def plotPhi(self, niveis = 20, figsize=(6,4), mostrarMalha=False, mostrarCG=False, num=1, taxa=None, unit='l.u.', cbarKwargs=dict({}), **kwargs):
        '''Plota a distribuição da função empenamento na seção transversal
        ## Parâmetros:
            níveis: (int) número de divisões da curva de níveis
            mostrarMalha: (bool) desenha a malha utilizada   
            num: número da figura
            taxa: taxa de torção, caso seja fornecido mostra o grafico do deslocamento 
            unit: unidade de comprimento utilizada
        ## Retorna:
            ax: plot da função empenamento na seção'''
        import matplotlib.tri as tri 
        import matplotlib.pyplot as plt

        phi = np.array(self.displacements).reshape(-1)
        if taxa != None:
            phi = taxa*phi
            label = f'Warping Displacement $\\alpha\\varphi^h$ [{unit}]'+'\n'+ f'$\\alpha$ = {taxa:.3e} [rad/{unit}]'
        else:
            label = f'Warping Function $\\varphi^h$ [{unit}^2]'
        Y   = np.array(self._nodesCoords)[:,1].reshape(-1)
        Z   = np.array(self._nodesCoords)[:,2].reshape(-1)

        fig, ax = plt.subplots(num=num, figsize=figsize)
        # Determinando a triangulação adequadas para o pósprocessamento 
        match self.nNodes:
            case 3:
                t0 = np.array(self._elementsNodes)[:,[0,1,2]]
                conect = t0
            case 6:
                t0 = np.array(self._elementsNodes)[:,[5,4,2]]
                t1 = np.array(self._elementsNodes)[:,[0,4,5]]
                t2 = np.array(self._elementsNodes)[:,[0,3,4]]
                t3 = np.array(self._elementsNodes)[:,[3,1,4]]
                conect = np.vstack([t0,t1,t2,t3])
            case 4:
                t0 = np.array(self._elementsNodes)[:,[0,1,2]]
                t1 = np.array(self._elementsNodes)[:,[2,3,0]]
                conect = np.vstack([t0,t1])
            case 9:
                t0 = np.array(self._elementsNodes)[:,[7,8,3]]
                t1 = np.array(self._elementsNodes)[:,[0,8,7]]
                t2 = np.array(self._elementsNodes)[:,[0,4,8]]
                t3 = np.array(self._elementsNodes)[:,[4,1,8]]
                t4 = np.array(self._elementsNodes)[:,[1,5,8]]
                t5 = np.array(self._elementsNodes)[:,[8,5,2]]
                t6 = np.array(self._elementsNodes)[:,[8,2,6]]
                t7 = np.array(self._elementsNodes)[:,[8,6,3]]
                conect = np.vstack([t0,t1,t2,t3,t4,t5,t6,t7])


        triang  = tri.Triangulation(Y,Z,conect)
        contour = ax.tricontourf(triang, phi, cmap='jet',levels=niveis, **kwargs)
        cbar    = plt.colorbar(contour,  **cbarKwargs)
        # cbar.set_label(r'Warping Function $\varphi$', fontname='Times New Roman')
        cbar.set_label(label, fontname='DejaVu Serif',fontsize=10)
        # cbar.set_label(r'Função Empenamento $\varphi$ [m²/rad]', fontname='Times New Roman',fontsize=16)

        if mostrarMalha:
            from matplotlib.patches import Polygon            
            # corPreench  = '#59c1f9'

            for elem in self._elems:
                # Coordenadas dos nos
                nodesCoords = elem.getNCs()
                
                # Determinando a ordem dos vertices
                if self.nNodes == 9:
                    # Caso Q9
                    inds = [0,4,1,5,2,6,3,7]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                elif self.nNodes == 6:
                    # Caso T6
                    inds = [0,3,1,4,2,5]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                else:
                    # Caso elementos lineares
                    Ys = np.array(nodesCoords)[:,1]
                    Zs = np.array(nodesCoords)[:,2]

                # Vertices
                verts = np.vstack([Ys,Zs]).T

                # poly        = Polygon(verts,fc=corPreench,alpha=.2,ec='k')
                poly        = Polygon(verts,fill=False,alpha=1,ec='k',lw=.3)
                ax.add_patch(poly)
                

        # Mostra o Centro de Torção
        ax.scatter(self.propsArea['Y_CT'], self.propsArea['Z_CT'], marker='^', edgecolors='r',facecolors='w', s=100, label='CT')
        # Mostra o Centroide
        if mostrarCG:
            ax.scatter(self.propsArea['Y_CG_pond'], self.propsArea['Z_CG_pond'], marker='o', edgecolors='b',facecolors='w', s=100, label='CG')

        # ax.set_title(f'Seção Transversal',fontsize=14,fontweight='bold')
        # ax.set_xlim([np.min(Y),np.max(Y)])
        # ax.set_ylim([np.min(Z),np.max(Z)])
        # ax.axis('square')
        ax.set_aspect('equal', adjustable='box')
        # ax.set_aspect('equal')
        
        # ax.set_xlabel(f'$Y$ [{unit}]',usetex=True,fontsize=20)
        # ax.set_ylabel(f'$Z$ [{unit}]',usetex=True,fontsize=20)
        # ax.tick_params(axis='y',labelsize=12)
        # ax.tick_params(axis='x',labelsize=12)
        
        # plt.legend(loc='center right')
        fig.tight_layout()
        return fig, ax
    def plotSec(self, corPreench:str='#59c1f9'):
        '''Plota a seção transversal'''
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon

        ax2         = plt.axes()
        for elem in self._elems:
            # Coordenadas dos nos
            nodesCoords = elem.getNCs()
            Ys = np.array(nodesCoords)[:,1]
            Zs = np.array(nodesCoords)[:,2]
            verts = np.vstack([Ys,Zs]).T

            # Vertices
            poly        = Polygon(verts,fc=corPreench,alpha=.2,ec='k')
            ax2.add_patch(poly)
            
        Ys = np.array(self._nodesCoords)[:,1]
        Zs = np.array(self._nodesCoords)[:,2]

        ax2.set_aspect('equal')
        ax2.set_title(f'Seção Transversal',fontsize=14,fontweight='bold')
        ax2.set_xlim([np.min(Ys),np.max(Ys)])
        ax2.set_ylim([np.min(Zs),np.max(Zs)])
        # plt.tight_layout()
        return 
    def toGid(self, filename:str):
        '''Gera arquivos de pós-processamento para leitura no GiD.'''

        nelem  = self.nElems        # número de elementos
        nnode  = self.nNodes        # número de nós por elemento
        npnod  = self.totalNodes    # numeros total de nós

        if self.tipoElem == 'Tri':
            eletyp = 'Triangle'
        else:
            eletyp = 'Quadrilateral'

        msh_file = filename + '.flavia.msh'
        res_file = filename + '.flavia.res'

        # Arquivo da malha
        with open(msh_file,'w') as file:
            file.write('# \n')
            file.write(f'MESH dimension 2   Elemtype {eletyp}   Nnode {nnode:.0f} \n')
            file.write('coordinates \n')
            for i in range(npnod): #i = 1 : npnod
                file.write('{:6d} {:12.5f} {:12.5f} \n'.format(i+1,self._nodesCoords[i,1],self._nodesCoords[i,2]))
            file.write('end coordinates \n \n')

            file.write('elements \n')
            for i in range(nelem): # i=1:nelem
                # file.write('{:6d} {} \n'.format(i+1,self._elementsNodes[i].__str__()[1:-1]))
                print(f'{i+1:6d}',*np.array(self._elementsNodes[i][:])+1, file=file)
            file.write('end elements \n \n')
        # Mensagem no console
        print(f'Arquivo Salvo: {msh_file}')

        # "result name" "analysis name" step_value result_type result_location "location name"
        # Arquivo dos Resultados (Empenamento)
        phi = self.displacements.reshape(-1)
        with open(res_file,'w') as file:
            file.write('Gid Post Results File 1.0 \n')
            file.write('# \n')
            file.write('Result "Warping" "analysis name" 1 Scalar OnNodes \n')
            file.write('ComponentNames "Warping Function" \n')
            file.write('Values \n')
            for i in range(self.totalNodes):# i=1:nelem
                file.write(f'{i+1:6d} {phi[i]:12.5f} \n')
            file.write('End Values \n')
        
        # Mensagem no console
        print(f'Arquivo Salvo: {res_file}\n')
    def plotTensoesCis(self, grau:int = 2, modo:str = 'vetorial', mostrarMalha=True, taxa=None, probeElems=None, unid='Mpa',
                       vetoresUnit:bool = True, POS:str = 'nos', passoVetores=1, figsize=(6,4), mises=False, 
                       cbarKwargs=dict(),**kwargs):
        '''Plota as tensões cisalhantes normalizadas por G*alpha para 
            visualização do fluxo de cisalhamento
        kwargs:
            modo : str = 'vetorial'
                    'vetorial : mostra o campo vetorial das tensões
                    'escalar' : mostra o campo escalar de cada componente
            mostrarMalha: bool = True
                Mostrar o a malha no plot
            taxa : float = None
                Taxa de torção. Se não fornecida as tensões serão dadas por unidade de alpha.
            unid : str = 'Mpa'
                Unidade de medida da tensão
            vetoresUnit : bool = True
                    utilizar vetores normalizados no modo vetorial
            POS  : bool = True
                    utilizar vetores normalizados no modo vetorial

            grau : str = 'nos'
                    posição do pontos onde a tensão é calculada
                    'nos'  : nos nós de cada elemento
                    'gauss': pontos de Gauss
            **kwargs : kwargs para a função de plot
                
        '''
        
        pontos, _   = self._elems[0].quadratura(grau)   # Calculando nos pontos gaussianos
        # pontos, _   = [-1, 0, 1]                # Calculando nos nós
        # pontos      = [-1, 0, 1] if POS == 'nos' else self._elems[0].quadratura(grau)
        
        tau_xy      = np.array([])
        tau_xz      = np.array([])
        Y_gauss     = np.array([])  # Posições dos pontos gaussianos
        Z_gauss     = np.array([])
        Y_CT = self.propsArea['Y_CT']
        Z_CT = self.propsArea['Z_CT']

        # CALCULO DAS TENSÕES CISALHANTES
        for e,elem in enumerate(self._elems):

            # Graus de liberdade do elemento
            elementDof = np.array(self._elementsNodes)[e, :]
            # elementDof = self.elementNodes[e, :]
            # print(f'elementDof = {elementDof}')

            # Encontrando as coordenadas dos nós e.r. ao centro de torção
            ys = self._nodesCoords[elementDof,1] - Y_CT
            zs = self._nodesCoords[elementDof,2] - Z_CT
            # print(zs)

            # Encontrando os valores nodais de phi pro elemento
            phi = self.displacements[elementDof]
            
            # Propriedades do material
            G   = self.G[e]

            # Quadratura para o elemento; Varre os pontos gaussianos
            if self.tipoElem == 'Tri':
                for i, coords in enumerate(pontos):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,2,3)

                    # Funções de interpolação e Jacobiano
                    N = elem.getShapeFun(coords)
                    # print(f'{coords = }')
                    J, By, Bz = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                    
                    # tau_xy Normalizado _______________________________
                    z_gp    = N@zs                   # Posição do ponto gaussiano e.r. ao CT
                    txy     = By@phi - z_gp
                    tau_xy  = np.append(tau_xy, txy)
                    Z_gauss = np.append(Z_gauss, z_gp)

                    # tau_xz Normalizado _______________________________
                    y_gp    = N@ys                   # Posição do ponto gaussiano e.r. ao CT
                    txz     = Bz@phi + y_gp
                    tau_xz  = np.append(tau_xz, txz)
                    Y_gauss = np.append(Y_gauss, y_gp)

            elif self.tipoElem == 'Quad':
                for j, coord_j in enumerate(pontos):
                    for i, coord_i in enumerate(pontos):
                        #     i: indice
                        # coord: coordenadas naturais do ponto (zeta1,2,3)
                        coords = [coord_i, coord_j]
                        # print(coords)
                        # print(f'coord_i: {coord_i}')
                        # print(f'coord_j: {coord_j}')
                        # print('')
                        
                        # Funções de interpolação e Jacobiano
                        N = elem.getShapeFun(coords)
                        J, By, Bz = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                            
                        # tau_xy  _______________________________
                        Z       = N@zs
                        z_gp    = Z                   # Posição do ponto gaussiano e.r. ao CT
                        txy     = (By@phi - z_gp)*G
                        tau_xy  = np.append(tau_xy, txy)
                        Z_gauss = np.append(Z_gauss, Z)

                        # tau_xz  _______________________________
                        Y       = N@ys
                        y_gp    = Y                   # Posição do ponto gaussiano e.r. ao CT
                        txz     = (Bz@phi + y_gp)*G
                        tau_xz  = np.append(tau_xz, txz)
                        Y_gauss = np.append(Y_gauss, Y)

                #     print('Fim i')
                # print('Fim j')
                
        ##  -------------------------  PLOTS  -------------------------  ##
        fontkwargs = dict(fontname='Liberation Serif', usetex=True, fontsize=18)
        # As tensões dependem de alpha
        if taxa != None:
            if mises:
                tau_xy = np.sqrt(3)*taxa*tau_xy
                tau_xz = np.sqrt(3)*taxa*tau_xz
                labelBar = r'$\sqrt{3}| \tau | $' + f' [{unid}]'
            else:
                tau_xy = taxa*tau_xy
                tau_xz = taxa*tau_xz
                labelBar = r'$| \tau | $' + f' [{unid}]'                
        else:
            if mises:
                labelBar = r'$\sqrt{3}| \tau | / \alpha $' + fr' $[\mathrm{{{unid}}}\cdot\mathrm{{mm/rad}}]$'
            else:    
                labelBar = r'$| \tau | / \alpha $' + fr' $[\mathrm{{{unid}}}\cdot\mathrm{{mm/rad}}]$'
            
        # Seleciona o tipo de plot: vetorial ou escalar
        match modo:
            case 'vetorial':
                ## MOSTRA AS TENSÕES COMO CAMPO VETORIAL
                # Formato para o grafico vetorial
                Y       = Y_gauss
                Z       = Z_gauss
                match grau:
                    case 1:
                        npts = 1
                    case 2:
                        npts = 4
                    case 3|4:
                        npts = 9
                    case 5:
                        npts = 16
                
                fig, ax = plt.subplots(figsize=figsize, num='Tensoes Vetorial')
                
                mag     = np.sqrt(tau_xy**2 + tau_xz**2)
                    
                print(f'Max. magnitude: {np.max(mag)}')
                print(f'Med. magnitude: {np.mean(mag)}')
                print(f'Min. magnitude: {np.min(mag)}')
                if probeElems is not None:
                    print(f'Elem ({probeElems[0]: 4d}): {mag[npts*(probeElems[0])-1]:.6g}')
                    print(f'Elem ({probeElems[1]: 4d}): {mag[npts*(probeElems[1])-1]:.6g}')
            
                (U, V) = (tau_xy/mag, tau_xz/mag) if vetoresUnit else (tau_xy, tau_xz)
                    
                
                plot    = ax.quiver(Y[::passoVetores]+Y_CT, Z[::passoVetores]+Z_CT,
                                    U[::passoVetores], V[::passoVetores], 
                                    mag[::passoVetores], cmap='jet', pivot='mid', units='xy', **kwargs)
                cbar    = plt.colorbar(plot, **cbarKwargs)
                cbar.set_label(labelBar, **fontkwargs)

                ax.set_xlabel(r'$Y$', **fontkwargs)
                ax.set_ylabel(r'$Z$', **fontkwargs)
                # ax.set_title('Shear Stress Distribuition ($yz$-plane)'  , **fontkwargs)
                ax.set_aspect('equal', adjustable='box')

            case 'escalar':
                ## COMPONENTES COMO CAMPOS ESCALARES
                # import matplotlib.tri as tri 
                from matplotlib import gridspec, tri

                # Pondo os vetores em formato n_elem x n_gp
                Y_gauss     = np.reshape(Y_gauss,[self.nElems,-1]) + Y_CT
                Z_gauss     = np.reshape(Z_gauss,[self.nElems,-1]) + Z_CT
                coords      = [Y_gauss, Z_gauss]
                
                # Y   = np.array(self._nodesCoords)[:,1].reshape(-1)
                # Z   = np.array(self._nodesCoords)[:,2].reshape(-1)
                # YY, ZZ = Y_gauss, Z_gauss

                fig, ax = plt.subplots(figsize=(2*figsize[0],figsize[1]), ncols=2, sharey='row', num='Tensoes Escalar')                
                triang  = tri.Triangulation(Y_gauss.reshape(-1),Z_gauss.reshape(-1))

                ## Limites da geometria
                Y_min, Y_max = np.min(self._nodesCoords[:,0]),np.max(self._nodesCoords[:,0])
                Z_min, Z_max = np.min(self._nodesCoords[:,1]),np.max(self._nodesCoords[:,1])

                # Máximos e mínimos das tensões
                maxPlot     = np.max([tau_xy, tau_xz]) 
                minPlot     = np.min([tau_xy, tau_xz]) 
                color_kw    = dict(vmax=maxPlot, vmin=minPlot)
                # print(maxPlot, minPlot)

                ##  ------------------------ Tau_xy ------------------------  ##
                plot    = ax[0].tricontourf(triang, tau_xy, cmap='jet', **color_kw, **kwargs)
                # t = tau_xy.reshape(Y_gauss.shape)
                # plot    = ax[0].scatter(Y_gauss, Z_gauss, c=t, cmap='jet', **color_kw, **kwargs)
                # cbar    = plt.colorbar(plot, ax=axes[0], fraction=0.046, pad=0.04)
                # cbar.set_label(r'$\ \tau_{xy}  / G \alpha $'      , **fontkwargs)

                ax[0].set_title(r'$\tau_{xy}$', usetex=True, fontsize=24)
                ax[0].set_xlabel(r'$Y$'          , **fontkwargs)
                ax[0].set_ylabel(r'$Z$'          , **fontkwargs)
                # ax[0].set_xlim([Y_min, Y_max])
                # ax[0].set_ylim([Z_min, Z_max])
                ax[0].axis('square')


                ##  ------------------------ Tau_xz ------------------------  ##
                plot    = ax[1].tricontourf(triang, tau_xz, cmap='jet', **color_kw, **kwargs)
                # t = tau_xz.reshape(Y_gauss.shape)
                # plot    = ax[1].scatter(Y_gauss, Z_gauss, c=t, cmap='jet', **color_kw, **kwargs)
                cbar    = plt.colorbar(plot, ax=ax, fraction=0.035, pad=0.04, use_gridspec=True, **cbarKwargs)
                # cbar    = plt.colorbar(plot, ax=ax, fraction=0.035, pad=0.04, use_gridspec=True, shrink=.7)
                # cbar    = plt.colorbar(plot, ax=axes, fraction=0.035, pad=0.04, use_gridspec=True)
                # cbar    = plt.colorbar(plot, ax=axes, fraction=0.035, pad=0.04)
                cbar.set_label(labelBar , **fontkwargs)

                ax[1].set_title(r'$\tau_{xz}$', usetex=True, fontsize=24)
                ax[1].set_xlabel(r'$Y$'               , **fontkwargs)
                # axes[1].set_ylabel(r'$Z$'               , **fontkwargs)
                # ax[1].set_xlim([Y_min, Y_max])
                # ax[1].set_ylim([Z_min, Z_max])
                ax[1].axis('square')
                # plt.subplots_adjust(hspace=1)
                # fig.tight_layout()


            case _:
                raise Exception('Modo inválido')    

        # Caso queira mostrar a malha
        if mostrarMalha:
            from matplotlib.patches import Polygon
            from matplotlib.axes._axes import Axes
            for elem in self._elems:
                # Coordenadas dos nos
                nodesCoords = elem.getNCs()
                
                # Determinando a ordem dos vertices
                if self.nNodes == 9:
                    # Caso Q9
                    inds = [0,4,1,5,2,6,3,7]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                elif self.nNodes == 6:
                    # Caso T6
                    inds = [0,3,1,4,2,5]
                    Ys = np.array(nodesCoords)[inds,1]
                    Zs = np.array(nodesCoords)[inds,2]
                else:
                    # Caso elementos lineares
                    Ys = np.array(nodesCoords)[:,1]
                    Zs = np.array(nodesCoords)[:,2]

                # Vertices
                verts = np.vstack([Ys,Zs]).T

                # Polígonos
                poly    = Polygon(verts,fill=False,alpha=1,ec='k',lw=.3)
                if type(ax) == Axes:
                    # Caso o plote tenha um só axes (caso vetorial)
                    ax.add_patch(poly)
                elif type(ax[0]) == Axes:
                    # Caso o plote tenha mais de um axes (caso escalar)
                    poly2    = Polygon(verts,fill=False,alpha=1,ec='k',lw=.3)
                    
                    ax[0].add_patch(poly)
                    ax[1].add_patch(poly2)
                    # {axes.add_patch(poly) for axes in ax}
                else:
                    raise Exception('Erro no plot')
                
                
        return fig, ax
        
    

    ##  -----------------  Métodos Privados -----------------  ## 
    ## CONSTANTE DE TORÇÃO
    def __calcJphi(self, grau=2):
        '''Calcula a constante de torção'''
        pontos, pesos = self._elems[0].quadratura(grau)
        E0 = self.E[0] # np.min(self.E)
        G0 = self.G[0] # np.min(G)
        
        Jphi          = 0
        Jphi_pond_G   = 0
        Jphi_pond_E   = 0
        for e,elem in enumerate(self._elems):

            # Graus de liberdade do elemento
            elementDof = np.array(self._elementsNodes)[e, :]

            # Encontrando as coordenadas dos nós no sistema do CT
            ys = np.array(self._nodesCoords[elementDof,1]) - self.propsArea['Y_CT']
            zs = np.array(self._nodesCoords[elementDof,2]) - self.propsArea['Z_CT']

            # Encontrando o vetor de empenamento do elemento
            phi = self.displacements[elementDof].reshape(-1,1)
            
            
            # Propriedades do material
            G   = self.G[e]
            E   = self.E[e]

            # Quadratura para o elemento; Varre os pontos gaussianos
            if self.tipoElem == 'Tri':
                for i, coords in enumerate(pontos):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,zeta2,zeta3)

                    # Funções de interpolação e Jacobiano
                    N = elem.getShapeFun(coords)
                    J, By, Bz = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss


                    # Função Intengrada _______________
                    Fy = N@ys.reshape(-1,1)
                    Fz = N@zs.reshape(-1,1)

                    FF = Fy@Fy + Fz@Fz + Fy@Bz@phi - Fz@By@phi

                    Jphi += pesos[i]*J*FF
                    Jphi_pond_G += pesos[i]*J*FF * G/G0
                    Jphi_pond_E += pesos[i]*J*FF * E/E0                  
                    # print(f'{e = :4d}, {i = :4d}, {float(Jphi) = :.4f}, {pesos[i]*J*FF}')
                # print(' ')
                      

            elif self.tipoElem == 'Quad':
                for j, coord_j in enumerate(pontos):
                    for i, coord_i in enumerate(pontos):
                        #     i: indice
                        # coord: coordenadas naturais do ponto (eta_i,zeta_j)
                        coords = [coord_i, coord_j]
                        # print(f'{i = : 2d}, {j = : 2d},\n    {coords = }')
                        
                        # Funções de interpolação e Jacobiano
                        N = elem.getShapeFun(coords)
                        J, By, Bz = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss


                        # Função Intengrada _______________
                        Fy = N@ys.reshape(-1,1)
                        Fz = N@zs.reshape(-1,1)

                        FF = Fy@Fy + Fz@Fz + Fy@Bz@phi - Fz@By@phi

                        Jphi += pesos[i]*pesos[j]*J*FF
                        Jphi_pond_G += pesos[i]*pesos[j]*J*FF * G/G0
                        Jphi_pond_E += pesos[i]*pesos[j]*J*FF * E/E0
                        # print(f'{e = :4d}, {i = :4d}, {j = :4d}, {float(Jphi) = :.4f}, {pesos[i]*pesos[j]*J*FF}')
                # print(' ')
                        
            else:
                raise Exception('Tipo de elemento não reconhecido')
                
        self.propsArea.update(Jphi        = float(Jphi),
                              Jphi_pond_G = float(Jphi_pond_G),
                              Jphi_pond_E = float(Jphi_pond_E))
        # print(f'{J = }, \n{By = }, \n{Bz = }')
        # print(f'{pontos = }, \n{pesos = }')
        return
    ##  PROPRIEDADES DE ÁREA
    def __calcpropsArea(self): 
        '''Calcula as propriedades de área de um seção composta'''

        pontos, pesos = self._elems[0].quadratura(2)
        A, Qy, Qz, Iyy, Izz, Iyz = 0,0,0,0,0,0
        A_pond, Qy_pond , Qz_pond , Izz_pond, Iyy_pond, Iyz_pond = 0,0,0,0,0,0

        
        E0 = self.E[0] # np.min(self.E)
        G0 = self.G[0] # np.min(G)
        

        # for e in range(self.nElems):
        for e,elem in enumerate(self._elems):

            # Graus de liberdade do elemento
            elementDof = np.array(self._elementsNodes)[e, :]
            # print(f'elementDof = {elementDof}')

            # Encontrando as coordenadas dos nós
            ys = self._nodesCoords[elementDof,1]; #print(ys)
            zs = self._nodesCoords[elementDof,2]



            # Quadratura para o elemento; Varre os pontos gaussianos
            if self.tipoElem == 'Tri':
                for i, coords in enumerate(pontos):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,2,3)

                    # Funções de interpolação e Jacobiano
                    N = elem.getShapeFun(coords); #print(N)
                    J, _, _ = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                    
                    # ÁREA _______________________________
                    F = 1
                    A   += pesos[i]*J*F
                    if np.isclose(A, 0):
                        print('Elemento com área quase nula')
                        print(f'{e = }, {elem = }')

                    # MOMENTOS ESTÁTICOS DE ÁREA _________
                    Fy = float(N@ys.reshape(-1,1))
                    Fz = float(N@zs.reshape(-1,1))
                    
                    Qz  += pesos[i]*J*Fy
                    Qy  += pesos[i]*J*Fz

                    # MOMENTOS DE INÉRCIA DE ÁREA ________
                    Izz += pesos[i]*J*(Fy*Fy)

                    Iyy += pesos[i]*J*(Fz*Fz)

                    # PRODUTO DE INÉRCIA _________________
                    Iyz += pesos[i]*J*(Fy*Fz)
                    
                    # PROPRIEDADES PONDERADAS ____________
                    # Módulos elásticos
                    E_pond = self.E[e] / E0
                    G_pond = self.G[e] / G0

                    ## PONDERAÇÕES (POR E)
                    A_pond      += pesos[i]*J*F        * E_pond
                    Qz_pond     += pesos[i]*J*Fy       * E_pond
                    Qy_pond     += pesos[i]*J*Fz       * E_pond
                    Izz_pond    += pesos[i]*J*(Fy*Fy)  * E_pond
                    Iyy_pond    += pesos[i]*J*(Fz*Fz)  * E_pond
                    Iyz_pond    += pesos[i]*J*(Fy*Fz)  * E_pond                    

            elif self.tipoElem == 'Quad':
                for j, coord_j in enumerate(pontos):
                    for i, coord_i in enumerate(pontos):
                        #     i: indice
                        # coord: coordenadas naturais do ponto (zeta1,2,3)
                        coords = [coord_i, coord_j]

                        # Funções de interpolação e Jacobiano
                        N = elem.getShapeFun(coords); #print(N)
                        J, _, _ = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                        
                        # ÁREA _______________________________
                        F    = 1
                        A   += pesos[i]*pesos[j]*J*F

                        # MOMENTOS ESTÁTICOS DE ÁREA _________
                        Fy   = N@ys.reshape(-1,1)
                        Fz   = N@zs.reshape(-1,1)
                        # O produto resulta em um array 1x1, passando para float    
                        Fy   = float(Fy)
                        Fz   = float(Fz)
                        Qz  += pesos[i]*pesos[j]*J*Fy
                        Qy  += pesos[i]*pesos[j]*J*Fz

                        # MOMENTOS DE INÉRCIA DE ÁREA ________
                        Izz += pesos[i]*pesos[j]*J*(Fy*Fy)

                        Iyy += pesos[i]*pesos[j]*J*(Fz*Fz)

                        # PRODUTO DE INÉRCIA _________________
                        Iyz += pesos[i]*pesos[j]*J*(Fy*Fz)

                        # PROPRIEDADES PONDERADAS ____________
                        # Módulos elásticos
                        E_pond = self.E[e] / E0
                        G_pond = self.G[e] / G0

                        ## PONDERAÇÕES (POR E)
                        A_pond      += pesos[i]*pesos[j]*J*F        * E_pond
                        Qz_pond     += pesos[i]*pesos[j]*J*Fy       * E_pond
                        Qy_pond     += pesos[i]*pesos[j]*J*Fz       * E_pond
                        Izz_pond    += pesos[i]*pesos[j]*J*(Fy*Fy)  * E_pond
                        Iyy_pond    += pesos[i]*pesos[j]*J*(Fz*Fz)  * E_pond
                        Iyz_pond    += pesos[i]*pesos[j]*J*(Fy*Fz)  * E_pond

        ## Determinação do centroide ponderado (sempre por E)
        Y_CG = float(Qz_pond / A_pond)
        Z_CG = float(Qy_pond / A_pond)
        ## Propriedades de flexão no Centroide (teorema dos eixos paralelos)
        Iyy_cent = Iyy_pond - A_pond * Z_CG**2
        Izz_cent = Izz_pond - A_pond * Y_CG**2
        Iyz_cent = Iyz_pond - A_pond * Y_CG*Z_CG

        # Sem poderação
        # self.propsArea.update(A=A, Qy=float(Qy), Qz=float(Qz), Iyy=float(Iyy), Izz=float(Izz), Iyz=float(Iyz), Y_CG=Y_CG, Z_CG=Z_CG)

        # Ponderação por E
        self.propsArea.update(
                                A_pond      = float(A_pond), 
                                Qy_pond     = float(Qy_pond), 
                                Qz_pond     = float(Qz_pond), 
                                Iyy_pond    = float(Iyy_pond), 
                                Izz_pond    = float(Izz_pond), 
                                Iyz_pond    = float(Iyz_pond), 
                                Y_CG_pond   = Y_CG, 
                                Z_CG_pond   = Z_CG
                              )

        # Propriedades de flexão (com respeito ao centroide)
        self.propsArea.update(
                                Iyy_pond_cent = float(Iyy_cent), 
                                Izz_pond_cent = float(Izz_cent), 
                                Iyz_pond_cent = float(Iyz_cent)
                              )
        return
    ##  PROPRIEDADES DE INERCIA
    def __calcpropsInercia(self): 
        '''Calcula as propriedades de inercia de um seção composta'''

        pontos, pesos = self._elems[0].quadratura(2)
        rho_l       , Qy_rho    , Qz_rho        = 0,0,0     #       ρ dA,     ρ*y dA,     ρ*z dA
        Izz_rho     , Iyy_rho   , Iyz_rho       = 0,0,0     #   ρ*z^2 dA,   ρ*y^2 dA,   ρ*y*z dA
        phi_rho     , phi_phi_rho               = 0,0       #     ρ*φ dA,   ρ*φ^2 dA
        phi_y_rho   , phi_z_rho                 = 0,0       #   ρ*φ*y dA,   ρ*φ*z dA
        phi_yy_rho  , phi_yz_rho, phi_zz_rho    = 0,0,0     # ρ*φ*y^2 dA, ρ*φ*z^2 dA, ρ*φ*y*z dA
        
        Qy_cent_pond, Qz_cent_pond,Iyz_cent_pond = 0,0,0     # E/E0*y dA,  E/E0*y dA, E/E0*y*z dA
        
        # Obtendo a posição do centroide
        Y_CG = self.propsArea['Y_CG_pond']
        Z_CG = self.propsArea['Z_CG_pond']
        
        # Obtendo a posição do centro de torção
        Y_CT = self.propsArea['Y_CT']
        Z_CT = self.propsArea['Z_CT']

        

        # for e in range(self.nElems):
        for e,elem in enumerate(self._elems):

            # Graus de liberdade do elemento
            elementDof = np.array(self._elementsNodes)[e, :]
            # print(f'elementDof = {elementDof}')

            # Encontrando as coordenadas dos nós no referencial local (centroide da seção)
            ys = self._nodesCoords[elementDof,1] - Y_CG
            zs = self._nodesCoords[elementDof,2] - Z_CG
            
            # Coordenadas e.r. ao CT
            # yts = self._nodesCoords[elementDof,1] - Y_CT
            # zts = self._nodesCoords[elementDof,2] - Z_CT
            
            # Determinando os valores nodais da funcao empenamento
            phi = self.displacements[elementDof]

            # Propriedades do elemento
            rho     = self.rho[e]
            E_pond  = self.E[e] / self.E[0]

            # Quadratura para o elemento; Varre os pontos gaussianos
            if self.tipoElem == 'Tri':
                for i, coords in enumerate(pontos):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,2,3)

                    # Funções de interpolação e Jacobiano
                    N = elem.getShapeFun(coords); #print(N)
                    J, _, _ = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                    WJ = pesos[i]*J

                    # INTERPOLAÇÃO DE Y, Z e PHI
                    Fy   = N@ys.reshape(-1,1)
                    Fz   = N@zs.reshape(-1,1)
                    Fphi = N@phi.reshape(-1,1)
                    # Fty  = N@yts.reshape(-1,1)
                    # Ftz  = N@zts.reshape(-1,1)

                    ## PROPRIEDADES 
                    rho_l           += WJ*1         * rho           #       ρ dA
                    Qz_rho          += WJ*Fy        * rho           #     ρ*y dA
                    Qy_rho          += WJ*Fz        * rho           #     ρ*z dA
                    Izz_rho         += WJ*Fy@Fy     * rho           #   ρ*z^2 dA
                    Iyy_rho         += WJ*Fz@Fz     * rho           #   ρ*y^2 dA
                    Iyz_rho         += WJ*Fy@Fz     * rho           #   ρ*y*z dA
                    
                    ## PROPRIEDADES DEPENDENTES DE PHI
                    phi_rho         += WJ*Fphi         * rho        #     ρ*φ dA
                    phi_phi_rho     += WJ*Fphi@Fphi    * rho        #   ρ*φ^2 dA
                    phi_y_rho       += WJ*Fy@Fphi      * rho        #   ρ*φ*y dA
                    phi_z_rho       += WJ*Fz@Fphi      * rho        #   ρ*φ*z dA
                    phi_yy_rho      += WJ*Fy@Fy@Fphi   * rho        # ρ*φ*y^2 dA
                    phi_zz_rho      += WJ*Fz@Fz@Fphi   * rho        # ρ*φ*z^2 dA
                    phi_yz_rho      += WJ*Fy@Fz@Fphi   * rho        # ρ*φ*y*z dA
                    
                    ## ATUALIZAÇÃO DOS MOMENTOS ESTÁTICOS (NO CENTROIDE PONDERADO)
                    Qz_cent_pond    += WJ*Fy         * E_pond   #   E/E0*y dA
                    Qy_cent_pond    += WJ*Fz         * E_pond   #   E/E0*z dA
                    Iyz_cent_pond   += WJ*Fy@Fz      * E_pond   #      y*z dA

            elif self.tipoElem == 'Quad':
                for j, coord_j in enumerate(pontos):
                    for i, coord_i in enumerate(pontos):
                        #     i: indice
                        # coord: coordenadas naturais do ponto (zeta1,2,3)
                        coords = [coord_i, coord_j]

                        # Funções de interpolação e Jacobiano
                        N = elem.getShapeFun(coords); #print(N)
                        J, _, _ = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                        WJ = pesos[i]*pesos[j]*J  # Produto dos pesos e jacobiano
                        
                        # INTERPOLAÇÃO DE Y, Z e PHI
                        Fy   = N@ys.reshape(-1,1)
                        Fz   = N@zs.reshape(-1,1)
                        Fphi = N@phi.reshape(-1,1)
                        # Fty  = N@yts.reshape(-1,1)
                        # Ftz  = N@zts.reshape(-1,1)

                        ## PROPRIEDADES 
                        rho_l           += WJ*1         * rho           #       ρ dA
                        Qy_rho          += WJ*Fy        * rho           #     ρ*y dA
                        Qz_rho          += WJ*Fz        * rho           #     ρ*z dA
                        Izz_rho         += WJ*Fy@Fy     * rho           #   ρ*z^2 dA
                        Iyy_rho         += WJ*Fz@Fz     * rho           #   ρ*y^2 dA
                        Iyz_rho         += WJ*Fy@Fz     * rho           #   ρ*y*z dA
                        
                        ## PROPRIEDADES DEPENDENTES DE PHI
                        phi_rho         += WJ*Fphi         * rho        #     ρ*φ dA
                        phi_phi_rho     += WJ*Fphi@Fphi    * rho        #   ρ*φ^2 dA
                        phi_y_rho       += WJ*Fy@Fphi      * rho        #   ρ*φ*y dA
                        phi_z_rho       += WJ*Fz@Fphi      * rho        #   ρ*φ*z dA
                        phi_yy_rho      += WJ*Fy@Fy@Fphi   * rho        # ρ*φ*y^2 dA
                        phi_zz_rho      += WJ*Fz@Fz@Fphi   * rho        # ρ*φ*z^2 dA
                        phi_yz_rho      += WJ*Fy@Fz@Fphi   * rho        # ρ*φ*y*z dA
                        
                        ## ATUALIZAÇÃO DOS MOMENTOS ESTÁTICOS (NO CENTROIDE PONDERADO)
                        Qz_cent_pond    += WJ*Fy         * E_pond   #     E/E0*y dA
                        Qy_cent_pond    += WJ*Fz         * E_pond   #     E/E0*z dA
                        Iyz_cent_pond   += WJ*Fy@Fz      * E_pond   #   E/E0*y*z dA
                        
                        
        # ATUALIZANDO AS PROPRIEDADES DE INERCIA
        self.propsInercia.update(   
                                    rho_l         = float(rho_l      ),
                                    Qy_rho        = float(Qy_rho     ),
                                    Qz_rho        = float(Qz_rho     ),
                                    Iyy_rho       = float(Iyy_rho    ),
                                    Izz_rho       = float(Izz_rho    ),
                                    Iyz_rho       = float(Iyz_rho    ),
                                    I_phi_rho     = float(phi_rho    ),
                                    I_phi_phi_rho = float(phi_phi_rho),
                                    I_phi_y_rho   = float(phi_y_rho  ),
                                    I_phi_z_rho   = float(phi_z_rho  ),
                                    I_phi_yy_rho  = float(phi_yy_rho ),
                                    I_phi_zz_rho  = float(phi_zz_rho ),
                                    I_phi_yz_rho  = float(phi_yz_rho )
                                )
        # ATUALIZANDO AS PROPRIEDADES DE AREA
        self.propsArea.update(
                                    Qy_cent_pond  = float(Qy_cent_pond),
                                    Qz_cent_pond  = float(Qz_cent_pond),
                                    Iyz_cent_pond = float(Iyz_cent_pond),
        )
        return
    ##  PROCESSAMENTO 
    def __rigidezGlobal (self, grau=4):
        '''Monta a matriz de rigidez Global
        ------------------------------------------------
        USO:
            Cheme este método na instâcia do problema:
        EX:
            >>> elementNodes = numpy.array([[0,1,2], [2,3,0]])
            >>> nodeCoords   = numpy.array([[0,0],[1,0],[1,1],[0,1]])
            >>> problema     = Torcao(elementNodes, nodeCoords)
            >>> problema.rigidezGlobal()'''
        
        # Loops para os elementos
        for e, elem in enumerate(self._elems):

            #print(f'Elem.: {e}')
    
            # Graus de liberdade do elemento
            elementDof = np.array(self._elementsNodes)[e, :]

            # Encontrando as coordenadas dos nós
            ys = self._nodesCoords[elementDof,1]
            zs = self._nodesCoords[elementDof,2]

            # Matriz de rigidez local do elemento
            kLocal = self.G[e]*elem.getKlocal(grau)
            
            # Montagem da matriz de rigidez global
            cols, lins     = np.meshgrid(range(self.nNodes),range(self.nNodes))
            dofCol, dofLin = np.meshgrid(elementDof,elementDof)
            self.stiffness[dofLin,dofCol] += kLocal[lins,cols] 

            # Montagem do Vetor de Forças Global
            fLocal = self.G[e]*elem.getFlocal(grau)

            self.forces[dofLin,:] += fLocal[lins,:]

        return
    #   CONDIÇÕES DE CONTORNO
    def __condContorno (self, prescribedDof):
        '''Calcula os GDLs livres, dados os GDLs restritos
        ------------------------------------------------
        USO:
            Crie um numpy array com os nós restritos e passe-o 
            para este método na instância do problema.
        EX:
            prescribedDof = np.array([[1, 3, 4]])
            problema.condContorno(prescribedDof)'''
        
        # Nós com deslocamento prescrito
        self.prescribedDof = prescribedDof # ajustando pra python
        
        # Nós com deslocamento livre
        self.activeDof = np.setdiff1d(self._elementsNodes, self.prescribedDof)
        return
    def __solucao(self, grau):
        '''Soluciona [K]{d} = {f} para os nós livres
        ------------------------------------------------
        USO:
            Chame este método na instancia do problema.
        EX:
            problema.solucao()'''
        
        from scipy.linalg import solve as sol
        
        livres = self.activeDof # nós sem restriçao
        gr1,gr2 = np.meshgrid(livres,livres)

        stiffTemp  = self.stiffness[gr1,gr2]
        forcesTemp = self.forces[livres]

        ################################################################
        # Solução do sistema
        dispTemp =  sol(stiffTemp, forcesTemp)
        
        # Preenchendo o vetor de deslocamentos com os desloc. livres
        j = 0
        for i in self.activeDof:
            self.phiStar[i, 0] = dispTemp[j]
            j += 1
               
        # Solucionando o vetor de forças
        self.forces = self.stiffness@self.phiStar

        ### Solução do centro de torção
        self.__solucaoCT(grau)

        ### Correção da solução
        Y_CT = self.propsArea['Y_CT']
        Z_CT = self.propsArea['Z_CT']
        c_CT = self.propsArea['c_CT']
        ys   = np.array(self._nodesCoords[:,1]).reshape(-1,1)
        zs   = np.array(self._nodesCoords[:,2]).reshape(-1,1)

        self.displacements = self.phiStar + Y_CT*zs - Z_CT*ys - c_CT #

        ### Calculo da Constante de Toção
        self.__calcJphi(grau)
        
        ### Calculo do empenamento restrito?
        # self.restrained = 
        
        return
    def __solucaoCT(self, grau):
        '''Encontra a posição do centro de torção'''
        pontos, pesos = self._elems[0].quadratura(grau)                     ############################
        Iphi,            Iphiy,      Iphiz = 0,0,0
        Iphi_pond , Iphiy_pond, Iphiz_pond = 0,0,0
        E0 = self.E[0]

        # Integrais de phi*
        for e,elem in enumerate(self._elems):

            # Graus de liberdade do elemento
            elementDof = np.array(self._elementsNodes)[e, :]
            # print(f'elementDof = {elementDof}')

            # Encontrando as coordenadas dos nós
            ys = self._nodesCoords[elementDof,1]
            zs = self._nodesCoords[elementDof,2]

            # Encontrando os valores nodais de phi* pro elemento
            phiStar = self.phiStar[elementDof]

            # Integrais de phi*
            if self.tipoElem == 'Tri':
                for i, coords in enumerate(pontos):
                    #     i: indice
                    # coord: coordenadas naturais do ponto (zeta1,2,3)

                    # Funções de interpolação e Jacobiano
                    N = elem.getShapeFun(coords)
                    J, _, _ = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                    
                    # phi  _______________________________
                    Fphi  = N@phiStar.reshape(-1,1)
                    # Iphi += pesos[i]*J*Fphi

                    # Y*phi _______________________________
                    Fy     = N@ys.reshape(-1,1)
                    # Iphiy += pesos[i]*J*Fy@Fphi

                    # Z*phi _______________________________
                    Fz     = N@zs.reshape(-1,1)
                    # Iphiz += pesos[i]*J*Fz@Fphi
                    
                    ## PONDERAÇÃO POR E
                    E_pond = self.E[e] / E0                    

                    ## PONDERAÇÕES (POR E)
                    Iphi_pond    += pesos[i]*J*Fphi    * E_pond
                    Iphiy_pond   += pesos[i]*J*Fy@Fphi * E_pond
                    Iphiz_pond   += pesos[i]*J*Fz@Fphi * E_pond

            elif self.tipoElem == 'Quad':
                for j, coord_j in enumerate(pontos):
                    for i, coord_i in enumerate(pontos):
                        #     i: indice
                        # coord: coordenadas naturais do ponto (zeta1,2,3)
                        coords = [coord_i, coord_j]
                        N = elem.getShapeFun(coords)
                        J, _, _ = elem.getJacob(coords) # (ys,zs): nós do elemento, coords: coordenadas zeta dos pnt gauss
                        
                        # phi  _______________________________
                        Fphi  = N@phiStar.reshape(-1,1)
                        # Iphi += pesos[i]*pesos[j]*J*Fphi

                        # Y*phi _______________________________
                        Fy     = N@ys.reshape(-1,1)
                        # Iphiy += pesos[i]*pesos[j]*J*Fy@Fphi

                        # Z*phi _______________________________
                        Fz     = N@zs.reshape(-1,1)
                        # Iphiz += pesos[i]*pesos[j]*J*Fz@Fphi

                        ## PONDERAÇÃO POR E
                        E_pond = self.E[e] / E0

                        ## PONDERAÇÕES (POR E)
                        Iphi_pond    += pesos[i]*pesos[j]*J*Fphi    * E_pond
                        Iphiy_pond   += pesos[i]*pesos[j]*J*Fy@Fphi * E_pond
                        Iphiz_pond   += pesos[i]*pesos[j]*J*Fz@Fphi * E_pond
        
        # Montagem do sistema de equações
        A   = self.propsArea['A_pond']
        Qy  = self.propsArea['Qy_pond']
        Qz  = self.propsArea['Qz_pond']
        Izz = self.propsArea['Izz_pond']
        Iyy = self.propsArea['Iyy_pond']
        Iyz = self.propsArea['Iyz_pond']

        M = np.array([[ -Qy, Qz, A],     # Matriz dos coeficientes [editada 10/04/25]
                      [-Iyz,Izz,Qz],
                      [-Iyy,Iyz,Qy]])
        # M = np.array([[ -Qz, Qy, A],     # Matriz dos coeficientes [anterior]
        #               [-Iyz,Iyy,Qy],
        #               [-Izz,Iyz,Qz]])
        # b = np.array([Iphi,Iphiy,Iphiz]) # Vetor dos termos independentes
        b = np.array([Iphi_pond,Iphiy_pond,Iphiz_pond]) # Vetor dos termos independentes

        # Vetor Solução
        x = np.linalg.inv(M)@b.reshape(-1,1)
        Y_CT, Z_CT, c_CT = x

        self.propsArea.update(Y_CT=float(Y_CT), Z_CT=float(Z_CT), c_CT=float(c_CT))
    
        return











