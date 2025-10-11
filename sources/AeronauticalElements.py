
##  ================================================================  ##
##                  CLASSES DE VIGAS 'AERONAUTICAS'                   ##
##       DESATUALIZADO -- NECESSITA HERANÇAS DA CLASSE DE VIGAS       ##
##  ================================================================  ##

# Elementos 'aeronauticos' são modelos de viga especiais para facilitar a
# construção de modelos de simulações de voo. 
from StructuralElements import VigaAssim3D

class VigaFuselagem():
    def __init__(self, *, secao:Secao2D, Xa0:float, comp:float, r_ext:float, esp:float, nElem:int):
        '''
        Modelo de fuselagem: seção circular oca. 
        Linha de referência: 'barriga', Za == 0
        Argumentos (Keyword):
        ---
        secao   : objeto <Secao> da seção transversal
        Xa0_vet : vetor posição inicial do nariz
        comp    : comprimento
        r_ext   : raio externo
        esp     : espessura da parade
        nElem   : número de elemento ao longo do comprimento
        '''

        ## Atributos relacionados à seção transversal
        self.rho_ls    = np.ones(nElem) * secao.propArea['rho_l']
        self.Qy_rhos   = np.ones(nElem) * secao.propArea['Qy_rho']
        self.Qz_rhos   = np.ones(nElem) * secao.propArea['Qz_rho']

        ## Atributos
        self.comp      = comp
        self.r_ext     = r_ext
        self.esp       = esp

        self.nElem     = nElem
        self.nNodes    = nElem +1
        self.elemNodes = np.concatenate([
            np.arange(self.nElem).reshape(-1,1),
            np.arange(self.nElem).reshape(-1,1)+1],
            axis=1) # Conectividade Local

        ## Criação da malha de vigas
        self.Xas = np.linspace(Xa0, Xa0+comp, self.nNodes)
        self.Yas = np.zeros(self.nNodes)
        self.Zas = .5*r_ext * np.ones(self.nNodes)
        self.Ls  = np.sqrt( np.diff(self.Xas)**2+np.diff(self.Yas)**2+np.diff(self.Zas)**2 )

class VigaAsa():
    def __init__(self, *, secao:Secao2D, Xa0_vet_raiz:NDArray[Shape['3, 1'], Float], env:float, corda:float, direc:str, nElem:int):
        '''
        Modelo de Asa.
        Argumentos (Keyword):
        ---
        secao         : objeto <Secao> da seção transversal
        Xa0_vet_raiz  : vetor posição inicial da raíz da asa
        env           : envergadura da asa
        corda         : dimensão da corda
        direc         : direção da ponta da asa a partir da raíz ('+y' ou '-y')
        nElem         : número de elemento ao longo do comprimento
        '''

        ## Atributos relacionados à seção transversal
        self.rho_ls    = np.ones(nElem) * secao.propArea['rho_l']
        self.Qy_rhos   = np.ones(nElem) * secao.propArea['Qy_rho']
        self.Qz_rhos   = np.ones(nElem) * secao.propArea['Qz_rho']

        ## Atributos
        self.env       = env
        self.corda     = corda
        self.direc     = direc

        self.nElem     = nElem
        self.nNodes    = nElem +1
        self.elemNodes = np.concatenate([
            np.arange(self.nElem).reshape(-1,1),
            np.arange(self.nElem).reshape(-1,1)+1],
            axis=1) # Conectividade Local

        ## Criação da malha de vigas
        X0, Y0, Z0 = Xa0_vet_raiz
        Yf = Y0 + env if self.direc == '+y' else Y0 - env
        self.Xas   = X0 * np.ones(self.nNodes)
        self.Yas   = np.linspace(Y0,Yf,self.nNodes)
        self.Zas   = Z0 * np.ones(self.nNodes)
        self.Ls    = np.sqrt( np.diff(self.Xas)**2+np.diff(self.Yas)**2+np.diff(self.Zas)**2 )

class VigaEmpHor(VigaAsa):
    def __init__(self, *, secao:Secao2D, Xa0_vet_raiz:NDArray[Shape['3, 1'], Float], env:float, corda:float, direc:str, nElem:int):
        '''
        Modelo de empenagem vertical.
        Argumentos (Keyword):
        ---
        secao         : objeto <Secao> da seção transversal
        Xa0_vet_raiz  : vetor posição inicial da raíz da asa
        env           : envergadura da asa
        corda         : dimensão da corda
        direc         : direção da ponta da asa a partir da raíz ('+y' ou '-y')
        nElem       
        '''
        super().__init__(secao=secao, Xa0_vet_raiz=Xa0_vet_raiz, env=env, corda=corda, direc=direc, nElem=nElem)

class VigaEmpVer(VigaAsa):
    def __init__(self, *, secao:Secao2D, Xa0_vet_raiz:NDArray[Shape['3, 1'], Float], env:float, corda:float, nElem:int):
        '''
        Modelo de empenagem vertical.
        Argumentos (Keyword):
        ---
        secao         : objeto <Secao> da seção transversal
        Xa0_vet_raiz  : vetor posição inicial da raíz da asa
        env           : envergadura da asa
        corda         : dimensão da corda
        nElem         : número de elemento ao longo do comprimento
        '''
        super().__init__(secao=secao, Xa0_vet_raiz=Xa0_vet_raiz, env=env, corda=corda, direc='+z', nElem=nElem)

        ## Criação da malha de vigas
        X0, Y0, Z0 = Xa0_vet_raiz
        Zf = Z0 + env if self.direc == '+z' else Y0 - env
        self.Xas   = X0 * np.ones(self.nNodes)
        self.Yas   = Y0 * np.ones(self.nNodes) 
        self.Zas   = np.linspace(Z0,Zf,self.nNodes)
        self.Ls    = np.sqrt( np.diff(self.Xas)**2+np.diff(self.Yas)**2+np.diff(self.Zas)**2 )