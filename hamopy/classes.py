# Copyright (C) 2014  Simon Rouchier
# 
# This file is part of Hamopy.
#
# Hamopy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Hamopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Hamopy.  If not, see <http://www.gnu.org/licenses/>.
# 
# To contact the author : s.rouchier@gmail.com

from __future__ import division

import numpy  as np
import pandas as pd
from . import ham_library as ham

from pylab             import sqrt
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.sparse      import coo_matrix, hstack, vstack

class Material(object):
    """
    The :class:`Material` class includes methods for the definition of all heat
    and moisture properties of a material, and for calling these properties
    as functions of the capillary pressure and temperature
    
    The following arguments can be given at the definition of a new
    :class:`Material` object :
    
    :name: Name of the material (string)
    :rho: Dry density (float)
    :cp: Heat capacity (float)
    
    Additionnaly, the following methods must be called prior to the integration
    of the material into a :class:`Mesh` as to define its properties :
    
    :set_density: Set the material's dry density rho [kg/m3]
    :set_capacity: Set the heat capacity :math:`c_p` and its eventual
                   dependency to the temperature
    :set_conduc: Set the heat conductivity :math:`\lambda` and its eventual
                 dependency to moisture content and temperature
    :set_isotherm: Set the profile of the sorption isotherm
    :set_perm_vapor: Set the profile of water vapour permeability
                     :math:`\delta_p` [s]
    :set_perm_liquid: (optional) Set the profile of liquid permeability
                      :math:`k_l` [s]
    :set_perm_air: (optional) Set the value of air permeability
                   :math:`k_\mathit{air}` [m2]
    
    Default values are rho :math:`= 1000`, :math:`c_p = 1000`,
    :math:`k_\mathit{air} = 0` and :math:`k_l = 0`
    
    Once defined, the properties can be called as functions of P and T with
    the following methods
    
    :cp: Heat capacity (can be a function of T)
    :conduc: Heat conductivity (can be a function of P and T)
    :w: Moisture content
    :c_w: Derivative of the moisture content
    :delta_p: Vapour permeability
    :k_l: Liquid permeability
    """
    def __init__(self, name = 'new material', rho = 1000., cp = 1000.):
        
        # Initialisation : density and heat capacity may be given here
        self.name = name
        self.rho  = rho
        self.cp_0 = cp
        self.cp_t = 0
        self.k_air = 0

    def set_density(self, rho):
        """
        Register the density of the dry material [kg/m3]
        
        One argument (float) needed
        """
        
        self.rho = rho
    
    def set_capacity(self, cp_0, cp_t = 0):
        """
        Register the heat capacity of the material [J/(kg.K)]
        :math:`c_p = c_\mathit{p,0} + T(^\circ C) c_\mathit{p,t}`
        
        :cp_0: Value of :math:`c_p` at :math:`T=0^\circ C`
        :cp_t: (optional) Variation of :math:`c_p` with the temperature.
               Default is 0
        """
        
        self.cp_0 = cp_0
        self.cp_t = cp_t
    
    def cp(self, T = 273.15):
        """
        Heat capacity [J/(kg.K)]
        
        The temperature T is an optional argument
        """
        
        return self.cp_0 + (T-273.15)*self.cp_t
    
    def set_conduc(self, lambda_0, lambda_m = 0, lambda_t = 0):
        """
        Register the heat conductivity of the material [W/(m.K)]
        :math:`\lambda = \lambda_0 + w/1000 \lambda_m + T(^\circ C) \lambda_t`
        
        :lambda_0: Value of :math:`\lambda` at :math:`T=0^\circ C`
        :lambda_m: (optional) Variation of :math:`\lambda` with moisture content
        :lambda_t: (optional) Variation of :math:`\lambda` with the temperature
        """
        
        self.lambda_0 = lambda_0
        self.lambda_m = lambda_m
        self.lambda_t = lambda_t
        
    def conduc(self, p_c = -np.inf, T = 293.15):
        """
        Heat conductivity [W/(m.K)]
        
        The capillary pressure p_c and the temperature T are optional arguments
        """
        
        w = self.w(p_c)
        return self.lambda_0 + w/1000.*self.lambda_m + self.lambda_t * (T-273.15)
    
    def set_isotherm(self, method, **kwargs):   
        """
        Register the sorption isotherm [kg/m3]
        
        :method: String indicating how the isotherm is defined. Available
                 methods are 'vangenuchten' and 'polynomial'
        :kwargs: Dictionary of parameters for the definition of the profile.
                 See the wiki of the hamopy project for instructions regarding
                 the definition of this object.
        """
        
        self.w_method = method        
        
        if method == 'vangenuchten':
            # Mono- or multimodal van Genuchten law
            self.w_sat  = kwargs['w_sat']
            if isinstance(kwargs['m'], float):
                self.w_l      = np.array([kwargs['l']])
                self.w_alpha  = np.array([kwargs['alpha']])
                self.w_m      = np.array([kwargs['m']])
            else:
                self.w_l      = np.array(kwargs['l'])
                self.w_alpha  = np.array(kwargs['alpha'])
                self.w_m      = np.array(kwargs['m'])
                
        elif method == 'polynomial':
            # Polynomial expression of w [kg/m3] as a function of HR
            self.w_poly = np.polyfit(kwargs['HR'], kwargs['W'], 3)
            # Interpolation points
            self.w_p1, self.w_p2, self.w_p3, self.w_p4 = kwargs['W']
            
        elif method == 'slope':
            
            self.xi_poly = np.polyfit(kwargs['HR'], kwargs['XI'], 2)
            # Interpolation points
            self.xi_p1, self.xi_p2, self.xi_p3 = kwargs['XI']
            
    def w(self, p_c, T = 293.15):
        """
        Moisture content [kg/m3]
        
        The capillary pressure p_c is required, and T is an optional argument
        """
        
        if self.w_method == 'vangenuchten':
        
            w = np.zeros(np.shape(p_c))
            n = 1./(1-self.w_m)
            for i in range(np.size(self.w_l)):
                w += self.w_sat * self.w_l[i] * \
                (1.+(self.w_alpha[i]*abs(p_c))**n[i])**(-self.w_m[i])
            
        elif self.w_method == 'polynomial':

            w = np.polyval(self.w_poly, ham.HR(p_c, T))
            
        elif self.w_method == 'slope':
            
            w = np.polyval(np.polyint(self.xi_poly), ham.HR(p_c, T))
            
        return w
        
    def c_w(self, p_c, T = 293.15):
        """
        Derivative of the moisture content over the capillary pressure
        
        :math:`c_w = \partial w / \partial p_c`
        
        The capillary pressure p_c is required, and T is an optional argument
        """
        
        if self.w_method == 'vangenuchten':
            
            c_w = np.zeros(np.shape(p_c))
            n = 1./(1-self.w_m)
            for i in range(np.size(self.w_l)):
                c_w += self.w_sat*self.w_l[i]*self.w_m[i]*n[i] / abs(p_c) * \
                (abs(p_c)*self.w_alpha[i])**n[i] * \
                (1.+(self.w_alpha[i]*abs(p_c))**n[i])**(-self.w_m[i]-1)
                
        elif self.w_method == 'polynomial':
            
            HR  = ham.HR(p_c, T)
            xi  = np.polyval(np.polyder(self.w_poly), HR)
            c_w = xi * HR / (ham.rho_liq*ham.Rv*T)
            
        elif self.w_method == 'slope':
            
            HR  = ham.HR(p_c, T)
            xi  = np.polyval(self.xi_poly, HR)
            c_w = xi * HR / (ham.rho_liq*ham.Rv*T)
            
        return c_w
    
    def set_perm_vapor(self, method, **kwargs):
        """
        Register the water vapor permeability [s]
        
        :method: String indicating how the permeability is defined. Available
                 methods are 'schirmer' and 'interp'
        :kwargs: Dictionary of parameters for the definition of the profile.
                 See the wiki of the hamopy project for instructions regarding
                 the definition of this object.
        """
        
        self.dp_method = method
        
        if method == 'schirmer':
            
            # dp_mu = permeability of the dry material, p = evolution avec HR
            self.dp_mu = kwargs['mu']
            self.dp_p  = kwargs['p']
            
        elif method == 'interp':
            
            # delta_p is defined by interpolation between measured points
            #self.dp_p1, self.dp_p2 = kwargs['dp']
            #self.dp_f = UnivariateSpline(kwargs['HR'], np.log10(kwargs['dp']), k=1 )
            self.dp_f = interp1d(np.array(kwargs['HR']), np.log10(kwargs['dp']), fill_value = 'extrapolate')
            
        elif method == 'interp_mu':
            
            self.mu_f = interp1d(np.array(kwargs['HR']), np.array(kwargs['MU']), fill_value = 'extrapolate')
            
    def delta_p(self, p_c, T = 293.15):
        """
        Water vapour permeability [s]
        
        The capillary pressure p_c is required, and T is an optional argument
        """
        
        if self.dp_method == 'schirmer':
        
            w  = self.w(p_c)
            dp = 26.1e-6 / (self.dp_mu*ham.Rv*T) * (1-w/self.w_sat) / \
            ( (1-self.dp_p)*(1-w/self.w_sat)**2 + self.dp_p )
            
        elif self.dp_method == 'interp':
            
            dp = 10 ** self.dp_f( ham.HR(p_c, T).ravel() )
            dp = dp.reshape(np.shape(p_c))
            
        elif self.dp_method == 'interp_mu':
            
            hr = ham.HR(p_c, T)
            mu = self.mu_f(hr)
            dp = 26.1e-6 / (mu*ham.Rv*T)
        
        return dp
    
    def set_perm_liquid(self, method, **kwargs):
        """
        Register the liquid permeability [s] (optional)
        
        :method: String indicating how the permeability is defined. Available
                 methods are 'durner' and 'exp'
        :kwargs: Dictionary of parameters for the definition of the profile
                 See the wiki of the hamopy project for instructions regarding
                 the definition of this object.
        """
        
        self.kl_method = method
        
        if method == 'durner':
        
            self.K_sat  = kwargs['K_sat']
            self.tau_kl = kwargs['tau']
            if isinstance(kwargs['l'], float):
                self.l_kl       = np.array([kwargs['l']])
                self.alpha_kl   = np.array([kwargs['alpha']])
                self.m_kl       = np.array([kwargs['m']])
            else:
                self.l_kl       = np.array(kwargs['l'])
                self.alpha_kl   = np.array(kwargs['alpha'])
                self.m_kl       = np.array(kwargs['m'])
                
        elif method == 'exp':
            
            self.a_kl = np.array(kwargs['a'])
            
        elif method == 'exp2':
            
            self.a_kl = np.array(kwargs['a'])
            self.w0_kl = np.array(kwargs['w0'])
            
        elif method == 'interp':
            
            self.logpsuc = np.log10( -np.array(kwargs['PC']) ) # positif
            self.logkl   = np.log10(  np.array(kwargs['KL']) )
            self.kl_f = UnivariateSpline(self.logpsuc, self.logkl, k=1)

            
    def k_l(self, p_c):
        """
        Liquid permeability [s]
        
        The capillary pressure p_c is required
        """
        
        if not hasattr(self, 'kl_method'):
            # Returns zero if no water permeability has been defined
            return np.zeros(np.shape(p_c))
        
        elif self.kl_method == 'durner':
            # Mono- or multimodal Durner law
            numerateur   = np.zeros(np.shape(p_c))
            s            = np.zeros(np.shape(p_c))
            denominateur = np.sum(self.l_kl * self.alpha_kl)
            n = 1./(1-self.m_kl)
            for i in range(np.size(self.l_kl)):
                s_i = (1+(self.alpha_kl[i]*abs(p_c))**n[i])**(-self.m_kl[i])
                s  += self.l_kl[i] * s_i
                numerateur += self.l_kl[i]*self.alpha_kl[i] * \
                (1- (1-s_i**(1/self.m_kl[i]))**self.m_kl[i] )
                
            return self.K_sat * s**self.tau_kl * (numerateur/denominateur)**2
            
        elif self.kl_method == 'exp':
            # Exponential law (used in some Hamstad benchmarks)
            w = self.w(p_c)
            lnkl = np.zeros(np.shape(p_c))
            for i in range(len(self.a_kl)):
                lnkl += self.a_kl[i] * (w/ham.rho_liq)**i
            
            return np.exp(lnkl)
            
        elif self.kl_method == 'exp2':
            
            w = self.w(p_c)
            lnkl = np.zeros(np.shape(p_c))
            for i in range(len(self.a_kl)):
                lnkl += self.a_kl[i] * (w - self.w0_kl)**i
            
            return np.exp(lnkl)
            
        elif self.kl_method == 'interp':
            
            log10psuc = np.log10(-p_c)
            kl = 10 ** self.kl_f( log10psuc.ravel() )
            
            return kl.reshape(np.shape(p_c))

    def set_perm_air(self, k_air):
        """
        Register the air permeability of the dry material [m2].
        One argument (float) needed
        """
        
        self.k_air = k_air


class FiniteElement(object):
    """
    The :class:`FiniteElement` object stores all information regarding the node
    coordinates, size and integration points of one finite element of the mesh.
    
    A list of :class:`FiniteElement` is created upon initialisation of the mesh
    and stored as the mesh.e object.
    
    Example : mesh.e[4].GL_W is an array containing the weights of the 4th
    element's integration points for the Gauss-Legendre integration scheme.
    """
    
    def __init__(self, x_min, x_max, indices, nbr_nodes):
        
        self.x       = np.round( np.linspace(x_min, x_max, 4) , 8)
        self.l       = x_max - x_min
        self.indices = indices
        
        # Calculate all variables regarding the numerical integration scheme
        self.gauss_legendre(x_min, x_max)
        
        # Register where the element is located in a square matrix
        self.global_indices(indices, nbr_nodes)

    def gauss_legendre(self, x_min, x_max):
        """
        Calculate all variables regarding the numerical integration scheme
        
        Input: edge coordinates of the element
        """
        
        # Local coordinate and weight of each integration point of the element
        GL_x = 1./35 * np.array([-sqrt(525+70*sqrt(30)), -sqrt(525-70*sqrt(30)), sqrt(525-70*sqrt(30)), sqrt(525+70*sqrt(30))])
        GL_w = 1./36 * np.array([18-sqrt(30)           , 18+sqrt(30)           , 18+sqrt(30)          , 18-sqrt(30)          ])
        # Global coordinate and weight of each integration point of the element
        self.GL_X = (x_max-x_min)/2 * GL_x + (x_max+x_min)/2
        self.GL_W = (x_max-x_min)/2 * GL_w
        
        # Interpolation functions evaluated at each integration point
        # GL_N[i,j] = value of the i-th interp. function at the j-th integ. point
        # GL_B[i,j] = derivative of the i-th interp. function at the j-th integ. point
        s  = (1+GL_x) / 2
        L1 = s
        L2 = 1 - s
        N1 = L1/2 * (3*L1-1) * (3*L1-2)
        N2 = 9./2 * L1 * L2  * (3*L1-1)
        N3 = 9./2 * L1 * L2  * (3*L2-1)
        N4 = L2/2 * (3*L2-1) * (3*L2-2)
        dN1dx = 1./(x_max-x_min) * 0.5 * (-27*s**2 + 36*s - 11)
        dN2dx = 1./(x_max-x_min) * 4.5 * (  9*s**2 - 10*s + 2 )
        dN3dx = 1./(x_max-x_min) * 4.5 * ( -9*s**2 +  8*s - 1 )
        dN4dx = 1./(x_max-x_min) * 0.5 * ( 27*s**2 - 18*s + 2 )
        GL_N = np.array([N1, N2, N3, N4])
        GL_B = np.array([dN1dx, dN2dx, dN3dx, dN4dx])
        
        # Interpolation matrices evaluated at each integration point
        # GL_NtN[:,:,i] = N*tN matrix evaluated at the i-th integration point
        # GL_BtB[:,:,i] = B*tB matrix evaluated at the i-th integration point
        NtN = np.zeros((4,4,4), dtype=float)
        BtB = np.zeros((4,4,4), dtype=float)
        BtN = np.zeros((4,4,4), dtype=float)
        for i in range(4):
            NtN[:,:,i] = np.outer(GL_N[:,i], GL_N[:,i])
            BtB[:,:,i] = np.outer(GL_B[:,i], GL_B[:,i])
            BtN[:,:,i] = np.outer(GL_B[:,i], GL_N[:,i])
        
        # This is what we save
        self.GL_N   = GL_N
        self.GL_B   = GL_B
        self.GL_NtN = NtN
        self.GL_BtB = BtB
        self.GL_BtN = BtN
    
    def global_indices(self, indices, nbr_nodes):
        """
        Register where the element is located in a square matrix
        
        This is later used for the assembling of the global matrices of the
        linear system
        """
        a = nbr_nodes * np.tile(indices, (4,1))
        b = np.tile(indices, (4,1)).T
        self.indices_globaux = (a + b).ravel()
        
        self.indices_row = np.tile(indices, 4)
        self.indices_col = np.sort(self.indices_row)


class Mesh(object):
    """
    Creation of the :class:`Mesh` object 
    
    Three lists of equal sizes are required as arguments :
    
    :materials: List of :class:`Material` : what each layer is made of
    :sizes: List of floats : width of each material layer
    :nbr_elements: List of integers : number of finite elements in each layer
    """
    def __init__(self, materials, sizes, nbr_elements):
        
        # Mesh initialisation
        self.materials = materials
        self.sizes     = sizes
        self.nbr_elem  = nbr_elements
        self.nbr_nodes = 3 * sum(nbr_elements) + 1

        # List of material indices of each finite element
        elem_material = np.array([])
        # List of the width of each finite element
        elem_size  = np.array([])
        for i in range(np.size(sizes)):
            foo = i * np.ones(nbr_elements[i])
            bar = sizes[i] / nbr_elements[i] * np.ones(nbr_elements[i])
            elem_material = np.concatenate((elem_material, foo), axis=0)
            elem_size     = np.concatenate((elem_size,     bar), axis=0)
        
        # Material index and size of each element of the mesh
        self.elem_material = [int(i) for i in elem_material]
        self.elem_size     = elem_size
        
        # Construction of the list of finite elements
        elements = []
        x_nodes  = np.zeros(3*sum(nbr_elements) + 1)
        for e in range(sum(self.nbr_elem)):
            x_min   = sum(self.elem_size[:e])
            x_max   = sum(self.elem_size[:e+1])
            indices = [3*e, 3*e+1, 3*e+2, 3*e+3]
            new_element = FiniteElement(x_min, x_max, indices, self.nbr_nodes)
            elements.append(new_element)
            x_nodes[indices] = new_element.x
        
        # The list of all elements and nodes is saved
        self.e = elements
        self.x = x_nodes
        
        # All information regarding the Gauss-Legendre integration scheme
        self.GL_X   = np.array( [self.e[_].GL_X   for _ in range(sum(self.nbr_elem))] )
        self.GL_W   = np.array( [self.e[_].GL_W   for _ in range(sum(self.nbr_elem))] )
        self.GL_N   = np.array( [self.e[_].GL_N   for _ in range(sum(self.nbr_elem))] )
        self.GL_B   = np.array( [self.e[_].GL_B   for _ in range(sum(self.nbr_elem))] )
        self.GL_NtN = np.array( [self.e[_].GL_NtN for _ in range(sum(self.nbr_elem))] )
        self.GL_BtB = np.array( [self.e[_].GL_BtB for _ in range(sum(self.nbr_elem))] )
        self.GL_BtN = np.array( [self.e[_].GL_BtN for _ in range(sum(self.nbr_elem))] )
        
        # These arrays are later used for the assembling of elementary matrices
        self.indices_row = np.concatenate(([self.e[i].indices_row for i in range(sum(self.nbr_elem))]))
        self.indices_col = np.concatenate(([self.e[i].indices_col for i in range(sum(self.nbr_elem))]))
        
        """
        # I don't remember what this is for - probably obsolete
        foo = 3*np.array( range(np.sum(self.nbr_elem)) )
        foo = foo[:, np.newaxis]
        self.indices_s = foo + np.array([0, 1, 2, 3])
        """
        
        # Equivalent air permeability of the entire wall
        if np.any( [materials[_].k_air == 0 for _ in range(len(materials))] ):
            self.C_air = 0
        else:
            D = sum(self.sizes)
            K_eq = D / sum([self.sizes[_] / self.materials[_].k_air for _ in range(len(self.materials))])
            self.C_air = K_eq / D * ham.rho_air / ham.mu_air
        # Matrice int(BtN) sur l'ensemble du systeme
        bar = np.sum( self.GL_W.reshape((sum(self.nbr_elem),1,1,4)) * self.GL_BtN, axis=-1)
        self.air_matrix = coo_matrix((bar.ravel(), (self.indices_row,self.indices_col)))
        
    def replace_materials(self, mat):
        """
        Small method to replace one of the materials without remeshing
        
        Allows saving time during sensitivity analyses and such
        """
        self.materials = mat
        
    def find_material(self, x):
        """
        Finds in which material the coordinate x belongs
        """
        
        a = np.cumsum(self.sizes)
        b = np.where(x<a)
        c = b[0][0]
        
        return self.materials[c]
        
    def system_matrices(self, P, T):
        """
        Method of :class:`Mesh` called by hamopy.iteration
        
        Calculates the capacity and conductivity matrices C and K of the global
        linearised system for the updating of P and T. This is where the
        storage and transport coefficients of the transport equations are
        involved.
        
        C and K are sparse matrices
        """
        
        # Value of P and T at the GL integration points
        f_p = interp1d(self.x, np.log10(-P))
        f_t = interp1d(self.x, T)
        GL_P = -10**f_p(self.GL_X)
        GL_T = f_t(self.GL_X)
        
        c_mm, c_mh, c_hh, c_hm = (np.zeros(np.shape(self.GL_X)) for _ in range(4))
        k_mm, k_mh, k_hh, k_hm = (np.zeros(np.shape(self.GL_X)) for _ in range(4))
        
        # Loop over the subdomains (one or several per material)
        for sd in range(np.size(self.materials)):
            
            # Pick which elements of the mesh belong to the current subdomain
            m = self.materials[sd]
            mask = (np.array(self.elem_material) == sd)
            
            p = GL_P[mask]
            t = GL_T[mask]
            
            # Capacities
            c_mm[mask] = m.c_w(p)
            c_mh[mask] = np.zeros(np.shape(p))
            c_hh[mask] = m.cp(t) * m.rho + ham.cp_liq * m.w(p)
            c_hm[mask] = ham.cp_liq * (t-273.15) * m.c_w(p)
            
            # Permeabilities
            G_P  = ham.p_v(p, t) / (ham.rho_liq * ham.Rv * t)
            G_T  = G_P / t * (ham.rho_liq*ham.l_lv - p)
            k_mm[mask] = G_P * m.delta_p(p,t) + m.k_l(p)
            k_mh[mask] = G_T * m.delta_p(p,t)
            k_hh[mask] = ham.l_lv * m.delta_p(p,t) * G_T + m.conduc(p,t)
            k_hm[mask] = ham.l_lv * m.delta_p(p,t) * G_P + ham.cp_liq * m.k_l(p) * (t-273.15)
        
        # Numerical integration of the elementary matrices
        E = np.sum(self.nbr_elem)
        def integration(c, NtN):
            return np.sum( (self.GL_W*c).reshape((E,1,1,4)) * NtN, axis=-1)
        C_mm_e, C_mh_e, C_hh_e, C_hm_e = \
            [integration(i, self.GL_NtN) for i in [c_mm, c_mh, c_hh, c_hm]]
        K_mm_e, K_mh_e, K_hh_e, K_hm_e = \
            [integration(i, self.GL_BtB) for i in [k_mm, k_mh, k_hh, k_hm]]
            
        # Assembling of the elementary matrices
        def assemblage(C):
            return coo_matrix((C.ravel(), (self.indices_row,self.indices_col)))
        C_mm, C_mh, C_hh, C_hm = \
            [assemblage(I) for I in [C_mm_e, C_mh_e, C_hh_e, C_hm_e] ]
        K_mm, K_mh, K_hh, K_hm = \
            [assemblage(I) for I in [K_mm_e, K_mh_e, K_hh_e, K_hm_e] ]
        
        # Construction of the global matrices
        C = vstack([ hstack([C_mm, C_mh]), hstack([C_hm, C_hh]) ])
        K = vstack([ hstack([K_mm, K_mh]), hstack([K_hm, K_hh]) ])
        
        return C, K
        
    
    def system_airflux(self, P, T, clim, time):
        """
        Method of :class:`Mesh` called by hamopy.iteration
        
        Construction of the source terms, added to the moisture and heat
        transfer equations, caused by air transfer in the wall
        """
        
        # Air flow [kg/m2s]
        g_air = self.C_air * (clim[0].P_air(time) - clim[1].P_air(time))
        
        omega_m = ham.p_v(P, T) * g_air / (ham.Rv * T * ham.rho_air)
        omega_h = ham.l_lv * omega_m + ham.cp_air * (T-273.15) * g_air
        
        # Hygric and thermal source terms
        QM = self.air_matrix * omega_m
        QH = self.air_matrix * omega_h
        
        return np.concatenate((QM, QH))

    def system_matrices_thermo(self, T):
        """
        Method of :class:`Mesh` called by hamopy.iteration_thermo
        
        Calculates the capacity and conductivity matrices C and K of the global
        linearised system for the updating of T. This is where the storage and
        transport coefficients of the transport equations are involved.
        
        C and K are sparse matrices
        """
        
        # Value of T at the GL integration points
        f_t = interp1d(self.x, T)
        GL_T = f_t(self.GL_X)
        
        c_hh, k_hh = (np.zeros(np.shape(self.GL_X)) for _ in range(2))
        
        # Loop over the subdomains (one or several per material)
        for sd in range(np.size(self.materials)):
            
            # Pick which elements of the mesh belong to the current subdomain
            m = self.materials[sd]
            mask = (np.array(self.elem_material) == sd)

            # Capacity and conductivity
            c_hh[mask] = m.cp(GL_T[mask]) * m.rho
            k_hh[mask] = m.conduc(GL_T[mask])
        
        # Numerical integration of the elementary matrices
        E = np.sum(self.nbr_elem)
        def integration(c, NtN):
            return np.sum( (self.GL_W*c).reshape((E,1,1,4)) * NtN, axis=-1)
        C_hh_e = integration(c_hh, self.GL_NtN)
        K_hh_e = integration(k_hh, self.GL_BtB)
            
        # Assembling of the elementary matrices
        def assemblage(C):
            return coo_matrix((C.ravel(), (self.indices_row,self.indices_col)))

        # Construction of the global matrices
        C = assemblage(C_hh_e)
        K = assemblage(K_hh_e)
        
        return C, K
        
    def system_matrices_hygro(self, P):
        """
        Method of :class:`Mesh` called by hamopy.iteration_thermo
        
        Calculates the capacity and conductivity matrices C and K of the global
        linearised system for the updating of T. This is where the storage and
        transport coefficients of the transport equations are involved.
        
        C and K are sparse matrices
        """
        
        # Value of T at the GL integration points
        f_p = interp1d(self.x, np.log10(-P))
        GL_P = -10**f_p(self.GL_X)
        
        c_mm, k_mm = (np.zeros(np.shape(self.GL_X)) for _ in range(2))
        
        # Loop over the subdomains (one or several per material)
        for sd in range(np.size(self.materials)):
            
            # Pick which elements of the mesh belong to the current subdomain
            m = self.materials[sd]
            mask = (np.array(self.elem_material) == sd)

            p = GL_P[mask]
            t = 293.15
            
            # Capacity
            c_mm[mask] = m.c_w(p)
            # Permeability
            G_P  = ham.p_v(p, t) / (ham.rho_liq * ham.Rv * t)
            k_mm[mask] = G_P * m.delta_p(p,t) + m.k_l(p)
        
        # Numerical integration of the elementary matrices
        E = np.sum(self.nbr_elem)
        def integration(c, NtN):
            return np.sum( (self.GL_W*c).reshape((E,1,1,4)) * NtN, axis=-1)
        C_mm_e = integration(c_mm, self.GL_NtN)
        K_mm_e = integration(k_mm, self.GL_BtB)
            
        # Assembling of the elementary matrices
        def assemblage(C):
            return coo_matrix((C.ravel(), (self.indices_row,self.indices_col)))

        # Construction of the global matrices
        C = assemblage(C_mm_e)
        K = assemblage(K_mm_e)
        
        return C, K
        
    def system_boundary(self, P, T, clim, t):
        """
        Method of :class:`Mesh` called by hamopy.iteration
        
        Calculates the vector F of boundary conditions, and the matrix dFdU
        involved in the Newton-Raphson iteration scheme
        """
        
        N    = self.nbr_nodes
        F    = np.zeros(2*N)
        data = []
        
        g_air = self.C_air * (clim[0].P_air(t) - clim[1].P_air(t))
        
        # Loop over each boundary
        ind = [0, N - 1]
        for i in range(2):
            
            # Moisture and rain
            E  = clim[i].h_m(t) * ( clim[i].p_v(t) - ham.p_v(P[ind[i]], T[ind[i]]) )
            R  = clim[i].g_l(t)
            
            # Heat
            H  = clim[i].h_t(t) * ( clim[i].T_eq(t) - T[ind[i]] )
            CR = R * ham.cp_liq * ( clim[i].T(t) - 273.15 )
            LE = E * ham.l_lv
            
            # Air
            A  = g_air * (-1)**i * (clim[i].p_v(t)/clim[i].T(t) - ham.p_v(P[ind[i]], T[ind[i]])/T[ind[i]] ) / (ham.Rv * ham.rho_air)
            HA = A * ham.l_lv + ham.cp_air * (clim[i].T(t)-T[ind[i]]) * g_air * (-1)**i
            
            # Construction of F
            F[ind[i]]   = E + A + R
            F[ind[i]+N] = H + HA + CR + LE
            
            # Construction of dFdU
            data.extend([ -clim[i].h_m(t) * ham.p_v(P[ind[i]],T[ind[i]]) / ( ham.rho_liq*ham.Rv*T[ind[i]] ),
                          clim[i].h_m(t) * ham.p_v(P[ind[i]],T[ind[i]]) / ( ham.Rv*T[ind[i]]**2) * ( P[ind[i]]/ham.rho_liq - ham.l_lv ) ])
            data.extend([ ham.l_lv * data[-2],
                          ham.l_lv * data[-1] - clim[i].h_t(t) ])
            
        row  = [0, 0, N, N, N-1,   N-1, 2*N-1, 2*N-1]
        col  = [0, N, 0, N, N-1, 2*N-1,   N-1, 2*N-1]
        dFdU = coo_matrix( (data, (row,col)) )
        
        return F, dFdU
        
    def system_boundary_thermo(self, T, clim, t):
        """
        Method of :class:`Mesh` called by hamopy.iteration_thermo
        
        Calculates the vector F of boundary conditions, and the matrix dFdU
        involved in the Newton-Raphson iteration scheme
        """
        
        N    = self.nbr_nodes
        F    = np.zeros(N)
        data = []
        
        # Loop over each boundary
        ind = [0, N - 1]
        for i in range(2):
            
            F[ind[i]] = clim[i].h_t(t) * ( clim[i].T_eq(t) - T[ind[i]] )
            data.extend([-clim[i].h_t(t)])
            
        row = [0, N-1]
        col = [0, N-1]
        dFdU = coo_matrix( (data, (row,col)) )
        
        return F, dFdU
        
    def system_boundary_hygro(self, P, clim, t):
        """
        Method of :class:`Mesh` called by hamopy.iteration_thermo
        
        Calculates the vector F of boundary conditions, and the matrix dFdU
        involved in the Newton-Raphson iteration scheme
        """
        
        N    = self.nbr_nodes
        F    = np.zeros(N)
        data = []
        
        # Loop over each boundary
        ind = [0, N - 1]
        for i in range(2):
            
            E  = clim[i].h_m(t) * ( clim[i].p_v(t) - ham.p_v(P[ind[i]], 293.15) )
            R  = clim[i].g_l(t)
            F[ind[i]]   = E + R
            data.extend([ -clim[i].h_m(t) * ham.p_v(P[ind[i]],293.15) / ( ham.rho_liq*ham.Rv*293.15 ) ])
            
        row = [0, N-1]
        col = [0, N-1]
        dFdU = coo_matrix( (data, (row,col)) )
        
        return F, dFdU
        
    def system_conserv(self, P, T):
        """        
        Method of :class:`Mesh` called by hamopy.iteration
        
        Calculates the S matrix for ensuring heat and mass conservation
        between time steps
        """
        
        # Value of T at the GL integration points
        f_p = interp1d(self.x, np.log10(-P))
        f_t = interp1d(self.x, T)
        GL_P = -10**f_p(self.GL_X)
        GL_T = f_t(self.GL_X)
        
        s_m, s_h = (np.zeros(np.shape(self.GL_X)) for _ in range(2))
        
        # Loop over the subdomains (one or several per material)
        for sd in range(np.size(self.materials)):
            
            # Pick which elements of the mesh belong to the current subdomain
            m = self.materials[sd]
            mask = (np.array(self.elem_material) == sd)
            
            p = GL_P[mask]
            t = GL_T[mask]
            
            s_m[mask] = m.w(p)
            s_h[mask] = (t-273.15) * (m.rho*m.cp(t) + ham.cp_liq*m.w(p))
            
        # Integration of the elementary vectors
        E = np.sum(self.nbr_elem)
        def integration(s):
            return np.sum( (self.GL_W*s).reshape((E,1,4)) * self.GL_N, axis=-1)
        [S_m_e, S_h_e] = [integration(_) for _ in [s_m, s_h]]
        
        # Assembling of the elementary vectors
        def assemblage_vecteur(S_e):
            S = np.zeros(self.nbr_nodes)
            for i in range(sum(self.nbr_elem)):
                S[self.e[i].indices] += S_e[i]
            return S
        S_m, S_h = [assemblage_vecteur(_) for _ in [S_m_e, S_h_e]]
        
        # Global S vector
        S = np.concatenate((S_m,S_h), axis = 0)
        
        return S
        
    """
    def assemblage_vecteur(self, S_elem):
        
        S = np.zeros(self.nbr_nodes)
        
        for i in range(sum(self.nbr_elem)):
            S[self.e[i].indices] += getattr(self.e[i], S_elem)
            
        return S
    """

class Boundary:
    """
    The :class:`Boundary` class includes methods to call the value of all field
    variables as functions of the time of simulation. The following arguments
    must be given at the initialisation of a new boundary:
    
    :method: String denoting the type of BC. Must be either 'Fourier' to define
             a transfer coefficient BC, or 'Dirichlet' to impose T and P
    :kwargs: Dictionary containing all numerical values of the BC
    
    If any boundary value varies in time, the kwargs dict must include a 'file'
    key pointing to the location of the .txt file, and optionally a 'delimiter'
    key (default is tab).
    
    A constant value at the boundary is given with a float associated to the
    corresponding key of the kwargs dictionary. A variable value is given with
    the string, which is the header of the value in the file denoted by 'file'.
    The following is a list of keys and their expected associated values:
    
    :T: Temperature
    :T_eq: (optional) equivalent temperature including effects of radiation
    :HR: Relative humidity
    :p_v: Vapor pressure. Either HR or p_v must be given.
    :h_t: (optional) Heat transfer coefficient. Default is 5 [W/m2K]. This
          value is not used in case of Dirichlet BC
    :h_m: (optional) Moisture transfer coefficient. Default is
          :math:`7.45 \times 10^{-9} h_t`
    :g_l: (optional) Rain [kg/(m2.s)]
    :P_air: (optional) Air pressure [Pa]
    """
    def __init__(self, method='Fourier', **kwargs):
        
        # kwargs must have a 'T' key, and either 'HR' or 'p_v'
        ok = 'T' in kwargs.keys() and 'HR' in kwargs.keys() or 'p_v' in kwargs.keys()
        if not ok:
            raise Exception('T and HR or p_v must be given in the boundary definition')
        
        if 'file' in kwargs.keys():

            if 'time' in kwargs.keys():
                if 'delimiter' not in kwargs.keys():
                    data0 = pd.read_csv(kwargs['file'], delimiter='\t')
                    del kwargs['file']
                else:
                    data0 = pd.read_csv(kwargs['file'], delimiter=kwargs['delimiter'])
                    del kwargs['file'], kwargs['delimiter']
            else:
                raise Exception('Time argument missing in the boundary definition')
            
            # If any value is variable, all other values are saved as constant
            # arrays
            kwargs_new = {}
            taille_du_fichier = np.shape( data0[kwargs['time']] )
            for key in kwargs:
                if type(kwargs[key]) == str:
                    kwargs_new[key] = np.array( data0[kwargs[key]] )
                else:
                    kwargs_new[key] = kwargs[key] * np.ones(taille_du_fichier)
               
            # If T has been given in C, it is switched to K 
            if kwargs_new['T'].mean() < 200:
                kwargs_new['T'] += 273.15
            if 'T_eq' in kwargs_new.keys():
                kwargs_new['T_eq'] += 273.15

        else:
            # All values are constant or given in the initial dictionary
            kwargs_new = kwargs
        
            # If T has been given in C, it is switched to K      
            if kwargs_new['T'] < 200:
                kwargs_new['T'] += 273.15
                if 'T_eq' in kwargs_new.keys():
                    kwargs_new['T_eq'] += 273.15

        
        # If no equivalent temperature is given, it is set to that of air
        if 'T_eq' not in kwargs_new.keys():
            kwargs_new['T_eq'] = kwargs_new['T']
            
        # If no rain is given, the value of g_l is 0
        if 'g_l' not in kwargs_new.keys():
            kwargs_new['g_l'] = 0.
        
        # If HR is not given, p_v is expected
        if 'HR' not in kwargs_new.keys():
            kwargs_new['HR'] = kwargs_new['p_v'] / ham.p_sat(kwargs_new['T'])
        
        # If p_v is not given, HR is expected
        if 'p_v' not in kwargs_new.keys():
            kwargs_new['p_v'] = kwargs_new['HR'] * ham.p_sat(kwargs_new['T'])
        
        # If h_t is not given, default value is 5 W/m2K
        if 'h_t' not in kwargs_new.keys():
            kwargs_new['h_t'] = 5.
        
        # If h_m is not given, default value is 7.45e-9 * h_t
        if 'h_m' not in kwargs_new.keys():
            kwargs_new['h_m'] = 7.45e-9 * kwargs_new['h_t']
        
        # If no air pressure is given, the value of P_air is 0
        if 'P_air' not in kwargs_new.keys():
            kwargs_new['P_air'] = 0.
         
        if method == 'Dirichlet':
            kwargs_new['h_t'] = 1000.
            kwargs_new['h_m'] = 1e-4
        
        self.type = method
        self.data = kwargs_new
        self.constant = not 'time' in self.data.keys()
        
    def has_constant(self, var):
        return isinstance(self.data[var], float) or isinstance(self.data[var], int)

    def T(self, t):
        """ Outside temperature [K] """
        if self.has_constant('T'):
            return self.data['T']
        else:
            out = self.data['T'][0]
            f = interp1d(self.data['time'], self.data['T'], bounds_error=False, fill_value=out)
            return f(t)
    
    def HR(self, t):
        """ Outside relative humidity [-] """
        if self.has_constant('HR'):
            return self.data['HR']
        else:
            out = self.data['HR'][0]
            f = interp1d(self.data['time'], self.data['HR'], bounds_error=False, fill_value=out)
            return f(t)
    
    def T_eq(self, t):
        """ Outside equivalent temperature [K] """
        if self.has_constant('T_eq'):
            return self.data['T_eq']
        else:
            out = self.data['T_eq'][0]
            f = interp1d(self.data['time'], self.data['T_eq'], bounds_error=False, fill_value=out)
            return f(t)
    
    def p_v(self, t):
        """ Outside vapor pressure [Pa] """
        if self.has_constant('p_v'):
            return self.data['p_v']
        else:
            out = self.data['p_v'][0]
            f = interp1d(self.data['time'], self.data['p_v'], bounds_error=False, fill_value=out)
            return f(t)
    
    def p_c(self, t):
        return ham.p_c(self.HR(t), self.T(t))
    
    def h_t(self, t):
        """ Heat transfer coefficient [W/(m2K)]"""
        if self.has_constant('h_t'):
            return self.data['h_t']
        else:
            out = self.data['h_t'][0]
            f = interp1d(self.data['time'], self.data['h_t'], bounds_error=False, fill_value=out)
            return f(t)
    
    def h_m(self, t):
        """ Moisture transfer coefficient [s/m] """
        if self.has_constant('h_m'):
            return self.data['h_m']
        else:
            out = self.data['h_m'][0]
            f = interp1d(self.data['time'], self.data['h_m'], bounds_error=False, fill_value=out)
            return f(t)
    
    def g_l(self, t):
        """ Rain [kg/(m2s)] """
        if self.has_constant('g_l'):
            return self.data['g_l']
        else:
            out = self.data['g_l'][0]
            f = interp1d(self.data['time'], self.data['g_l'], bounds_error=False, fill_value=out)
            return f(t)
    
    def P_air(self, t):
        """ Air pressure [Pa] """
        if self.has_constant('P_air'):
            return self.data['P_air']
        else:
            out = self.data['P_air'][0]
            f = interp1d(self.data['time'], self.data['P_air'], bounds_error=False, fill_value=out)
            return f(t)

class Time:
    """
    The :class:`Time` class stores information regarding the duration of the
    simulation and the discretisation in time. Two arguments are demanded upon
    instantiation of the object:
    
    :method: String, either 'constant' or 'variable', indicating whether the
             time step must adapt to the speed of convergence
    :kwargs: Dictionary containing all numerical values
    
    The kwargs dict must contain the following keys, pointing to their
    respective values:
    
    :delta_t: Time step size [s]
    :t_max: Total simulation time [s]
    :iter_max: (required if method=='variable') Maximum number of iterations
               per time step
    :delta_min: (required if method=='variable') Minimum allowed time step,
                under which the simulation fails
    :delta_max: (required if method=='variable') Maximum time step size
    """
    
    def __init__(self, method, **kwargs):
        self.method = method
        self.vec    = [0]
        self.t      = kwargs['delta_t']
        self.delta  = kwargs['delta_t']
        self.end    = kwargs['t_max']
        self.stop   = False
        
        if method == 'variable':
            
            if 'iter_max' in kwargs.keys():
                self.iter_max = kwargs['iter_max']
            else:
                self.iter_max = 12
                
            if 'delta_min' in kwargs.keys():
                self.delta_min = kwargs['delta_min']
            else:
                self.delta_min = 1e-3
            
            if 'delta_max' in kwargs.keys():
                self.delta_max = kwargs['delta_max']
            else:
                self.delta_max = 900
            
        elif method == 'constant':
            pass
        
        else:
            raise Exception('Time method must be either constant or variable')
        
        
    def next_step(self, nbr_iter):        
        
        if self.method == 'variable':
            self.delta *= np.min( [self.iter_max/(2.*nbr_iter), 2.])
            self.delta  = np.min( [self.delta, self.delta_max] )
            if self.delta < self.delta_min:
                self.stop = True
        
        self.t += self.delta
    

    def try_again(self, value, logobj):
        """
        Called by calcul() if the maximum number of iterations has been reached
        in the current time step.
        """
        
        if self.method == 'constant':
            self.stop = True
            
        elif self.method == 'variable':
            if logobj is not None:
                logobj.write("Convergence not reached: time step reduced by %i \n" % value)
            self.delta   /= value
            self.t       -= (value-1.) * self.delta
            self.vec[-1]  = self.t
            if self.delta < self.delta_min:
                self.stop = True
            elif logobj is not None:
                logobj.write("t = %s \n" % self.t)
                
    def restart(self):
        self.vec = [0]
        self.t   = self.delta
