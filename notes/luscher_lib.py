import numpy as np
import scipy
import scipy.integrate

class Kcot_luscher(object):
    """
    Calculate the scattering phase shifts using finite volume formula
    """
    def __init__(self, L, m_1, m_2, m_pi, d=[0,0,0]):
        """
        Input parameters
        L : int.
            length of the box in lattice unit
        m_1 (m_2) : float
            mass of the particle in lattice unit
        m_pi : float
            pion mass (reference scale) in lattice unit
        d : list
            boosted momentum ex. d=[1,0,0], d=[2,0,0]
        """
        self.L = L
        self.m_1  = m_2 # in lattice unit
        self.m_2  = m_1 # in lattice unit
        self.m_pi = m_pi # in lattice unit

        self.d = d
        self.z00 = Luscher_zeta(L, d)

    def gamma(self, k2, d, Lval):
        p = np.sqrt((np.dot(d,d))) * 2.0*np.pi/float(Lval)
        Wcm = np.sqrt(k2 + self.m_1**2) +  np.sqrt(k2 + self.m_1**2)
        Wl = np.sqrt(Wcm**2 + p**2)
        v = p/Wl
        gamma_val = 1.0/np.sqrt(1.0 - v**2)
        return gamma_val

    def kcot(self, k2):
        """
        Return kcot/mpi

        Input parameter 
        k2 : float (in lattice unit)
        """
        if self.L == 'inf':
            return ( -  np.sqrt(-k2))

        if self.d == [0,0,0]:
            gamma_val = 1.0
            z00 = self.z00.z00(k2)
        else:
            gamma_val = self.gamma(k2, self.d, self.L)
            z00 = self.z00.z00d(k2, gamma)
        return 2.0*z00/(self.L*np.sqrt(np.pi))/gamma_val/self.m_pi

    def kcot_from_dE(self, deltaE):
        """
        Return kcot/mpi

        Input parameter 
        delta_E : float (in lattice unit)
            energy shift delta E
        """
        k2 = self.tok2(deltaE)
        if self.L == 'inf':
            return ( -  np.sqrt(-k2))

        if self.d == [0,0,0]:
            gamma_val = 1.0
            z00 = self.z00.z00(k2)
        else:
            gamma_val = self.gamma(k2, self.d, self.L)
            z00 = self.z00.z00d(k2, gamma)
        return 2.0*z00/(self.L*np.sqrt(np.pi))/gamma_val/self.m_pi


    def tok2(self, deltaE):
        Cval = lambda dE: 0.5*dE**2 + (self.m_1 + self.m_2)*dE + self.m_1*self.m_2
        conv_k2 = lambda dE: (Cval(dE)**2 - self.m_1**2 * self.m_2**2)/(self.m_1**2 + self.m_2**2 + 2.0*Cval(dE))
        return conv_k2(deltaE)


class Luscher_zeta(object):
    """
    Evaluate zeta function in Luscher's formula for S-wave state (with boosted frame)
    based on Ref. M. Luscher, Nucl. Phys. B354 (1991) 531-578.
    """
    def __init__(self, L, d = [0,0,0], Nmax=20, lcut=3):
        self.L = L

        self.d = d
        self.r_vec = np.array([[ix,iy,iz] for iz in range(-Nmax//2,Nmax//2)
            for iy in range(-Nmax//2,Nmax//2) for ix in range(-Nmax//2,Nmax//2)])
        self.lcut = lcut

        n2_list = np.array([ix**2 + iy**2 + iz**2 for iz in range(-Nmax//2, Nmax//2)
                            for iy in range(-Nmax//2, Nmax//2) for ix in range(-Nmax//2, Nmax//2)])

        self.n2_small = n2_list[n2_list < lcut**2]
        self.n2_large = n2_list[n2_list >= lcut**2]
        self.n2_ne0   = n2_list[n2_list >= 1] 

    def z00(self, k2, tmin = 0.0e0, tcut = 1.0e0, tmax = np.inf):
        q2 = k2 * self.L**2/((2.0*np.pi)**2)

        integ_I_low = lambda t, q2: 0.5*np.pi * (np.exp(t*q2) - 1.0)/(t**1.5) \
                + 0.5*np.pi * np.sum(np.exp(t*q2 - np.pi**2 * self.n2_ne0/t))/(t**1.5) \
                - (1.0/(np.sqrt(4.0*np.pi))) * np.sum(np.exp((q2 - self.n2_small) * t))
                
        I_low, I_low_err = scipy.integrate.quad(integ_I_low,tmin,tcut,args=(q2))
        
        integ_I_high = lambda t, q2: (1.0/(np.sqrt(4.0*np.pi))) \
                    * np.sum(np.exp((q2 - self.n2_large) * t)) - 0.5 * np.pi /(t**1.5)

        I_high, I_high_err = scipy.integrate.quad(integ_I_high,tcut,tmax, args=(q2))

        z_0 = (1.0/np.sqrt(4.0*np.pi)) * np.sum(1.0/(self.n2_small - q2))

        return z_0 + I_low + I_high

    def z00d(self, k2, gamma = 1.0, tmin = 0.0e0, tcut = 1.0e0, tmax = np.inf):
        q2 = k2 * self.L**2/((2.0*np.pi)**2)
        
        gv = np.array([gamma, 1.0, 1.0])
        gvinv = np.array([1.0/gamma, 1.0, 1.0])

        r2_list = np.sum(((self.r_vec[:,:] + 0.5*np.array(self.d))*gvinv)**2, axis=1)
     
        r2_small = r2_list[r2_list < self.lcut**2]
        r2_large = r2_list[r2_list >= self.lcut**2]
        
        n2_list = np.sum(((self.r_vec[:,:] + 0.5*np.array(self.d))*gv)**2, axis=1)
        n2_ne0   = n2_list[n2_list >= 1.0**2] 
        
        integ_I_low = lambda t, q2: 0.5*np.pi * gamma * (np.exp(t*q2) - 1.0)/(t**1.5) \
                + 0.5*np.pi * gamma * np.sum(np.exp(t*q2 - np.pi**2 * n2_ne0/t))/(t**1.5) \
                - (1.0/(np.sqrt(4.0*np.pi))) * np.sum(np.exp((q2 - r2_small) * t))

        integ_I_high = lambda t, q2: (1.0/(np.sqrt(4.0*np.pi))) \
                * np.sum(np.exp((q2 - r2_large) * t)) - 0.5 * gamma * np.pi / (t**1.5)

                
        I_low, I_low_err = scipy.integrate.quad(integ_I_low,tmin,tcut,args=(q2))
        I_high, I_high_err = scipy.integrate.quad(integ_I_high,tcut,tmax, args=(q2))

        z_0 = (1.0/np.sqrt(4.0*np.pi)) * np.sum(1.0/(r2_small - q2))
        
        return z_0 + I_low + I_high

