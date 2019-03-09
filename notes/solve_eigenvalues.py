# solving eigen values using HAL QCD pot.
import pickle

import os
import scipy
import numpy as np

import pickle
import scipy.constants as sc
hbarc = ( sc.hbar * sc.speed_of_light
         / sc.mega / sc.electron_volt / sc.femto )
ainv = 2.194e3 # lattice cutoff in GeV
lat_unit = hbarc/ainv # lattice spacing in fm

mxi = 0.665 # mass of Xi in lattice unit
mpi = 510/ainv # mass of pion in lattice unit
m_red = 0.5 * mxi # reduced mass

from scipy.sparse.linalg import LinearOperator


if os.path.isfile('pkls/v0_4gauss_func_prm_av_jk.pkl'):
    with open('pkls/v0_4gauss_func_prm_av_jk.pkl', 'rb') as fin:
        prm_v0_4gauss_avs, prm_v0_4gauss_jks = pickle.load(fin)
else:
    print('Run "Fitting and Scattering phase shift.ipynb"')

if os.path.isfile('pkls/v2_2gauss_func_prm_av_jk.pkl'):
    with open('pkls/v2_2gauss_func_prm_av_jk.pkl', 'rb') as fin:
        prm_v2_2gauss_avs, prm_v2_2gauss_jks = pickle.load(fin)
else:
    print('Run "Fitting and Scattering phase shift.ipynb"')


V4gauss = lambda p, x: (p[0]*np.exp(-p[1]*x**2) + p[2]*np.exp(-p[3]*x**2)
                        + p[4]*np.exp(-p[5]*x**2) + p[6]*np.exp(-p[7]*x**2))
                        
V2gauss = lambda p, x: p[0]*np.exp(-p[1]*(x-p[2])**2) + p[3]*np.exp(-p[4]*(x-p[5])**2)

it0 = 13
bin_num =  20


def A1_projection(wave_in):
    Ns = round(len(wave_in)**(1/3))
    wave = wave_in.reshape(Ns,Ns,Ns)
    wave_tmp1 = (wave[:,:,:] + np.roll(wave,-1,0)[::-1,:,:]
                + np.roll(wave,-1,1)[:,::-1,:]
                + np.roll(wave,-1,2)[:,:,::-1]
                + np.roll(np.roll(wave,-1,0),-1,1)[::-1,::-1,:]
                + np.roll(np.roll(wave,-1,1),-1,2)[:,::-1,::-1]
                + np.roll(np.roll(wave,-1,2),-1,0)[::-1,:,::-1]
                + np.roll(np.roll(np.roll(wave,-1,0),-1,1),-1,2)[::-1,::-1,::-1])/8.0
    wave_tmp2 = (wave_tmp1 
                + np.swapaxes(wave_tmp1,0,1)
                + np.swapaxes(wave_tmp1,1,2)
                + np.swapaxes(wave_tmp1,2,0)
                + np.swapaxes(np.swapaxes(wave_tmp1,0,1),1,2)
                + np.swapaxes(np.swapaxes(wave_tmp1,0,2),2,1))/6.0e0

    return wave_tmp2.flatten()



def solve_eigen(ibin, L=40, Nev=5, it0=13):
    print(ibin)
    rs = np.array([np.sqrt(ix**2+iy**2+iz**2) + 1.0e-5
                  for iz in range(-L//2,L//2) for iy in range(-L//2,L//2)
                  for ix in range(-L//2,L//2)]).reshape(L,L,L)
    rs = np.roll(rs, (L//2,L//2,L//2), (0, 1, 2)).flatten()
    
    pot_v0 = V4gauss(prm_v0_4gauss_jks[('n2lo',it0)][ibin,:], rs)
    pot_v2 = V2gauss(prm_v2_2gauss_jks[it0][ibin,:], rs)

    lap = lambda vec: - 6.0*vec + ( np.roll(vec,+1,0) + np.roll(vec,-1,0) 
                                  + np.roll(vec,+1,1) + np.roll(vec,-1,1) 
                                  + np.roll(vec,+1,2) + np.roll(vec,-1,2) )

    Vol = L**3

    H = LinearOperator((Vol,Vol), 
        matvec = lambda vec: (pot_v0 * vec 
          + lap(vec.reshape(L,L,L)).flatten()*(pot_v2 - 1/(2.0*m_red))), 
        dtype='float64')

    vals, vecs = scipy.sparse.linalg.eigs(H, which='SM', k = Nev)
    return vals, vecs

vals_jk = {}
vecs_jk = {}

for L in [40, 48, 64]:
    print(L)
    for ibin in range(bin_num):
        print(ibin)
        _vals, _vecs = solve_eigen(ibin, L=L)
        vals_jk[L].append(_vals)
        vecs_jk[L].append(_vecs)

    vals_jk[L] = np.array(vals_jk[L])
    vecs_jk[L] = np.array(vecs_jk[L])

for L in [40, 48, 64]:
    vals_jk[L] = np.array(vals_jk[L])
    vecs_jk[L] = np.array(vecs_jk[L])

eigens_jk = {}
for L in [40, 48, 64]:
    eigens_jk[L] = vals_jk[L][:,0]

eigen_scat = vals_jk[64][:,1]

with open('pkls/finite_volume_spectra_n2lo.pkl', 'wb') as fout:
    pickle.dump((eigens_jk, eigen_scat), fout)
