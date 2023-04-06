# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 00:01:52 2023

@author: nbarbour13
"""

import numpy as np
import matplotlib.pyplot as plt

# rk4 derivative function
def frk3(g,ik,ids,mcut,chit):
    M = g.shape[0]
    gtmp = np.zeros(M, dtype=complex)
    gtmp[0] = -ik*g[1]/np.sqrt(2)  
    gtmp[1] = -ik*(g[2] + np.sqrt(2)*g[0])
    gtmp[2] = -ik*(np.sqrt(3/2)*g[3] + g[1])  + chit  
    gtmp[3:M-1] = -ik*(np.sqrt(ids[4:M]/2)*g[4:M] + np.sqrt(ids[3:M-1]/2)*g[2:M-2]) - g[3:M-1]*(ids[3:M-1]/mcut)**6
    gtmp[M-1] = -ik*(np.sqrt((M-1)/2)*g[M-2]) - g[M-1]*((M-1)/mcut)**6
    return gtmp

# standard rk4 timestepping function
def RK43(g,ik,ids,M,chit,dt):
    k1 = frk3(g          ,ik,ids,M,chit)
    k2 = frk3(g + dt*k1/2,ik,ids,M,chit)
    k3 = frk3(g + dt*k2/2,ik,ids,M,chit)
    k4 = frk3(g + dt*k3  ,ik,ids,M,chit)
    return g + dt*(k1 + 2*k2 + 2*k3 + k4)/6

# Generate a dataset. returns G and chi(t) for all t
def run_sim(simparams,seedint,chis=None):
    G = np.zeros((simparams['M'],simparams['nt']+1),dtype=complex)
    np.random.seed(seedint)
    if chis is None:
        chis = np.random.rand(simparams['nt'])
    for j in np.arange(simparams['nt']):
        G[:,j+1] = RK43(G[:,j], simparams['ik'], simparams['m'], simparams['simcut'], chis[j], simparams['dt'])
    return G, chis

def absq(g):
    return np.abs(g)**2

#%%
simparams = {
    "M": 100,  # Total number of moments for the high resolution simulation
    "simcut": 64, # hyperviscosity -> order 1 at simmcut, (m/simcut)**6 
    "L": 7*np.pi, # box size in spatial z
    "nt": 200000, 
    "t": 0, # initial time
    "dt": 0.01
    }
simparams["m"] = np.arange(simparams["M"])  # labels Hermite index
simparams["ik"] = complex(0,1)*2*np.pi/simparams["L"] # precalculate Fourier derivative 

seedint = 737
G, chis = run_sim(simparams, seedint)

# Plot the dataset
g2 = absq(G)
g2 = np.transpose(g2)  
#g2 = np.flip(g2,axis=0)
t = 10000
plt.imshow(g2[:t,:],cmap='jet',aspect=1,extent=[0,simparams['M'],t*simparams['dt'],0])
plt.colorbar()
plt.xlabel('m')
plt.ylabel('t')
plt.title('$|g|^2$, L=%d$\pi$, dt=%.02f' % (int(simparams['L']/np.pi), simparams['dt']))
plt.show()