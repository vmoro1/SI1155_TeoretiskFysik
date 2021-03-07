#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 13:42:05 2020

@author: viggomoro
"""

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

a=1.0e-9           # well width a=1 nm
hbar=1.0546e-34    # Plancks constant
m=9.1094e-31       # electron mass
e=1.6022e-19       # electron charge=-e
c=2.0*m/hbar**2    # constant in Schrödinger equation
N=1000             # number of mesh points
dx=a/N             # step length
dx2=dx**2          # step length squared

EeV = 4 * (hbar*np.pi/a)**2/(2.0*m)/e        # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5 . Energiegenvärdena kommer gen av E1*n^2
E = EeV*e          # input energy in J

def V(x):
    y = 0.0
    #y=e*5.*x/a # triangular potential
    #if x>0. and x<1. : y=e*5.*(x/a-1.) # finite triangular potential
    return y
    
def sol(N, EeV=(hbar*np.pi/a)**2/(2.0*m)/e):
    x = 0.0            
    psi = 0.0
    dpsi = 1.0 
    dx=a/N             
    dx2=dx**2
    E = EeV * e
    
    for i in range(N):
        d2psi = c*(V(x)-E)*psi
        d2psinew = c*(V(x+dx)-E)*psi
        psi += dpsi*dx + 0.5*d2psi*dx2
        dpsi += 0.5*(d2psi+d2psinew)*dx
        x += dx
    return abs(psi)


N = [1000, 10000, 100000, 1000000]
psi_vals = []
for n in N:
    psi_vals.append(sol(n))
    
x = []
y = []    
for i in range(10**3,(10**6)+1):    
    x.append(i)
    y.append((10**-9)/(i**2))

plt.loglog(N, psi_vals, marker='*', linestyle='-', label='Integrationsfelet')
plt.loglog(x,y, label='(10^-9) / (N^2)') 
plt.legend(loc='best')
plt.xlabel('Antalet indelningspunkter') 
plt.ylabel('Integratoinsfelet (abs(psi(a)))')


def hitta_energier(E_left, E_right):
    while True:
        E = (E_right + E_left)/2
        psi = sol(1000, E)
        if abs(psi) < 1e-11:
            return E
        elif psi > 0:
            E_left = E
        else:
            E_right = E
            
        
# print(hitta_energier(3.38, 3.39))      






x = 0            # initial value of position x
psi = 0.0           # wave function at initial position
dpsi = 1.0          # derivative of wave function at initial position
x_tab = []          # list to store positions for plot
psi_tab = []        # list to store wave function for plot
x_tab.append(x/a)
psi_tab.append(psi)

N = 1000
for i in range(N):
    d2psi = c*(V(x)-E)*psi
    d2psinew = c*(V(x+dx)-E)*psi
    psi += dpsi*dx + 0.5*d2psi*dx2
    dpsi += 0.5*(d2psi+d2psinew)*dx
    x += dx
    x_tab.append(x/a)
    psi_tab.append(psi)
  
# print(max(psi_tab))     
# print(psi, dpsi)

# plt.close()
# plt.figure(num=None, figsize=(8,8), dpi=80, facecolor='w', edgecolor='k')
# plt.plot(x_tab, psi_tab, linewidth=1, color='red')
# plt.xlim(0, 1)


# limit=1.e-9
# plt.ylim(0, limit)
# plt.ylim(-limit, limit)
# plt.autoscale(False)
# plt.xlabel('x/a')
# plt.ylabel('$\psi$ för E2, det första exciterade tillståndet')
# # #plt.savefig('psi.pdf')
# plt.show()

           
           



           
           