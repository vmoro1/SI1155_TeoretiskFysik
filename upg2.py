#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 08:21:51 2020

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
N=100000             # number of mesh points
dx=5*a/N             # step length
dx2=dx**2          # step length squared

EeV = -4          # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5
E = EeV*e          # input energy in J


# potential energy function
def V(x):
    y = 0.0
    #y=e*5.*x/a # triangular potential
    if x>0. and x<1. : y=e*5.*(x/a-1.) # finite triangular potential
    return y

# initial values and lists
x = -a             # initial value of position x
psi = 0.0           # wave function at initial position
dpsi = 1.0          # derivative of wave function at initial position
x_tab = []          # list to store positions for plot
psi_tab = []        # list to store wave function for plot
x_tab.append(x/a)
psi_tab.append(psi)

# for i in range(N) :
#     d2psi = c*(V(x)-E)*psi
#     d2psinew = c*(V(x+dx)-E)*psi
#     psi += dpsi*dx + 0.5*d2psi*dx2
#     dpsi += 0.5*(d2psi+d2psinew)*dx
#     x += dx
#     x_tab.append(x/a)
#     psi_tab.append(psi)


def sol(EeV):
    E = e*EeV
    x = -3*a
    psi = 0.0           
    dpsi = 1.0          
    for i in range(100000) :
        d2psi = c*(V(x)-E)*psi
        d2psinew = c*(V(x+dx)-E)*psi
        psi += dpsi*dx + 0.5*d2psi*dx2
        dpsi += 0.5*(d2psi+d2psinew)*dx
        x += dx
        x_tab.append(x/a)
        psi_tab.append(psi)
    return psi


def hitta_energier(E_left, E_right):
    for i in range(50):
        E = (E_right + E_left)/2
        psi = sol(E)
        if psi >= 0:
            E_left = E
        else:
            E_right = E
    return E

print(hitta_energier(-1, -0.1))

# print('E=',EeV,'eV , psi(x=a)=',psi)

# plt.close()
# plt.figure(num=None, figsize=(8,8), dpi=80, facecolor='w', edgecolor='k')
# plt.plot(x_tab, psi_tab, linewidth=1, color='red')
# # plt.xlim(0, 1)
# #limit=1.e-9
# #plt.ylim(0, limit)
# #plt.ylim(-limit, limit)
# #plt.autoscale(False)
# plt.xlabel('x/a')
# plt.ylabel('$\psi$')
# #plt.savefig('psi.pdf')
# plt.show()
