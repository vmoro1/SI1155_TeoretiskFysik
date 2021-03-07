#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:41:31 2020

@author: viggomoro
"""

import matplotlib.pyplot as plt

# E = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
# T = [0.31, 0.25, 0.21, 0.15, 0.12, 0.082, 0.048, 0.054, 0.040, 0.0091]
E = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
T = [0.25, 0.21, 0.15, 0.12, 0.082, 0.048, 0.054, 0.040, 0.0091]
for i in range(len(T)):
    x = T[i] / (16*E[i]*(1-E[i]))
    T[i] = x
  
kw = []    
for i in range(9):
    x = (2*(1-E[i]))**(0.5)
    kw.append(4*x)
       
   

plt.semilogy(kw, T)
plt.ylabel('T/[16E(V_0 − E)/V_0^2] ')
plt.xlabel('2kw')