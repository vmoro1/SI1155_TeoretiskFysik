
""" 
Matrix diagonalization
"""

import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt

a = np.arange(0.0, 1.01, 0.01)
E1=[]
E2=[]
E3=[]
E4=[]

for i in range(101) :
    x=a[i] # x är a
    y = 1 - x # y är b
    # print(x)
    
    H = np.array([[2*x + y, 0, 0, 0], [0, -y, 2*y, 0], [0, 2*y, -y, 0], [0, 0, 0, -2*x + y]])
    
    E=LA.eigvalsh(H) # eigvalsh beräknar egenvärden till hermiteska matriser
    
    E1.append(E[0])
    E2.append(E[1])
    E3.append(E[2])
    E4.append(E[3])

    # print(x,E[0])

# print(E)


plt.figure(num=None, figsize=(12,8), dpi=80, facecolor='w', edgecolor='k')
plt.plot(a, E1, 'r-', label='E_1')
plt.plot(a, E2, 'g-', label='E_2')
plt.plot(a, E3, 'b-', label='E_3')
plt.plot(a, E4, 'k-', label='E_4')
plt.legend(loc="upper left")
# plt.plot(a, E1, 'ro-', a, E2, 'go-', a, E3, 'bo-', a, E4, 'ko-')
plt.legend(loc='best')
plt.xlabel('$a$')
plt.ylabel('$E_n$')
# plt.savefig('eigenvalues.pdf')
plt.show()










