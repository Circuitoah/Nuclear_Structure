import numpy as np
from matplotlib import pyplot as plt

#constants
e = np.sqrt(1.4399764) #Mev fm
a  = 0.7 #fm 
hbar_mpi_c_squared = 2.044 # fm**2

hbar = 197.32 # 

# Potentals
V_0 = 54 # MeV
V_SO = 0.2 *V_0

mu  = 938

R = 1.2 * A**(1/3) # fm 
#------------------------------------------------------------------------------------------------------------------------

A = [109, 112, 147, 147, 150, 151, 156, 160 ]
Z = [53, 55, 69, 69, 71, 71, 73, 75]
Q = [.829, .977, 1.071, 1.139, 1.285, 1.255 , 1.015, 1250]
n = [1, 0 ,0 , 1 , 0, 0, 1, 1 ]
l = [2, 4 , 5, 2, 5, 5, 2, 2 ]
j = [5/2, 7/2, 11/2 , 3/2 , 11/2, 11/2 , 3/2, 3/2]
Gamma_exp = [1.09e-4 , 3.3e-5, 2.7, 3.6e-4, 1e-2, 1.2e-1, 1.65e-1, 8.7e-4]

#Constants etc----------------------------------------------------------------------------------------------------------


e = np.sqrt(1.4399764) #Mev fm
a  = 0.7 #fm 
hbar_mpi_c_squared = 2.044 # fm**2

hbar = 197.32 # 

# Potentals
V_0 = 54 # MeV
V_SO = 0.2 *V_0

mu  = 938 #Mev

R = 1.2 * A[0]**(1/3) # fm 



#Function----------------------------------------------------------------------------------------------------------------
def V_coul(r):

    if (r > R):
        result = (Z[0] * e**2)/ r

    else:
        result = (Z[0] * e**2) * (3  - (r/R)**2 )/(2 * R)

    return result

fws = lambda r : (1 + np.exp( (r - R)/a ))

dfws_dr = lambda r : - (  np.exp((r- R) / a )/a ) * ( 1 + np.exp((r - R)/a) )**2

l_dot_s  = lambda j, l : 1/2 *( j* (j + 1) - l * (l  + 1) - .5 * (1.5) )

V_ws = lambda r, j , l : - V_0 * fws(r) + V_SO * hbar_mpi_c_squared * (1/r) * dfws_dr(r) * l_dot_s(j , l)

V =lambda r , j , l :  V_ws(r,j, l) + V_coul(r) +  (hbar**2)/ (2 * r * mu**2)* (l * (l + 1 ))
#--------------------------------------------------------------------------------------------------------------------------

r_test = np.linspace(0.1, 100, 1000)

V_test = [V(i, j[0], l[0]) for i in r_test.tolist()]

plt.plot(r_test, V_test)
plt.show()