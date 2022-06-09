# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:33:45 2022

@author: matti
"""

import matplotlib.pyplot as plt
import numpy as np

c = 2.999E+8
MPc = 3.08567758e22
H0 = 2.17132e-18 #km/s/Mpc

ell = np.loadtxt("cells.txt")[:,0]
C_ell = np.loadtxt("cells.txt")[:,1]
"""SW = np.loadtxt("cells.txt")[:,2]
ISW = np.loadtxt("cells.txt")[:,3]
Doppler = np.loadtxt("cells.txt")[:,4]
Term4 = np.loadtxt("cells.txt")[:,5]"""

k = np.loadtxt("MPS.txt")[:,0]
MPS = np.loadtxt("MPS.txt")[:,1]
Theta6 = np.loadtxt("MPS.txt")[:,2]
Theta100 = np.loadtxt("MPS.txt")[:,3]
Theta200 = np.loadtxt("MPS.txt")[:,4]
Theta500 = np.loadtxt("MPS.txt")[:,5]
Theta1000 = np.loadtxt("MPS.txt")[:,6]

err_ell = np.loadtxt("planck_cell.txt")[:,0]
data = np.loadtxt("planck_cell.txt")[:,1]
data_p = np.loadtxt("planck_cell.txt")[:,2]
data_n = np.loadtxt("planck_cell.txt")[:,3]
error = np.array([data_p,data_n])

err_k = np.loadtxt("data_mps.txt")[:,0]
err_P = np.loadtxt("data_mps.txt")[:,1]
error2 = np.loadtxt("data_mps.txt")[:,2]
more_error_k = np.loadtxt("more_MPS.txt")[:,0]
more_error_p = np.loadtxt("more_MPS.txt")[:,1]
more_error = np.loadtxt("more_MPS.txt")[:,2]
k_eq = 4.372E-25*MPc/0.67
more_err_arr = np.array([np.zeros(len(more_error)),more_error])

print(more_err_arr)

plt.figure(0)
plt.plot(ell,C_ell)
plt.errorbar(err_ell,data, yerr=error,fmt='x')
#plt.yscale('log')
plt.xscale('log')
plt.xlabel('Multipole l')
plt.ylabel(' Temperature fluctuations [$\mu k^2$]')
plt.title('CMB Power Spectrum')
plt.show()


plt.figure(1)
plt.plot(k,MPS)
plt.errorbar(err_k,err_P,yerr=error2,fmt='x')
plt.errorbar(more_error_k,more_error_p,yerr=more_err_arr,fmt='x')

plt.yscale('log')
plt.xscale('log')
plt.title('Matter Power-Spectrum')
plt.xlabel('Wavenumber k [h/Mpc]')
plt.ylabel('P(k) [Mpc/h]$^{3}$')
plt.axvline(k_eq, color = 'c', linestyle='--')
plt.show()

k2 = k*c*0.7/(MPc*H0)
plt.figure(2)
plt.plot(k2,Theta6, label= 'l = 6')
plt.plot(k2,Theta100, label = 'l = 100')
plt.plot(k2,Theta200 , label = 'l = 200')
plt.plot(k2,Theta500, label = 'l = 500')
plt.plot(k2,Theta1000, label = 'l = 1000')
plt.ylabel('$\Theta_l(k)$')
plt.xlabel('ck/$H_0$')
plt.title('Transfer function')
plt.legend()
plt.show()

plt.figure(3)
plt.plot(k2,Theta6*Theta6/k2, label= 'l = 6')
plt.plot(k2,Theta100*Theta100/k2, label = 'l = 100')
plt.plot(k2,Theta200*Theta200/k2 , label = 'l = 200')
plt.plot(k2,Theta500*Theta500/k2, label = 'l = 500')
plt.plot(k2,Theta1000*Theta1000/k2, label = 'l = 1000')
plt.ylabel('$\Theta_l(k)^2 /k $')
plt.xlabel('ck/$H_0$')
plt.title('Spectrum integrand')
plt.legend()
plt.show()
"""
plt.figure(4)
plt.plot(ell,SW, 'r', label = 'SW')
plt.plot(ell,ISW, 'b', label = 'ISW')
plt.plot(ell,Doppler,'k', label = 'Doppler')
plt.plot(ell,Term4,'c', label = '4th term')
plt.yscale('log')
plt.xscale('log')
plt.title('Components of source function')
plt.legend()
plt.show()"""
