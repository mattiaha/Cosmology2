# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 13:35:56 2022

@author: matti
"""
import matplotlib.pyplot as plt
import numpy as np

c = 2.999E+8
MPc = 3.08567758e22
GYear = 365 * 24 * 3600 * 10**9
H0 = 2.17132e-18 #km/s/Mpc

x = np.loadtxt("cosmology.txt")[:,0]
eta = np.loadtxt("cosmology.txt")[:,1]
Hp = np.loadtxt("cosmology.txt")[:,2]
dHpdx = np.loadtxt("cosmology.txt")[:,3]
OmegaB = np.loadtxt("cosmology.txt")[:,4]
OmegaCDM = np.loadtxt("cosmology.txt")[:,5]
OmegaLambda = np.loadtxt("cosmology.txt")[:,6]
OmegaR = np.loadtxt("cosmology.txt")[:,7]
H = np.loadtxt("cosmology.txt")[:,8]
ddHpddx = np.loadtxt("cosmology.txt")[:,9]
time =  np.loadtxt("cosmology.txt")[:,10]
z = np.loadtxt("cosmology.txt")[:,11]
dL = np.loadtxt("cosmology.txt")[:,12]


plt.figure(0)
plt.plot(x,H, 'r', label = 'H(x)')
plt.yscale('log')
plt.xlabel('x')
plt.ylabel('H(x)/$H_{0}$ ')
plt.title('Evolution of Hubble parameter')
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')

plt.show()

plt.figure(1)
plt.plot(x,1/Hp*dHpdx, 'r', label = '1/$\mathcal{H}$((x) * d$\mathcal{H}$(x)/dx')
plt.plot(x,(1/Hp)*ddHpddx, 'm' , label = '1/$\mathcal{H}$(x) * dd$\mathcal{H}$(x)/ddx')
plt.xlabel('x')
plt.title('Evolution of derivatives of $\mathcal{H}$(x)')
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')


plt.legend()
plt.show()

plt.figure(2)
plt.plot(x,eta, 'm')
plt.title('Conformal time')
plt.xlabel('x')
plt.ylabel('$\eta$ (x) [Mpc]')
plt.yscale('log')
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')

plt.show()

plt.figure(3)
plt.plot(x,eta*Hp*H0*MPc/c, 'g')
plt.title('$\eta$(x) $\mathcal{H}$ (x) / c')
plt.xlabel(('x'))
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')

plt.show()

plt.figure(4)
plt.plot(x,time, 'r')
plt.ylabel('t(x) [GY]')
plt.xlabel('x')
plt.title('Age of the universe')
plt.yscale('log')

plt.axvline(-8.657, color = 'c', linestyle='--')

plt.axvline(-0.487, color = 'm' , linestyle='--')
plt.show()

plt.figure(5)
plt.title('Density functions')
plt.ylabel('$\Omega$ (x)')
plt.xlabel('x')
plt.plot(x,OmegaB+OmegaCDM, 'r',label = '$\Omega _b + \Omega _{CDM}$' )
plt.plot(x,OmegaR, 'k', label = '$\Omega _{\gamma}$')
plt.plot(x,OmegaLambda, 'g', label = '$\Omega _{\Lambda}$')
plt.legend()
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')
plt.axvline(-0.49, color = 'b' , linestyle='--')


plt.show()


z2 = np.loadtxt("supernova.txt")[:,0]
dL2 = np.loadtxt("supernova.txt")[:,1]
err1 = np.loadtxt("supernova.txt")[:,2]

plt.figure(6)
plt.title('Luminosity distance compared to supernova data')
plt.errorbar(z2,dL2, yerr= err1,color = 'b', fmt = 'o', label = 'Supernova data')
plt.xlabel('z')
plt.ylabel('$Distance$ [Gpc]')
plt.axis([1E-2,2.0,1E-2,10])

plt.plot(z,dL,color = 'r', label = 'Theoretical luminosity distance')
plt.show()
plt.figure(7)

plt.plot(x,Hp, 'b', label = 'Hp(x)')
plt.xlabel('x')
plt.ylabel('$\mathcal{H}$(x)/$H_{0}$ ')
plt.title('Evolution of scaled Hubble parameter')

plt.yscale('log')
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')

plt.show()
