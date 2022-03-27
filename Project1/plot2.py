# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 21:53:24 2022

@author: matti
"""
import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt("recombination.txt")[:,0]
Xe = np.loadtxt("recombination.txt")[:,1]
ne = np.loadtxt("recombination.txt")[:,2]
tau = np.loadtxt("recombination.txt")[:,3]
dtaudx = np.loadtxt("recombination.txt")[:,4]
ddtauddx = np.loadtxt("recombination.txt")[:,5]
g_tilde = np.loadtxt("recombination.txt")[:,6]
dgdx_tilde = np.loadtxt("recombination.txt")[:,7]
ddgddx_tilde = np.loadtxt("recombination.txt")[:,8]
Xe_saha = np.loadtxt("recombination.txt")[:,9]

plt.figure(0)
plt.plot(x,Xe, 'r', label = '$X_e$')
plt.plot(x, Xe_saha, 'b', linestyle ='dashed' , label = 'Saha approximation' )
plt.yscale('log')
plt.ylabel("$X_e$")
plt.legend()
plt.title('Fractional electron density')
plt.xlabel('x')
plt.axvline(-6.969, color = 'g' , linestyle='--')

plt.show()




plt.figure(1)
plt.plot(x,ne)
plt.title("ne")
plt.yscale('log')
plt.show()

plt.figure(2)
plt.plot(x,tau, 'r', label = '$ \u03C4 $' )
plt.plot(x, -dtaudx, 'g', label = '- d\u03C4 /dx')
plt.plot(x,ddtauddx, 'k', label = 'd$^2$\u03C4 /dx$^2$')
plt.xlabel('x')
plt.title('Optical depth')
plt.legend()
plt.axvline(-6.969, color = 'b' , linestyle='--')

plt.yscale('log')
plt.show()

plt.figure(3)
plt.plot(x,g_tilde, label = 'g')
plt.title('Visibility function')
plt.ylabel('g (x)')
plt.xlabel('x')
plt.show()

plt.figure(4)
plt.plot(x, dgdx_tilde, label = ' dg /dx')
plt.ylabel('dg/dx')
plt.xlabel('x')
plt.show()

plt.figure(5)
plt.plot(x,ddgddx_tilde, label = 'd$^2g/dx^2$')
plt.xlabel('x')
plt.ylabel('d$^2 $g/dx$^2$')
plt.show()


