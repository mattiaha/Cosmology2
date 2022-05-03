# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 16:05:25 2022

@author: matti
"""

import matplotlib.pyplot as plt
import numpy as np

d_cdm1 = np.loadtxt("perturbations_k0.001.txt")[:,1]
d_b1 = np.loadtxt("perturbations_k0.001.txt")[:,2]
v_cdm1 = np.loadtxt("perturbations_k0.001.txt")[:,3]
v_b1 = np.loadtxt("perturbations_k0.001.txt")[:,4]
theta0_1 = np.loadtxt("perturbations_k0.001.txt")[:,5]
theta1_1 = np.loadtxt("perturbations_k0.001.txt")[:,6]
Phi1 = np.loadtxt("perturbations_k0.001.txt")[:,7]
Psi1 = np.loadtxt("perturbations_k0.001.txt")[:,8]


x = np.loadtxt("perturbations_k0.01.txt")[:,0]
d_cdm2 = np.loadtxt("perturbations_k0.01.txt")[:,1]
d_b2 = np.loadtxt("perturbations_k0.01.txt")[:,2]
v_cdm2 = np.loadtxt("perturbations_k0.01.txt")[:,3]
v_b2 = np.loadtxt("perturbations_k0.01.txt")[:,4]
theta0_2 = np.loadtxt("perturbations_k0.01.txt")[:,5]
theta1_2 = np.loadtxt("perturbations_k0.01.txt")[:,6]
Phi2 = np.loadtxt("perturbations_k0.01.txt")[:,7]
Psi2 = np.loadtxt("perturbations_k0.01.txt")[:,8]

d_cdm3 = np.loadtxt("perturbations_k0.1.txt")[:,1]
d_b3 = np.loadtxt("perturbations_k0.1.txt")[:,2]
v_cdm3 = np.loadtxt("perturbations_k0.1.txt")[:,3]
v_b3 = np.loadtxt("perturbations_k0.1.txt")[:,4]
theta0_3 = np.loadtxt("perturbations_k0.1.txt")[:,5]
theta1_3 = np.loadtxt("perturbations_k0.1.txt")[:,6]
Phi3 = np.loadtxt("perturbations_k0.1.txt")[:,7]
Psi3 = np.loadtxt("perturbations_k0.1.txt")[:,8]


plt.figure(0)
plt.plot(x,d_cdm1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x,d_cdm2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x,d_cdm3, 'g', label = 'k = 0.1/Mpc')
plt.plot(x,d_b1, 'b',linestyle='dashed')
plt.plot(x,d_b2, 'r', linestyle='dashed')
plt.plot(x,d_b3, 'g', linestyle='dashed')
plt.yscale('log')
plt.legend()
plt.title('$\u03B4_{CMD}(x), \u03B4_{b}(x)$')
plt.xlabel('x')
plt.show()

plt.figure(1)
plt.plot(x,v_cdm1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x,v_cdm2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x,v_cdm3, 'g', label = 'k = 0.1/Mpc')
plt.plot(x,v_b1, 'b',linestyle='dashed')
plt.plot(x,v_b2, 'r', linestyle='dashed')
plt.plot(x,v_b3, 'g', linestyle='dashed')
plt.yscale('log')
plt.title('$v_{cdm}(x) , v_{b}(x)$')
plt.xlabel('x')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(x,4*theta0_1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x,4*theta0_2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x,4*theta0_3, 'g', label = 'k = 0.1/Mpc')
plt.legend()
plt.xlabel('x')
plt.title('$\delta_{\gamma}(x) = 4\Theta_0$(x)')
plt.show()


plt.figure(3)
plt.plot(x,-3*theta1_1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x,-3*theta1_2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x,-3*theta1_3, 'g', label = 'k = 0.1/Mpc')
plt.legend()
plt.title('$v_{\gamma}(x) = -3\Theta_1$(x)')
plt.xlabel('x')
plt.show()

plt.figure(4)
plt.plot(x,Phi1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x,Phi2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x,Phi3, 'g', label = 'k = 0.1/Mpc')
plt.xlabel('x')
plt.legend()
plt.title('Gravitational potential $\Phi$')
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')
plt.ylabel('$\Phi$(x)')
plt.show()

plt.figure(5)
plt.plot(x,Psi1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x,Psi2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x,Psi3, 'g', label = 'k = 0.1/Mpc')
plt.title('Gravitational potential $\Psi$')
plt.axvline(-8.657, color = 'c', linestyle='--')
plt.axvline(-0.2558, color = 'b' , linestyle='--')
plt.ylabel('$\Psi$(x)')
plt.legend()
plt.xlabel('x')
plt.show()

plt.figure(6)
plt.plot(x, Psi1+Phi1, 'b', label = 'k = 0.001/Mpc')
plt.plot(x, Psi2+Phi2, 'r', label = 'k = 0.01/Mpc')
plt.plot(x, Psi3+Phi3, 'g', label = 'k = 0.1/Mpc')
plt.ylabel('$\Psi (x) + \Phi (x)$')
plt.legend()
plt.title('Sum of potentials')
plt.xlabel('x')
plt.show()

