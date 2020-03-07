# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 20:58:56 2016

@author: Admin
"""
import timeit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import linalg as lin
from scipy.optimize import fsolve
from scipy.linalg import solve
from scipy.linalg import solve_banded
from scipy.special import jn_zeros #prvi parameter je order, drugi št. ničel
from scipy.special import jv #prvi order drugi argument
from scipy.special import beta
import scipy.sparse
from scipy.optimize import root
from scipy.integrate import quad
from scipy.integrate import complex_ode
from scipy.optimize import linprog
import scipy.optimize as opt

# Plot the surface.
def volterra(t,x,p): #zajci lisice
    return [p*x[0]*(1-x[1]),x[1]/p*(x[0]-1)]

def laser(t,x,p,q): #atomi fotoni
    return [q-p*x[0]*(x[1]+1),x[1]/p*(x[0]-1)]
def epidemija(t,x,p): #D B I
    return [-x[0]*x[1],x[0]*x[1]-p*x[1],p*x[1]]
def integriraj(zacetni,p,cas,tol=0.1):
    r = ode(volterra).set_integrator("dopri5")

    r.set_initial_value(zacetni,0).set_f_params(p)
    rezultati = []
    perioda = 0
    per = True
    for t in np.linspace(0.1,cas,1000):
        trenutno = r.integrate(t)
        if abs(trenutno[0]-zacetni[0])<tol and t>2:
            if abs(trenutno[1]-zacetni[1])<tol and per:
                perioda = t
                per = False
        rezultati.append(trenutno)
    return [np.asarray(rezultati),perioda]
barve = cm.get_cmap("rainbow")

def integriraj2(zacetni,p,q,cas,tol=0.1):
    r = ode(laser).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p,q)
    rezultati = []
    #perioda = 0
    #per = True
    for t in np.linspace(0.1,cas,10000):
        trenutno = r.integrate(t)
        #if abs(trenutno[0]-zacetni[0])<tol and t>2:
            #if abs(trenutno[1]-zacetni[1])<tol and per:
                #perioda = t
                #per = False
        rezultati.append(trenutno)
    return np.asarray(rezultati)

def hann(n, N):
    return np.sin(np.pi*n/(N-1))**2

def integriraj3(zacetni,p,cas,tol=0.1):
    r = ode(epidemija).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p)
    rezultati = []
    #perioda = 0
    #per = True
    for t in np.linspace(0.001,cas,10000):
        trenutno = r.integrate(t)
        #if abs(trenutno[0]-zacetni[0])<tol and t>2:
            #if abs(trenutno[1]-zacetni[1])<tol and per:
                #perioda = t
                #per = False
        rezultati.append(trenutno)
    return np.asarray(rezultati)
def epidemija2(t,x,alfa,beta1,beta2):
    return [-alfa*x[0]*x[1],alfa*x[0]*x[1]-beta1*x[1],beta1*x[1] - beta2*x[2],beta2*x[2]]
def integriraj4(zacetni,alfa,beta1,beta2,cas,tol=0.1):
    r = ode(epidemija2).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(alfa,beta1,beta2)
    rezultati = []
    #perioda = 0
    #per = True
    for t in np.linspace(0.001,cas,10000):
        trenutno = r.integrate(t)
        #if abs(trenutno[0]-zacetni[0])<tol and t>2:
            #if abs(trenutno[1]-zacetni[1])<tol and per:
                #perioda = t
                #per = False
        rezultati.append(trenutno)
    return np.asarray(rezultati)

tajm = 10
cas = np.linspace(0.001,tajm,10000)



rezultati=integriraj4((10,10,0,0),3,2,1,tajm)
plt.plot(cas,rezultati[:,0],label="Dovzetni")
plt.plot(cas,rezultati[:,1]+rezultati[:,2],label="Bolani=B1+B2")
plt.plot(cas,rezultati[:,3],label="Imuni")
plt.legend(loc="best")
plt.title(r"$\alpha=3, \beta_1=2, \beta_3=1\ (D_0,B_{10},B_{20},I_0) = (10,10,0,0)$")
plt.savefig("epidemija/dodatna6.pdf")
plt.show()


































     
"""    
N = len(rezultati[:,1])
fourier = [hann(n,N)*rezultati[:,0][n] for n in range(len(rezultati[:,0]))
"""
"""
pp =[0,0.05,0.01] #qji 

for p in pp:
    casi = np.linspace(0.1,10,10000)
    rezultati = integriraj2((1,2/0.03-1+p),0.03,2,10)
    plt.plot(casi,rezultati[:,0],label="x="+str(p),color=barve(p/0.05)) #zajci 0 lisice 1

plt.legend(loc="right")
#plt.yscale("log")
plt.title(r"$\ \ \ \ \ \ A(\tau), q=2,p=0.03, (A_0,F_0) = (1,q/p-1+x)$")
plt.savefig("laser/stac8.pdf")
plt.show()
"""






#pp = np.arange(1,5,1)
#rezultati = integriraj((2,2),1,9)
#perioda = rezultati[1]
#rezultati = rezultati[0]
#plt.plot(rezultati[:,0],rezultati[:,1])
#print(perioda)
"""
for p in pp:
    if p<=2:
        iksi = [x*0.5 for x in range(1,21)]
        casi = [integriraj((x,2),p,200)[1] for x in iksi]
    else:
        iksi = [x*0.5 for x in range(1,21)]
        casi = [integriraj((x,2),p,200)[1] for x in iksi]
    plt.plot(iksi,casi,label="p="+str(p),color=barve(p/5)) #zajci 0 lisice 1

#plt.legend(loc="best")
rezultati = integriraj((9,2),4,1000)[0]
plt.plot(rezultati[:,0],rezultati[:,1])
plt.title(r"$p=4, (z_0,l_0) = (9,2), t_{max} = 1000$")
plt.xlabel("Zajci")
plt.ylabel(r"lisice")
plt.savefig("volterra/kaos2.pdf")
plt.show()
"""