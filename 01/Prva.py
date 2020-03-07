# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 20:58:56 2016

@author: Admin
"""
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


def f(t,y,lamb,p): #prvi v y je v, drugi a, funkcija ki jo rešujemo
    return [y[1],-0.2*lamb/(p*(y[1]**(p-2) + (p-2)*y[1]**p))]
def integriraj(t,zac,lamb,p): #propagator
    r = ode(f).set_integrator("dopri5")
    r.set_initial_value([1,zac])
    r.set_f_params(lamb,p)
    return r.integrate(t)[0]
    
def integriraj2(zac,lamb,p): #ipčemo ničlo tega da dobimo zac pogoj
    return 3-integriraj(1,zac,lamb,p)
    
def dobilamb(lamb,zac,p): #iščemo ničlo tega da dobimo lambda
    return 5-quad(integriraj,0,1,(zac,lamb,p))[0]
def fuckcija(x,p): #iščemo nilo tega
    zac, lamb = x
    return (integriraj2(zac,lamb,p),dobilamb(lamb,zac,p))
def iteriraj(p): #poišče ničlo
    zac = 1
    lamb = 200
    koncna = 0
    print("sm hitro")
    zac, lamb = fsolve(fuckcija,(zac,lamb),(p))
    #print(root(fuckcija,(zac,lamb),(p),method="lm"))
    #return [zac,lamb]
        #zac = fsolve(integriraj2,zac,(lamb,p))[0]
        #lamb = fsolve(dobilamb,lamb,(zac,p))
""" ne deluje
def hitrost6(t,lamb,p):
    def pomozna(t,y):
        return f(t,y,lamb,p)
    x = [0,t/2,t,(1+t)/2,1]
    y = [[1,1,1,1,3],[0,0,0,0,0]]
    solution = solve_bvp(pomozna,bc,x,y)
def iteriraj2(p):
    lamb = 200
    
    def fk(x,lambd):
        return (5-quad())
    lamb = fsolve(fk,)
"""
    
def hitrost5(t,p,lamb):
    A = 0.1*lamb/(p*(2*p-1))
    B = 1 + (2*p-1)/(0.2*lamb)*((0.2*lamb/(2*p))**(2*p/(2*p-1)))
    return -(2*p-1)/(0.2*lamb)*(-0.2*lamb/(2*p)*t + A*(2*p-1))**(2*p/(2*p-1)) + B

def posp5(t,p,lamb):
    A = 0.1*lamb/(p*(2*p-1))
    B = 1 + (2*p-1)/(0.2*lamb)*((0.2*lamb/(2*p))**(2*p/(2*p-1)))
    alfa = (2*p-1)/(0.2*lamb)
    beta = 0.2*lamb/(2*p)
    gama = A*(2*p-1)
    delta = B
    return alfa*beta*2*p/(2*p-1)*(-beta * t + gama)**(1/(2*p-1))

def pot5(t,p,lamb):
    A = 0.1*lamb/(p*(2*p-1))
    B = 1 + (2*p-1)/(0.2*lamb)*((0.2*lamb/(2*p))**(2*p/(2*p-1)))
    alfa = (2*p-1)/(0.2*lamb)
    beta = 0.2*lamb/(2*p)
    gama = A*(2*p-1)
    delta = B
    return alfa*(2*p-1)/(beta*(4*p-1))*((-beta*t+gama)**((4*p-1)/(2*p-1)) - gama**((4*p-1)/(2*p-1))) + delta * t    
def najdi(p):
    def fn(lamb):
        return abs(5-quad(hitrost5,0,1,(p,lamb))[0])
    lambd = fsolve(fn,100)
    return lambd




def hitrost3(t,K): #+K K>0
    k = K**0.5
    m = [[1,1,0.1/K],[np.exp(k),np.exp(-k),0.1/K],[(np.exp(k)-1)/k,(1-np.exp(-k))/k,0.1/K]]
    vekt = np.matrix(m)**(-1)*np.matrix([[1],[3],[5]])
    A = float(vekt[0])
    B = float(vekt[1])
    C = float(vekt[2])
    print(C)
    return 0.1*C/K+A*np.exp(k*t)+B*np.exp(-k*t)
def posp3(t,K): #+K K>0
    k = K**0.5
    m = [[1,1,0.1/K],[np.exp(k),np.exp(-k),0.1/K],[(np.exp(k)-1)/k,(1-np.exp(-k))/k,0.1/K]]
    vekt = np.matrix(m)**(-1)*np.matrix([[1],[3],[5]])
    A = float(vekt[0])
    B = float(vekt[1])
    C = float(vekt[2])
    return A*k*np.exp(k*t)-B*k*np.exp(-k*t)
def pot3(t,K): #+K K>0
    k = K**0.5
    m = [[1,1,0.1/K],[np.exp(k),np.exp(-k),0.1/K],[(np.exp(k)-1)/k,(1-np.exp(-k))/k,0.1/K]]
    vekt = np.matrix(m)**(-1)*np.matrix([[1],[3],[5]])
    A = float(vekt[0])
    B = float(vekt[1])
    C = float(vekt[2])
    return 0.1*C/K*t+ A/k*(np.exp(k*t)-1)+B/k*(1-np.exp(-k*t))      
def hitrost4(t,K): #-K K>0
    k = K**0.5
    m = [[0,1,0.1/K],[np.sin(k),np.cos(k),0.1/K],[(1-np.cos(k))/k,np.sin(k)/k,0.1/K]]
    vekt = np.matrix(m)**(-1)*np.matrix([[1],[3],[5]])
    A = float(vekt[0])
    B = float(vekt[1])
    C = float(vekt[2])
    return 0.1*C/K+A*np.sin(k*t)+B*np.cos(k*t)
    
def posp4(t,K): #-K K>0
    k = K**0.5
    m = [[0,1,0.1/K],[np.sin(k),np.cos(k),0.1/K],[(1-np.cos(k))/k,np.sin(k)/k,0.1/K]]
    vekt = np.matrix(m)**(-1)*np.matrix([[1],[3],[5]])
    A = float(vekt[0])
    B = float(vekt[1])
    C = float(vekt[2])
    return A*k*np.cos(k*t)-B*k*np.sin(k*t)
def pot4(t,K): #-K K>0
    k = K**0.5
    m = [[0,1,0.1/K],[np.sin(k),np.cos(k),0.1/K],[(1-np.cos(k))/k,np.sin(k)/k,0.1/K]]
    vekt = np.matrix(m)**(-1)*np.matrix([[1],[3],[5]])
    A = float(vekt[0])
    B = float(vekt[1])
    C = float(vekt[2])
    return 0.1*C/K*t-A/k*(np.cos(k*t)-1)+B/k*np.sin(k*t)    
def hitrost(t,c): #v(t) pri ekstremalni
    #c = t0*v0/l0
    return 1.5*(c-1)/c*t*t - 3*(c-1)/c*t + 1
def pot(t,c): #x(t) pri ekstremalni
    #c = t0*v0/l0
    return 0.5*(c-1)/c*t*t*t - 1.5*(c-1)/c*t*t + t
def posp(t,c):
    return 3*(c-1)/c*t - 3*(c-1)/c
def hitrost2(t,c,vk): #v(t) pri poljubni
    return 3*(c*(vk+1)-2)/c*t*t +(6/c - 2*(vk+2))*t + 1
def pot2(t,c,vk): #x(t) pri poljubni
    return (c*(vk+1)-2)/c*t*t*t +(3/c - (vk+2))*t*t + t
def posp2(t,c,vk):
    return 6*(c*(vk+1)-2)/c*t +6/c - 2*(vk+2)
def hitrost6(t,c,v0):
    return -8*(1/c - v0)*t**2 + 8*(1/c-v0)*t + v0
def posp6(t,c,v0):
    return -16*(1/c-v0)*t + 8*(1/c-v0)
def pot6(t,c,v0):
    return -8/3*(1/c-v0)*t**3 +4*(1/c-v0)*t**2 + v0*t



def hitrost7(t,c):
    lamb = -3/c - 9
    mu = 21/c + 3
    A = mu/8 + 3*lamb/8 - 1/2
    B = 1
    C = -lamb/8 + mu/8 - 1/2
    D = mu/8 + lamb/8 + 1/2
    if t<=1:
        return -lamb/4 *t**2 + A*t + B
    elif t>1:
        tc = t-1
        return -mu/4*tc**2 + C*tc + D
        
def posp7(t,c):
    lamb = -3/c - 9
    mu = 21/c + 3
    A = mu/8 + 3*lamb/8 - 1/2
    B = 1
    C = -lamb/8 + mu/8 - 1/2
    D = mu/8 + lamb/8 + 1/2
    if t<=1:
        return -lamb/2 *t + A
    elif t>1:
        tc = t-1
        return -mu/2*tc + C
        
def pot7(t,c):
    lamb = -3/c - 9
    mu = 21/c + 3
    A = mu/8 + 3*lamb/8 - 1/2
    B = 1
    C = -lamb/8 + mu/8 - 1/2
    D = mu/8 + lamb/8 + 1/2
    if t<=1:
        return -lamb/12 *t**3 + A/2*t**2 + B*t
    elif t>1:
        tc = t-1
        return -mu/12*tc**3 + C/2*tc**2 + D*tc + pot7(1,c)            
    
#ATTEMTPI ZA LIHE
"""
t = np.linspace(0,1,60)
zac, lamb = iteriraj(3)
y=[integriraj(x,zac,lamb,3) for x in t]
plt.plot(t,y,label=r"$p=3$")
#zac, lamb = iteriraj(5)
#y=[integriraj(x,zac,lamb,5) for x in t]
#plt.plot(t,y,label=r"$p=5$")
plt.show()
"""
""" ATTEMPTI ZA SODE
zac, lamb = iteriraj(6)
y=[integriraj(x,zac,lamb,6) for x in t]
plt.plot(t,y,label=r"$p=6$")
zac, lamb = iteriraj(8)
y=[integriraj(x,zac,lamb,8) for x in t]
plt.plot(t,y,label=r"$p=8$")
zac, lamb = iteriraj(10)
y=[integriraj(x,zac,lamb,10) for x in t]

plt.plot(t,y,label=r"$p=10$")
"""
#plt.show()
#plt.plot(t,hitrost2(t,0.2,3))
""" P=2
#plt.plot(t,posp2(t,0.2,0.1),label=r"$v_k=0.1$")
#plt.plot(t,hitrost3(t,50),label=r"$K=50$")
#plt.plot(t,hitrost3(t,30),label=r"$K=30$")
#plt.plot(t,hitrost4(t,50),label=r"$K=-50$")
plt.plot(t,hitrost2(t,0.2,3),label=r"$K\equiv0$")
plt.plot(t,hitrost3(t,10),label=r"$K=10$")
plt.plot(t,hitrost3(t,50),label=r"$K=50$")
plt.plot(t,hitrost3(t,300),label=r"$K=300$")
plt.plot(t,hitrost3(t,600),label=r"$K=600$")
#plt.plot(t,pot4(t,50),label=r"$K=-50$")
#plt.plot(t,posp(t,5),label=r"$c=5$")

plt.plot(t,y)
"""


#ATTEMPTI ZA SODE DELNO ANAL.
t = np.linspace(0,2,100)

#lambd = najdi(0.5)
#y = [hitrost5(x,0.5,lambd) for x in t]
#plt.plot(t,y,label=r"$p=1$")
y = [pot7(x,0.2) for x in t]
plt.plot(t,y,label=r"$x(t)$")
y = [hitrost7(x,0.2) for x in t]
plt.plot(t,y,label=r"$v(t)$")
y = [posp7(x,0.2) for x in t]
plt.plot(t,y,label=r"$a(t)$")
#plt.plot(t,hitrost6(t,0.2,10),label=r"$v(0) = 10$")

plt.legend(loc="lower center")
plt.xlabel(r"$t$")
#plt.ylim([-1,20])
plt.title(r"Dva semaforja $c=0.2,v_k=0$")
plt.savefig("bonus.pdf")
plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    