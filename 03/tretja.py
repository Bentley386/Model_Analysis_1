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

ch = np.cosh(1)
sh = np.sinh(1)

def elektro(xx,n,m): #0,2,4,.... fiji  1,3,5... thete
    x = xx[0:len(xx):2]
    x = np.insert(x,0,0)
    y = xx[1:len(xx):2]
    y = np.insert(y,0,0)
    vsota = 0
    faktorji = [1]*len(x)
    for i in range(n):
        faktorji[i]=2
    for i in range(n,n+m):
        faktorji[i]=10
    #x je fi, y so thete  thet je en več
    for i in range(len(x)): #i=0,1,2,3..
        for j in range(i+1,len(x)): #j=1,2,3...
            fi1 = x[i]
            fi2 = x[j]
            th1 = y[i]
            th2 = y[j]
            medvs = (1-np.sin(th1)*np.sin(th2)*np.cos(fi1-fi2)-np.cos(th1)*np.cos(th2))**(-0.5)
            vsota+= medvs*faktorji[i]*faktorji[j]
    return vsota            

def elektro2(xx,n,m): #0,2,4,.... fiji  1,3,5... thete
    x = xx[0:len(xx):2]
    #x = np.insert(x,0,0)
    y = xx[1:len(xx):2]
    #y = np.insert(y,0,0)
    vsota = 0
    faktorji = [1]*len(x)
    for i in range(n):
        faktorji[i]=2
    for i in range(n,n+m):
        faktorji[i]=10
    #x je fi, y so thete  thet je en več
    for i in range(len(x)): #i=0,1,2,3..
        for j in range(i+1,len(x)): #j=1,2,3...
            fi1 = x[i]
            fi2 = x[j]
            th1 = y[i]
            th2 = y[j]
            medvs = (ch**2 *(np.cos(th1)**2+np.cos(th2)**2-2*np.cos(th1)*np.cos(th2)*np.cos(fi1-fi2)) + sh**2*(np.sin(th1)-np.sin(th2))**2)**(-0.5)
            vsota+= medvs*faktorji[i]*faktorji[j]
    return vsota      
    
    
    
    
"""
fiji = np.random.rand(9)
thete = np.random.rand(9)
oboje = []
for i in range(len(fiji)):
    oboje.append(fiji[i])
    oboje.append(thete[i])
oboje = np.array(oboje) 
print(opt.minimize(elektro,oboje,method="Nelder-Mead",args=(0,0),options={"maxiter":500000,"maxfev":500000}))    
"""
def timer(code):
    start_time = timeit.default_timer()
    exec(code)
    elapsed = timeit.default_timer() - start_time
    return elapsed
    
def semafor(xx,vmax,alfa,lamb):
    N = len(xx) #v0=vN = 1
    delta = 1/(N+1)
    #x = np.insert(xx,0,1)
    #x = np.append(x,1)
    f1 = 0.5*((xx[0]-1)/delta)**2 + 0.5*((3-xx[-1])/delta)**2 #tle
    for i in range(len(xx)-1):
        f1 += ((xx[i+1]-xx[i])/delta)**2
    f2=0
    for i in xx:
        f2 += np.exp((alfa*(i-vmax)))
    #f4 = 0
    #for i in xx:
        #f4 += np.exp(alfa*(0.5-i))
    #f3 = np.exp(-beta*(1+sum(xx)-1/delta)**2) #vez da se ustavi pri L=1
    #f4 = np.exp(beta*(-2-sum(xx)))#nad vmin=0
    f5 = -lamb*delta*(0.5+sum(xx))
    return f1+f2 +f5
    

def najdiniclo(lamb,alfa):
    hitrosti = opt.minimize(semafor,np.array([1]*50),method="Powell",args=(4,alfa,lamb)).x
    rezultat = 5.5 - 1/51*(0.5+sum(hitrosti))
    return rezultat
def naredi(alpha):
    lamb = fsolve(najdiniclo,1,args=(alpha),maxfev=10000)
    x = opt.minimize(semafor,np.array([1]*50),method="Powell",args=(4,alpha,lamb)).x
    x = np.insert(x,0,1)
    x = np.append(x,3)
    return x
#SEMAFOR
#x=naredi(0.5)
#plt.plot(np.linspace(0,1,len(x)),x,label=r"$\alpha=0.5$")    
"""
x=naredi(1)
plt.plot(np.linspace(0,1,len(x)),x,label=r"$\alpha=1$")
x=naredi(2)
plt.plot(np.linspace(0,1,len(x)),x,label=r"$\alpha=2$")
x=naredi(3)
plt.plot(np.linspace(0,1,len(x)),x,label=r"$\alpha=3$")
x=naredi(4)
plt.plot(np.linspace(0,1,len(x)),x,label=r"$\alpha=4$")
x=naredi(5)
plt.plot(np.linspace(0,1,len(x)),x,label=r"$\alpha=5$")
plt.title(r"$L=5.5, v_{max}=6,v_{N}=3$")
plt.ylabel("v")
plt.xlabel("t")
plt.legend(loc="best")
plt.savefig("semafor3.pdf")
plt.show()
"""
    
#THOMSON

x = np.random.rand(6) #24
n = 0
m = 0
x=opt.minimize(elektro,x,method="powell",args=(n,m))
print(x.success)
E = x.fun

x = x.x
x = np.insert(x,[0,0],[0,0])
fi = x[0:len(x):2]
th = x[1:len(x):2]
zz = np.cos(th)
yy = np.sin(th)*np.sin(fi)
xx = np.sin(th)*np.cos(fi)  


fig = plt.figure()
ax = fig.gca(projection='3d')
# Make data.
u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:200j] #u fi v theta
x = np.sin(v)*np.cos(u)
y = np.sin(v)*np.sin(u)
z = np.cos(v)
ax.plot_surface(x, y, z,alpha=0.15,linewidth=0.5,color="r")
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
ax.text(0.6,0.6,1.1,"N=4")
#pr = ax.scatter(xx[:n],yy[:n],zz[:n],color="blue")  
#pr = ax.scatter(xx[n:n+m],yy[n:n+m],zz[n:n+m],color="red")
#pr = ax.scatter(xx[n+m:],yy[n+m:],zz[n+m:],color="k") 
pr = ax.scatter(xx[n:],yy[n:],zz[n:],color="k") 
ax.set_aspect("equal")
ax.dist = 10.4
plt.savefig("2.pdf")

# Plot the surface.
"""
def minimiziraj(zacetni):
    x = opt.minimize(elektro,zacetni,method="Powell",args=(0,0))
    print(x.success)
    return x.fun
def minimiziraj2(zacetni):
    x = opt.minimize(elektro,zacetni,method="Nelder-Mead",args=(0,0),options={"maxiter":500000,"maxfev":500000})
    print(x.success)
    return x.fun
def minimiziraj3(zacetni):
    x = opt.minimize(elektro2,zacetni,method="Powell",args=(0,0))
    print(x.success)
    return x.fun
nji = range(2,20)
priblizki = [np.random.rand((i-1)*2) for i in nji]
priblizki2=[np.random.rand(2*i) for i in nji]
energije = [minimiziraj(i)/(np.sqrt(2)) for i in priblizki]
#casi = [timer("minimiziraj({0})".format(i.tolist())) for i in priblizki]
plt.plot(nji,energije,"ro",label="Krogla")
plt.xlim(2,20)
plt.xlabel("N")
plt.ylabel("E")
#plt.tick_params("y",colors="r")
energije2 = [minimiziraj3(i) for i in priblizki2]
#ax2 = plt.twinx()
#casi2 = [timer("minimiziraj2({0})".format(i.tolist())) for i in priblizki]
plt.plot(nji,energije2,"bo",label="Elipsoid")
plt.legend(loc="best")
#plt.tick_params("y",colors="b")
#ax2.set_ylabel("E(Nelder-Mead)",color="b")
plt.savefig("tom4.pdf")
plt.show()
"""
