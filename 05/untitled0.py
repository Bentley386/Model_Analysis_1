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

plt.rc("text",usetex=True)
plt.close("all")
# Plot the surface.
def ura(t,x,a,c): #x je i2, i,s23,s28
    x1 = a*x[3]*x[1] - c*x[0]*x[2]
    x2 = -2*a*x[3]*x[1] + 2*c*x[2]*x[0]
    x3 = -2*c*x[2]*x[0]
    x4 = -a*x[3]*x[1]
    return [x1,x2,x3,x4]
def integrirajura(zacetni,a,c,cas):
    r = ode(ura).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(a,c)
    rezultati = []
    inicialc = zacetni[0]
    perioda=0
    smoze = True
    for t in np.linspace(0.001,cas,500):
        trenutno = r.integrate(t)
        if trenutno[0] > 0.5 and smoze and t>5:
            perioda = t
            smoze = False
        rezultati.append(trenutno)
    return [np.asarray(rezultati),perioda]
def bineks(t,x,p,q,r): #binarna eksaktno, x = [A,A*,B=r*Azvezda]
    return [-p*x[0]*x[0] + q*x[0]*x[1],p*x[0]*x[0]-r*x[1]-q*x[0]*x[1],r*x[1]]
def bineksjac(t,x,p,q,r):
    return [[-p*x[0]+2*q*x[1],2*q*x[0],0],[2*p*x[0],-r,0],[0,r,0]]
def binstac(t,x,p,q,r): #isto stac. x =[A,B=pA^2]
    return [-p*x[0]**2 + q*p*x[0]**3/(r+q*x[0]),r*p*x[0]**2/(r+q*x[0])]
def binstacanal(x,t,p,q,r,A):
    return 4*q/r*(np.log(x)-np.log(r-3*q*x))-1/x + p*t/2 - A
def binstacanal2(A,t,p,q,r,x):
    return 4*q/r*(np.log(x/(r-3*q*x)))-1/x + p*t/2 - A
def binstacjac(t,x,p,q,r):
    return [[-p*x[0]+6*q*p/r*x[0]**2,0],[2*p*x[0],0]]
def drugaeks(cas,x,p,q,r,s,t): #x=[H2,Br2,HBr,Br,H]
    x1 = s*x[2]*x[4] - r*x[3]*x[0]
    x2 = q*x[3]**2 - p*x[1] - t*x[4]*x[1]
    x3 = r*x[3]*x[0] -s*x[2]*x[4]+t*x[4]*x[1]
    x4 = p*x[1]-q*x[3]**2 + t*x[4]*x[1]+s*x[2]*x[4]-r*x[3]*x[0]
    x5 = r*x[3]*x[0] - s*x[4]*x[2] - t*x[4]*x[1]
    return [x1,x2,x3,x4,x5]
def drugastac(cas,x,p,q,r,s,t): #x = h2,br2,hbr
    br = np.sqrt(p/q*x[1])
    h = r*br*x[0]/(s*x[2]+t*x[1])
    x1 = s*x[2]*h - r*br*x[0]
    x2 = q*br**2 - p*x[1] - t*h*x[1]
    x3 = r*br*x[0]-s*x[2]*h+t*h*x[1]
    return [x1,x2,x3]  

def integrirajstac2(zacetni,p,q,rr,s,tt,cas):
    r = ode(drugastac).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p,q,rr,s,tt)
    rezultati = []
    for t in np.linspace(0.001,cas,5000):
        trenutno = r.integrate(t)
        rezultati.append(trenutno)
    return np.asarray(rezultati)  
def integriraj2(zacetni,p,q,rr,s,tt,cas):
    r = ode(drugaeks).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p,q,rr,s,tt)
    rezultati = []
    for t in np.linspace(0.001,cas,5000):
        trenutno = r.integrate(t)
        rezultati.append(trenutno)
    return np.asarray(rezultati)    
def integriraj(zacetni,p,q,rr,cas,tol=0.1):
    r = ode(bineks).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p,q,rr)
    rezultati = []
    for t in np.linspace(0.001,cas,10000):
        trenutno = r.integrate(t)
        rezultati.append(trenutno)
    return np.asarray(rezultati)
def ravnovesje(zacetni,p,q,rr,cas,tol=0.1):
    r = ode(bineks).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p,q,rr)
    rezultati = []
    for t in np.linspace(0.001,cas,10000):
        trenutno = r.integrate(t)
        if abs(-p*trenutno[0]**2 + q*trenutno[0]*trenutno[1]) < 0.005 and t>1:
            print(t)
            return t
    return 0
def integrirajstac(zacetni,p,q,rr,cas,tol=0.1):
    r = ode(binstac).set_integrator("dopri5")
    r.set_initial_value(zacetni,0).set_f_params(p,q,rr)
    rezultati = []
    for t in np.linspace(0.001,cas,10000):
        trenutno = r.integrate(t)
        rezultati.append(trenutno)
    return np.asarray(rezultati)
def nicla(x,h2,br2,hbr,odvod): #x=[k,m]
    rez= x[0]*h2*np.sqrt(br2)/(x[1]+hbr/br2)-odvod
    return (rez,0)

barve = cm.get_cmap("autumn")
barve2 = cm.get_cmap("copper")
p,q,r,s,t = 1,1,1,1,2.5
tajm = 10
x = np.linspace(0.001,tajm,5000) #x=[H2,Br2,HBr,Br,H] m=t/s=2.5. h2/br2 = 100,1,0.01
"""
for i in np.linspace(1,10,10):
    rezultati = integrirajstac2((1,1,0),i,q,r,s,t,tajm)    
    #rezultati = integriraj2((1,1,0,0,0),p,q,r,s,i,tajm)    
    if i==1:
        plt.plot(x,rezultati[:,2],label=r"$HBr,\  p=1$",color=barve2(i/10),zorder=10)
    elif i==10:
        plt.plot(x,rezultati[:,2],label=r"$HBr,\  p=100$",color=barve2(i/10),zorder=10)
    else:
        plt.plot(x,rezultati[:,2],color=barve2(i/10),zorder=10)        
    if i==1:
        plt.plot(x,rezultati[:,0],label=r"$H_2,\  p=1$",color=barve(i/10))
    elif i==10:
        plt.plot(x,rezultati[:,0],label=r"$H_2,\  p=100$",color=barve(i/10))
    else:
        plt.plot(x,rezultati[:,0],color=barve(i/10))
"""
tajm = 1000
x = np.linspace(0.001,tajm,500)
casi = []

for i in range(5,21):
    print(i)
    casi.append(integrirajura((1,0,i,10),0.01,1000,tajm)[1])
plt.plot(range(5,21),casi)
"""
rezultati = integrirajura((1,0,20,10),0.01,1000,tajm)[0]
plt.plot(x,rezultati[:,0])
"""
plt.ylabel("Cas do reakcije")
plt.xlabel(r"Konc. tiosulfata na zacetku") #i2 i s2
plt.title(r"$c/a=100000, (I_2,I^-)(0)=(1,0)$")
#plt.legend(loc="best")
plt.savefig("druga/hitrost.pdf")
plt.show()









#odvodistac = [drugastac(x[i],(rezultati[i,0],rezultati[i,1],rezultati[i,2]),p,q,r,s,t)[2] for i in range(len(x))]
#kji = [fsolve(nicla,[1,1],(rezultati[i,0],rezultati[i,1],rezultati[i,2],odvodistac[i]))[0] for i in range(len(x))]
#mji = [fsolve(nicla,[1,1],(rezultati[i,0],rezultati[i,1],rezultati[i,2],odvodistac[i]))[1] for i in range(len(x))]
#plt.plot(x,kji,label="k pri stac. approx")
#plt.plot(x,mji,label="m pri stac. approx") 
#rezultati = integriraj2((1,1,0,0,0),p,q,r,s,t,tajm)
#odvodieks = [drugaeks(x[i],(rezultati[i,0],rezultati[i,1],rezultati[i,2],rezultati[i,3],rezultati[i,4]),p,q,r,s,t)[2] for i in range(len(x))]
#kji = [fsolve(nicla,[1,1],(rezultati[i,0],rezultati[i,1],rezultati[i,2],odvodieks[i]))[0] for i in range(len(x))]
#mji = [fsolve(nicla,[1,1],(rezultati[i,0],rezultati[i,1],rezultati[i,2],odvodieks[i]))[1] for i in range(len(x))]
#plt.plot(x,kji,label="k pri eksaktni integraciji")
#plt.plot(x,mji,label="m pri eksaktni integraciji")
#plt.plot(x,[5 for i in x],label="pravi k")
#plt.plot(x,[1 for i in x],label="pravi m")
#rezultati = integriraj2((1,1,0,np.sqrt(p/q),r*np.sqrt(p/q)),i,q,r,s,t,tajm)
#plt.plot(x,rezultati[:,0])
#rezultati = integrirajstac2((1,1,0),i,q,r,s,t,tajm)
#plt.plot(x,rezultati[:,0])
#plt.plot(x,rezultati[:,2],label="hbr")
#plt.plot(x,rezultati[:,2])



#plt.plot(cas,abs(rezultati1[:,0]-rezultati2[:,0])/rezultati1[:,0],color=barve(0),label=r"$r/q A(0) = 0.1$")
#plt.plot(cas,rezultati[:,0],label="eksaktno A(t)")
#plt.plot(cas,rezultati[:,1],label="eksaktno A*(t)")
#plt.plot(cas,rezultati[:,2],label="eksaktno B(t)")
#rezultati = integrirajstac((10,0),0.1,100,100,tajm)
#plt.plot(cas,rezultati[:,0],"--",label="stac. approx. A(t)")
#plt.plot(cas,[0.1*x**2/(100+100*x) for x in rezultati[:,0]],label="stac. approx. A*(t)")
#plt.plot(cas,rezultati[:,0],label="stac. approx. A(t)")
#plt.plot(cas,rezultati[:,1],"--",label="stac. approx. B(t)")










#rezultati1 = integriraj((1,0,0),0.1,100,100,50)
#rezultati2 = integrirajstac((1,0),0.1,100,100,50)
#plt.plot(cas,abs(rezultati1[:,0]-rezultati2[:,0])/rezultati1[:,0],color=barve(0.15),label=r"$q/p = 1000$")
#rezultati1 = integriraj((1,0,0),1,100,100,50)
#rezultati2 = integrirajstac((1,0),1,100,100,50)
#plt.plot(cas,abs(rezultati1[:,0]-rezultati2[:,0])/rezultati1[:,0],color=barve(0.3),label=r"$q/p = 100$")
#rezultati1 = integriraj((1,0,0),100,100,100,50)
#rezultati2 = integrirajstac((1,0),100,100,100,50)
#plt.plot(cas,abs(rezultati1[:,0]-rezultati2[:,0])/rezultati1[:,0],color=barve(0.45),label=r"$q/p=1$")
#rezultati1 = integriraj((1,0,0),10000,100,100,50)
#rezultati2 = integrirajstac((1,0),10000,100,100,50)
#plt.plot(cas,abs(rezultati1[:,0]-rezultati2[:,0])/rezultati1[:,0],color=barve(0.6),label=r"$q/p = 0.01$")
#rezultati1 = integriraj((1,0,0),100000,100,100,50)
#rezultati2 = integrirajstac((1,0),100000,100,100,50)
#plt.plot(cas,abs(rezultati1[:,0]-rezultati2[:,0])/rezultati1[:,0],color=barve(0.75),label=r"$q/p = 0.001$")
#rezultati = integriraj((10,0,0),0.1,100,100,tajm) #p q r
#rezultati2 = integrirajstac(1,0.1,1,0.1,50)