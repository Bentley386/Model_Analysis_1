# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 10:46:11 2017

@author: Admin
"""
import timeit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.integrate import ode
import matplotlib
import matplotlib.pyplot as plt
from scipy import linalg as lin
from scipy.optimize import fsolve
from scipy.linalg import solve
from scipy.linalg import solve_banded
from scipy.special import jn_zeros #prvi parameter je order, drugi št. ničel
from scipy.special import jv #prvi order drugi argument
from scipy.special import beta
import scipy.special as spec
import scipy.sparse
from scipy.optimize import root
from scipy.integrate import quad
from scipy.integrate import complex_ode
from scipy.optimize import linprog
import scipy.optimize as opt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.linalg import svd
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rc("text",usetex=True)
matplotlib.rcParams["text.latex.unicode"] = True
plt.close("all")
#izmerki = np.loadtxt("ledvice.dat",skiprows=1)
izmerki = np.loadtxt("korozija.txt")
#izmerki = np.loadtxt("farmakoloski.dat")
#x doza y odziv
def hikvadrat(izmerki, model, parametri, napake):
    rezid = (izmerki - model(izmerki,parametri))**2
    return rezid/(napake**2)
def prvip1(x,y0,a):
    return y0*x/(x+a)
def prvip(x,y0,a,p):
    return y0*x**p/(x**p + a**p)
def prvipwat(x,y0,a,p1,p2):
    return y0*x**p1/(x**p2 + a)
def prvip2(x,y0,a,p):
    return y0*x**p/(x**p + a)

def ledviceeno(t,c0,lamb,):
    return c0*np.exp(-lamb*t)
def ledviceenoA(t,c0,lamb,A):
    return c0*np.exp(-lamb*t)+A
def ledvicedvo(t,A,lamb1,B,lamb2):
    return A*np.exp(-lamb1*t) + B*np.exp(-lamb2*t)
def ledvicedvoA(t,A,lamb1,B,lamb2,C):
    return A*np.exp(-lamb1*t)+B*np.exp(-lamb2*t) + C
def ledvicedvo2(t,A,lamb1,B,lamb2):
    return A*np.exp(-lamb1*t**0.5) + B*np.exp(-lamb2*t**0.5)
def ledvicedvoA2(t,A,lamb1,B,lamb2,C):
    return A*np.exp(-lamb1*t**0.5)+B*np.exp(-lamb2*t**0.5) + C

def korozija(u,i0,ua,uc):
    return  i0*(np.exp(u/ua) - np.exp(-u/uc))
def korozija2(u,i0,ua,uc,u0):
    return  i0*(np.exp((u-u0)/ua) - np.exp(-(u-u0)/uc))

yy = izmerki[:,1]    
param,kovar = curve_fit(korozija,izmerki[:,0],izmerki[:,1],(100,100,100),sigma=np.ones(izmerki[:,0].size)*0.0001)    

x = izmerki[:,0]
#x = np.linspace(-100,100,1000)    
y = [korozija(u,param[0],param[1],param[2]) for u in x]

plt.plot(x,yy-y,"--.",label=r"$U_a = U_c = 100$")
param,kovar = curve_fit(korozija2,izmerki[:,0],izmerki[:,1],(100,100,100,10),sigma=np.ones(izmerki[:,0].size)*0.0001)    
y = [korozija2(u,param[0],param[1],param[2],param[3]) for u in x] 
hi2 = (y-yy)**2/(np.ones(x.size)*0.0001)**2
#print(np.sum(hi2)/(x.size-4))
plt.plot(x,yy-y,"--.",label=r"S parametrom $U_0$")
#plt.plot(izmerki[:,0],izmerki[:,1],".")
#plt.errorbar(izmerki[:,0],izmerki[:,1],yerr=0.0001,linestyle="none",color="k")
plt.ylabel("I")
plt.xlabel("U")
plt.title("Korozija - razlike")
plt.legend(loc="best")
plt.savefig("tretja/razlike.pdf")
#plt.show()














































































"""
def narediAiksi(izmerki,napake,funkcije):
    matrika = []
    for i,j in zip(izmerki,napake):
        vrstica = [f(i)/j for f in funkcije]
        matrika.append(vrstica)
    return np.matrix(matrika)
def ita(dekompozicija,b,i):
    return np.sum(dekompozicija[0][:,i]*b)/dekompozicija[1][i] * dekompozicija[2][i]
"""
"""
x = izmerki[:,0]
y = izmerki[:,1]
parametri,kovarianca = curve_fit(ledviceenoA,izmerki[:,0],izmerki[:,1],(100000,0.1,100))

yy = parametri[0]*np.exp(-parametri[1]*x) + parametri[2]
hi2 = np.sum((y-yy)**2)/(x.size-3)
print(hi2)








#en 0.3 malo zlomljeno,  0.2 0.3 malo cudno oba cez 1
x = np.linspace(0,3000,5000)
yy = izmerki[:,1]
parametri,kovarianca = curve_fit(ledvicedvo2,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1))
y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5)
plt.plot(x,y,label=r"$\lambda_{10} = \lambda_{20} = 0.1$")
parametri,kovarianca = curve_fit(ledvicedvo2,izmerki[:,0],izmerki[:,1],(100,0.1,100,3))
y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5)
plt.plot(x,y,label=r"$\lambda_{10} = 0.1, \lambda_{20} = 3$")
#parametri,kovarianca = curve_fit(ledvicedvo2,izmerki[:,0],izmerki[:,1],(100,3,100,3))
#y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5)
#plt.plot(x,y,label=r"$\lambda_{10} = \lambda_{20} = 3$")#parametri,kovarianca = curve_fit(ledvicedvoA2,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1,10))
parametri,kovarianca = curve_fit(ledvicedvoA2,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1,10))
y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5) + parametri[4]
plt.plot(x,y,label=r"Z aditivno konstanto")
plt.plot(izmerki[:,0],izmerki[:,1],".")


parametri,kovarianca = curve_fit(ledvicedvo2,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1))
y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5)
plt.plot(x,yy-y,"--.",label=r"Korenski čas $\lambda_{10} = \lambda_{20} = 0.1$")
parametri,kovarianca = curve_fit(ledvicedvo2,izmerki[:,0],izmerki[:,1],(100,0.1,100,3))
y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5)
plt.plot(x,yy-y,"--.",label=r"Korenski čas $\lambda_{10} = 0.1, \lambda_{20} = 3$")
parametri,kovarianca = curve_fit(ledvicedvoA2,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1,10))
y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5) + parametri[4]
plt.plot(x,yy-y,"--.",label=r"Korenski čas z aditivno konstanto")
parametri,kovarianca = curve_fit(ledvicedvo,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1))
y = parametri[0]*np.exp(-parametri[1]*x) + parametri[2]*np.exp(-parametri[3]*x)
plt.plot(x,yy-y,"--.",label=r"$\lambda_{10} = \lambda_{20} = 0.1$")
parametri,kovarianca = curve_fit(ledvicedvoA,izmerki[:,0],izmerki[:,1],(100,0.1,100,0.1,10))
y = parametri[0]*np.exp(-parametri[1]*x) + parametri[2]*np.exp(-parametri[3]*x) + parametri[4]
plt.plot(x,yy-y,"--.",label=r"Z aditivno konstanto")
plt.show()


plt.legend(loc="best")
plt.title("Rešitve za dvorazdelčni model")
plt.ylabel("c(t)")
plt.xlabel("t")
plt.savefig("druga/amb12.pdf")
plt.show()
"""
"""
for i in np.linspace(0.1,5,10):
    parametri,kovarianca = curve_fit(ledvicedvo2,izmerki[:,0],izmerki[:,1],(100,0.1,100,i))
    y = parametri[0]*np.exp(-parametri[1]*x**0.5) + parametri[2]*np.exp(-parametri[3]*x**0.5)
    plt.plot(x,y,label=r"$\lambda_{10} = \lambda_{20} = 0.1$")
plt.show()
"""



#mejna lamb = 0.2
"""
x = izmerki[:,0]
y = izmerki[:,1]
parametri,kovarianca = curve_fit(ledviceeno,izmerki[:,0],izmerki[:,1],(10000,0.2))
c0 = parametri[0]
lamb = parametri[1]
plt.plot(x,(y-parametri[0]*np.exp(-parametri[1]*x)),"--.",label=r"$A=0$")
parametri,kovarianca = curve_fit(ledviceenoA,izmerki[:,0],izmerki[:,1],(10000,0.2,100))
c0 = parametri[0]
lamb = parametri[1]
A = parametri[2]
hi2 = (izmerki[:,1]-(c0*np.exp(-lamb*izmerki[:,0])+A))**2
print(np.sum(hi2)/(izmerki[:,0].size - 3))
plt.plot(x,(y- parametri[0]*np.exp(-parametri[1]*x)-parametri[2]),"--.",label=r"$A \neq 0$")

plt.plot(izmerki[:,0],izmerki[:,1],".")
plt.title("Rešitve za enodelčni model - ostanki")
plt.xlabel("t")
plt.ylabel("c(t)")
plt.legend(loc="best")
plt.savefig("druga/lamb13.pdf")
plt.show()
"""































"""
funkcije = [lambda x: 1, lambda x: -x]
A = narediAiksi(izmerki[:,0],np.ones(izmerki[:,0].size),funkcije)
b = np.log(izmerki[:,1])
dekompozicija = svd(A,full_matrices=False)    
suma = np.zeros(2)
for i in range(2):
    suma+= ita(dekompozicija,b,i)
c0 = np.exp(suma[0])
lamb= suma[1]
hi2 = (izmerki[:,1]-c0*np.exp(-lamb*izmerki[:,0]))**2
print(np.sum(hi2)/(izmerki[:,0].size - 2))
plt.plot(x,c0*np.exp(-lamb*x),label=r"Metoda SVD")
"""








































    
    
    












































"""

    
plt.title("Razlike med modeli in izmerki")
   

x = izmerki[:,0]
param = curve_fit(prvip,izmerki[:,0],izmerki[:,1],p0=(100,20,1),sigma = np.ones(8)*3)[0]
y = param[0]*x**param[2]/(x**param[2] + param[1]**param[2])
plt.plot(izmerki[:,0],izmerki[:,1]-y,"r.")
plt.legend(loc="best")
plt.savefig("prva/razlike.pdf")
plt.show()
"""


#param,kovar = curve_fit(prvipwat,izmerki[:,0],izmerki[:,1],p0=(100,20,1,1),sigma=np.ones(8)*3)    
#hi2 = np.sum((izmerki[:,1]-(izmerki[:,0]**param[2]*param[0])/(izmerki[:,0]**param[3]+param[1]))**2/(np.ones(8)*3)**2) 
#param2,kovar2 = curve_fit(prvip2,izmerki[:,0],izmerki[:,1],p0=(100,20,1),sigma=np.ones(8)*3)    
#hi2 = np.sum((izmerki[:,1]-(izmerki[:,0]**param2[2]*param2[0])/(izmerki[:,0]**param2[2]+param2[1]))**2/(np.ones(8)*3)**2)    