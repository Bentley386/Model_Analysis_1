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
pi = np.pi
def timer(code):
    start_time = timeit.default_timer()
    exec(code)
    elapsed = timeit.default_timer() - start_time
    return elapsed

def volumenkonst1(N):
    koliko = 0
    for i in range(N):
        x = np.random.rand() - 0.5
        y = np.random.rand() - 0.5
        z = np.random.rand() - 0.5
        if ((x*x + y*y) < 0.25 and (x*x + z*z) < 0.25 and (z*z + y*y) < 0.25):
            koliko += 1
    return koliko/N
def volumenkonst2(N):
    koliko = 0
    for i in range(N):
        fi = np.random.rand()*2*pi
        theta = np.arccos(2*np.random.rand()-1)
        r = np.random.rand()**(1/3)
        if (r*np.sin(theta) < 1 and r*r*(np.cos(fi)**2 *np.sin(theta)**2 + np.cos(theta)**2) < 1 and r*r*(np.sin(fi)**2 *np.sin(theta)**2 + np.cos(theta)**2) < 1):
            koliko += 1
    return koliko/N
def generirajpresek(N):
    st = 0
    stevila = []
    while(st<N):
        x,y,z = np.random.rand(3) - 0.5
        if ((x*x + y*y) < 0.25 and (x*x + z*z) < 0.25 and (z*z + y*y) < 0.25):        
            st+=1
            stevila.append((x,y,z))
    return stevila
def vztrajnostnikonst(N):
    stevila = np.asarray(generirajpresek(N))
    return np.sum(stevila[:,0]**2 + stevila[:,1]**2)*0.58/N
def volumenr(N,p):
    r0 = np.sqrt(3/2)*0.5
    stevila = np.asarray(generirajpresek(N))
    volum = np.sum((np.sqrt(stevila[:,0]**2+stevila[:,1]**2 + stevila[:,2]**2)/r0)**p) / N
    volum2 =np.sum((np.sqrt(stevila[:,0]**2+stevila[:,1]**2 + stevila[:,2]**2)/r0)**(2*p)) / N
    return [0.586*volum,0.586/np.sqrt(N)*np.sqrt(volum2-volum**2)]
def vztrajnostnir(N,p):
    r0 = np.sqrt(3/2)*0.5
    stevila = np.asarray(generirajpresek(N))
    volum = np.sum((np.sqrt(stevila[:,0]**2+stevila[:,1]**2 + stevila[:,2]**2)/r0)**p * (stevila[:,0]**2+stevila[:,1]**2)) / N
    volum2 =np.sum((np.sqrt(stevila[:,0]**2+stevila[:,1]**2 + stevila[:,2]**2)/r0)**(2*p)*(stevila[:,0]**2+stevila[:,1]**2)**2) / N
    return [0.586*volum,0.586/np.sqrt(N)*np.sqrt(volum2-volum**2)]
#x = np.linspace(100,100000,100)
#y = [vztrajnostnikonst(int(i)) for i in x]

def delez(N,k=1):
    stevilo = 0
    for i in range(N):
        cos = 2*np.random.rand() - 1
        r = np.random.rand()**(1/3)
        x = -k*np.log(1-np.random.rand())
        if x>=(-r*cos + np.sqrt(r*r*cos*cos + 1 - r*r)):
            stevilo +=1
    return stevilo/N
def stsipanj():
    stevilo = 1
    lok = -0.5*np.log(1-np.random.rand())
    while lok>0 and lok<1:
        plus = np.random.rand()
        if plus>0.5:
            lok += -0.5*np.log(1-np.random.rand())
            stevilo += 1
        else:
            lok -= -0.5*np.log(1-np.random.rand())
            stevilo +=1
    return stevilo
def stizo():
    stevilo = 1
    lok = -0.5*np.log(1-np.random.rand())
    while lok>0 and lok<1:
        cos = 2*np.random.rand()-1
        lok += (-0.5*np.log(1-np.random.rand()))*cos
        stevilo += 1
    return stevilo

def izosipanja(d=1):
    stevilo = 1
    lok = -0.5*np.log(1-np.random.rand())
    while lok>0 and lok<d:
        cos = 2*np.random.rand()-1
        lok += (-0.5*np.log(1-np.random.rand()))*cos
        stevilo += 1
    if lok>=d:
        return 1
    else:
        return 0
def sipanja():
    stevilo = 1
    lok = -0.5*np.log(1-np.random.rand())
    while lok>0 and lok<1:
        plus = np.random.rand()
        if plus>0.5:
            lok += -0.5*np.log(1-np.random.rand())
            stevilo += 1
        else:
            lok -= -0.5*np.log(1-np.random.rand())
            stevilo +=1
    if lok>=1:
        return 1
    else:
        return 0
    #return stevilo
pi = np.pi
def izo():
    lok = -0.5*np.log(1-np.random.rand())
    if lok>=1:
        return 0
    while lok>0 and lok<1:
        cos = (2*np.random.rand()-1)
        razdalja = -0.5*np.log(1-np.random.rand())
        lok += razdalja*cos 
    return np.arccos(cos)
        
def prepustnost(N):
    prepusceni = 0
    for i in range(N):
        prepusceni+=sipanja()
    return prepusceni/N    
def izoprepustnost(N,d=1):
    prepusceni = 0
    for i in range(N):
        prepusceni+=izosipanja(d)
    return prepusceni/N     
def model(x,A,B,C):
    return A*x**B + C
    
#koti = [izo() for i in range(1000000)]    
if 0:
    dji = np.linspace(0.1,10,100)
    prepustnosti = [izoprepustnost(10000,d) for d in dji]
    plt.plot(dji,prepustnosti)
    plt.title("Prepustnost za izotropen primer")
    plt.xlabel("d")
    plt.ylabel("Prepustnost")
    plt.savefig("tretja/debelina2.pdf")
    plt.show()
if 0:
    vrednosti = np.asarray([izoprepustnost(10000) for i in range(200)])
    povprecje = np.sum(vrednosti)/200
    napaka = np.sqrt(np.sum((vrednosti - povprecje)**2)/200)
    print(povprecje)
    print(napaka)
if 0:
    nji = np.linspace(100,100000,100)    
    nji = [int(i) for i in nji]
    prepustnosti = [izoprepustnost(i) for i in nji]
    plt.plot(nji,prepustnosti)
    plt.title("Prepustnosti")
    plt.ylabel("Prepustnost")
    plt.xlabel("N")
    plt.savefig("tretja/izoprepustnost.pdf")
    plt.show()  
if 0:
    #bini = [0.1*i for i in range(1,101)] + [i for i in range(11,26)]
    #bini = [i for i in range(1,21)] + [i for i in range(21,60)]
    #bini2 = [0.1*i + 0.05 for i in range(1,101)] + [i+0.5 for i in range(11,26)]
    #bini2 = [i + 0.5 for i in range(1,21)] + [i+0.5 for i in range(21,60)]
    bini = np.linspace(0.01,pi,300)
    #stevila = [stizo() for i in range(1000000)]
    histogram = np.histogram(koti,bini,normed=True)
    plt.plot(bini[:-1],histogram[0])
    plt.title("Kotna porazdelitev, N=1000000")
    plt.xlabel("Kot [Rad]")    
    #plt.hist(stevila,bins=20,normed=True)
    plt.savefig("tretja/kotna.pdf")
    plt.show()

print(delez(10000000,1))
"""    
nji = np.linspace(0.1,10,100)
delezi = [delez(10000,i) for i in nji]    
parametri = curve_fit(model,nji,delezi,(1,0.5,0.1))[0]
fit = [model(x,parametri[0],parametri[1],parametri[2]) for x in nji]
plt.plot(nji,fit)
print(parametri)
plt.title(r"Delež fotonov, ki uide iz krogle")
plt.xlabel(r"N")
plt.ylabel("Delež")
plt.plot(nji,delezi)
#plt.savefig("druga/delezi.pdf")
plt.show()
"""


if 0:
    pji = np.linspace(0.5,5,100)
    volumni = [volumenr(10000,i) for i in pji]
    fig, ax1 = plt.subplots()
    ax1.plot(pji,[i[0] for i in volumni],"b")
    ax2 = ax1.twinx()
    ax2.plot(pji,[i[1] for i in volumni],"r")
    ax1.set_xlabel("p")
    ax1.set_ylabel("Masa",color="b")
    ax1.tick_params("y",colors="b")
    ax2.set_ylabel("Napaka",color="r")
    ax2.tick_params("y",colors="r")
    #ax2.set_yscale("log")
    #plt.ylabel("Volumen")
    plt.title(r"Mase preseka treh valjev, $\rho = (r/r_0)^p$")
    plt.tight_layout()
    
    plt.savefig("prva/volumenp.pdf")
    plt.show()
