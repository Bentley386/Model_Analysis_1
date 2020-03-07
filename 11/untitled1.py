# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 10:46:11 2017

@author: Admin
"""
import timeit
import matplotlib.animation as animation
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
#from scipy.special import beta
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





def izumiraj(N,betas,deltat,betar=0):
    populacija = [N]
    trenutna = N
    cas = [0]
    while trenutna > 0:
        if betar:
            trenutna += np.random.poisson(betar*deltat*trenutna) - np.random.poisson(betas*deltat*trenutna)
        else:
            trenutna -= np.random.poisson(betas*deltat*trenutna)
        cas.append(cas[-1]+deltat)
        if trenutna <= 0:
            populacija.append(0)
            return (cas,populacija)
        populacija.append(trenutna)        
    print("wtf")
    return (cas,populacija)    
def casizumrtja(N,betas,deltat,betar=0):
    cas = 0
    trenutna = N
    while trenutna > 0:
        if betar:
            trenutna += np.random.poisson(betar*deltat*trenutna) - np.random.poisson(betas*deltat*trenutna)            
        else:
            trenutna -= np.random.poisson(betas*deltat*trenutna)
        cas += deltat
    return cas
def povpinvar(funk,N,betas,deltat,n,betar):
    if betar:
        sample = np.asarray([funk(N,betas,deltat,betar) for i in range(n)])
    else:
        sample = np.asarray([funk(N,betas,deltat) for i in range(n)])
    povp = np.sum(sample)/n
    var = np.sum(sample-povp)/(n-1)
    return (povp,var)
#podatki = izumiraj(250,1,0.1)
def naredimatriko(betar,betas,delta,N):
    def funkcija(x,y):
        if x==y:
            return 1-(betar+betas)*delta*x
        elif x==(y+1):
            return betar*y*delta
        elif y==(x+1):
            return delta*betas*y
        else:
            return 0
    matrika = np.matrix(np.fromfunction(np.vectorize(funkcija),(N,N)))
    return matrika
def averagiraj(zacetni):
    iksi = np.asarray(range(len(zacetni)))
    return np.sum(iksi*zacetni)
def variraj(zacetni):
    iksi = np.asarray(range(len(zacetni)))
    def gaussovka(x,mu,sigma):
        return 1/(2*pi*sigma**2)**0.5 * np.exp(-(x-mu)**2/(2*sigma**2))
    parametri = curve_fit(gaussovka,iksi,np.array(zacetni.T).ravel(),(np.argmax(np.array(zacetni.T).ravel()),1))[0]
    return parametri[1]
def prehodna(N,betar,betas,delta,povprecje=False,spremljaj=False):
    if povprecje:
        povprecjee = []
        cass = []
    if spremljaj:
        spremljajoce = []
    cas = 0
    counter = 0
    if N==250:
        zacetni = [0]*250
        zacetni[249] = 1
        zacetni = np.matrix(zacetni).T
        matrika = naredimatriko(betar,betas,delta,250)
    elif N==25:
        zacetni = [0]*25
        zacetni[24] = 1
        zacetni = np.matrix(zacetni).T
        matrika = naredimatriko(betar,betas,delta,25)
    #zacetni = [0]*N
    #zacetni[-1] = 1
    #zacetni = np.matrix(zacetni).T
    #matrika = naredimatriko(betar,betas,delta,N)
    while zacetni[0]<0.99:
        if spremljaj and counter%250==0:
            spremljajoce.append(zacetni)
        if povprecje == True and counter%1000 == 0:
            cass.append(cas)
            povprecjee.append(variraj(zacetni))
        zacetni = matrika*zacetni
        cas += delta
        counter+=1
    if spremljaj:
        return spremljajoce
    if povprecje:
        return (cass,povprecjee)
    return cas        
#print(variraj([0.9,0.01,0.009,0.0000002,0.000000002]))
def analpovp(x):
    return 249*np.exp(-0.5*x)
def analpovp2(x):
    return 24*np.exp(-0.5*x)
def analvar(x):
    return -249*np.exp(-x)*(np.exp(-x)-1)
def analvar2(x):
    return -24*np.exp(-x)*(np.exp(-x)-1)    

def animate(t):
    global ax1
    global ax2
    plt.suptitle("N = {}".format(str(podatki[1][t])))
    ax1.clear()
    ax2.clear()
    ax1.plot(podatki[0][:(t+1)],podatki[1][:(t+1)],"k")
    ax1.set_xlabel("t")
    ax1.set_ylabel("N")
    kji = np.asarray(range(21))
    povp = podatki[1][t]*0.1
    ax2.set_title("Poissonova porazdelitev")
    ax2.plot(kji,povp**kji * np.exp(-povp)/(scipy.misc.factorial(kji)),"o--")    
    plt.tight_layout()
def drugaanimacija():
    frejmi = prehodna(250,0,1,0.0001,spremljaj=True)
    casi = [i*0.025 for i in range(len(frejmi))]
    print(len(casi))
    fig, ax = plt.subplots()
    def animiraj(t):
        ax.clear()
        plt.suptitle("t = {}".format(str(round(casi[t],2))))
        ax.plot(range(250),frejmi[t])
    ani = animation.FuncAnimation(fig,animiraj,range(0,len(casi)),interval=50)
    ani.save("prehodna.mp4")
def zajcilisice(z0,l0,alfa,beta,delta):
    zajci = [z0]
    lisice = [l0]
    cas = [0]
    while zajci[-1] > 0 and lisice[-1]>0:
        zajec = zajci[-1]
        lisica = lisice[-1]
        zajci.append(zajec + np.random.poisson(5*alfa*delta*zajec)-np.random.poisson(4*alfa*delta*zajec)-np.random.poisson(alfa/l0*delta*zajec*lisica))
        lisice.append(lisica + np.random.poisson(4*beta*delta*lisica)-np.random.poisson(5*beta*delta*lisica)+ np.random.poisson(beta/z0*delta*zajec*lisica))
        cas.append(cas[-1]+delta)
    return (zajci,lisice,cas)
if 0:
    fig, ax1 = plt.subplots()
    rezultat = zajcilisice(200,50,1,1,0.1)
    #ax1.plot(rezultat[2],rezultat[1])
    #ax1.plot(rezultat[2],rezultat[0])
    #ax1.set_xlabel("t")
    #ax1.set_title(r"$Z(t), L(t), \alpha=\beta=1, \Delta t = 0.001$")
    #plt.show()
    #plt.savefig("tretja/41.pdf")
    #ax1.clear()
    ax1.plot(rezultat[0],rezultat[1])
    ax1.set_xlabel("Z")
    ax1.set_ylabel("L")
    #ax1.set_title(r"$\alpha=\beta=1, \Delta t = 0.001$")
    #plt.savefig("tretja/42.pdf")
if 0:   
    #alfaz/alfal od 0.5 do 2
    cmap1 = plt.get_cmap("autumn")
    cmap2 = plt.get_cmap("summer")
    counter = 0
    for i in np.linspace(0.5, 2, 5):
        if counter==0:
            resitev = zajcilisice(200,50,i,1,0.001)
            plt.plot(resitev[2],resitev[0],label=r"Zajci, $\alpha / \beta = 0.5$",color = cmap1(counter))
            plt.plot(resitev[2],resitev[1],label=r"Lisice, $\alpha / \beta = 0.5$",color = cmap2(counter))
        resitev = zajcilisice(200,50,i,1,0.001)
        plt.plot(resitev[2],resitev[0],color = cmap1(counter))
        plt.plot(resitev[2],resitev[1],color = cmap2(counter/4))
        print(counter)
        if counter>0.34:
            resitev = zajcilisice(200,50,2,1,0.001)
            plt.plot(resitev[2],resitev[0],label=r"Zajci, $\alpha / \beta = 2$",color = cmap1(counter))
            plt.plot(resitev[2],resitev[1],label=r"Lisice, $\alpha / \beta = 2$",color = cmap2(counter))
        counter+=0.2
    plt.title(r"$Z(t),L(t), \Delta t = 0.001$")
    plt.xlabel("t")
    plt.legend(loc="best")
    plt.savefig("tretja/1.pdf")
    plt.show()
if 0:
    #alfaz/alfal od 0.5 do 2
    cmap = plt.get_cmap("hot")
    counter = 0
    for i in np.linspace(0.5, 2, 5):
        if counter==0:
            resitev = zajcilisice(200,50,i,1,0.001)
            plt.plot(resitev[0],resitev[1],label=r"$\alpha / \beta = 0.5 $",color = cmap(counter))
        resitev = zajcilisice(200,50,i,1,0.001)
        plt.plot(resitev[0],resitev[1],color = cmap1(counter))
        if counter>0.78:
            resitev = zajcilisice(200,50,2,1,0.001)
            plt.plot(resitev[0],resitev[1],label=r"$\alpha / \beta = 2$",color = cmap1(counter))
        counter+=0.2
    plt.title(r"$\Delta t = 0.001$")
    plt.xlabel("Z")
    plt.ylabel("L")
    plt.legend(loc="best")
    plt.savefig("tretja/2.pdf")
    plt.show()
if 0:
    #Končni za več delta t
    delte = np.linspace(0.0001,0.01,100)
    koncni = []
    for i in delte:
        casi = [zajcilisice(200,50,0.1,0.5,0.4,0.4,0.5,i) for i in range(10)]
        koncni.append(sum(casi)/10)
    plt.plot(delte,koncni)
    plt.title(r"$\gamma = 1, \Delta t = 0.001$")
    plt.xlabel(r"$\Delta t$")
    plt.ylabel(r"$t_{kon}$")
    plt.legend(loc="best")
    plt.savefig("tretja/delte.pdf")
    plt.show()    
#fig, (ax1,ax2) = plt.subplots(2)
#fig,ax1 = plt.subplots()
#ani = animation.FuncAnimation(fig,animate,range(0,len(podatki[0])),interval=250)    
#plt.show()
#ani.save("Poisson.mp4") 
if 1:
        #TO NUJNO POPRAVI ENKRAT
    resitev = prehodna(250,0,1,0.0001,True)
    plt.plot(resitev[0],resitev[1],label="N=250",color="b")
    x = np.linspace(0,resitev[0][-1],1000)
    plt.plot(x,analvar(x),label="N=250, analitično",color="k",ls="--")
    resitev = prehodna(25,0,1,0.0001,True)
    plt.plot(resitev[0],resitev[1],label="N=25",color="r")
    x = np.linspace(0,resitev[0][-1],1000)
    plt.plot(x,analvar2(x),label="N=25, analitično",color="m",ls="--")
    plt.legend(loc="best")
    plt.xlabel("t")
    plt.ylabel(r"$var$")
    plt.title(r"Odvisnost variance verjetnostne porazdelitve od časa $\beta_s=1, \beta_r = 0$")
    plt.savefig("druga/variance.pdf")
    plt.show()

if 0:
    cmap = plt.get_cmap("copper")
    for i in range(5):
        y = izumiraj(250,5,1,4)
        if i==3:
            plt.plot(y[0],y[1],color=cmap(i/4),label="N=250")
            continue
        plt.plot(y[0],y[1],color=cmap(i/4))
    cmap = plt.get_cmap("summer")
    for i in range(5):
        y = izumiraj(25,5,1,4)
        if i==3:
            plt.plot(y[0],y[1],color=cmap(i/4),label="N=25")
        plt.plot(y[0],y[1],color=cmap(i/4))
    plt.title(r"$\Delta t=0.01, \beta_r = 4, \beta_s = 5$")
    plt.xlabel("t")
    plt.ylabel("N")
    plt.legend(loc="best")
    #plt.savefig("prva/rojstva.pdf")
    plt.show()
if 0:
    povprecja = [povpinvar(casizumrtja,250,5,0.001,10,i)[0] for i in np.linspace(0.1,4,1000)]
    plt.plot(np.linspace(0.01,4,1000),povprecja,label="N=250")
    povprecja = [povpinvar(casizumrtja,25,5,0.001,10,i)[0] for i in np.linspace(0.1,4,1000)]
    plt.plot(np.linspace(0.01,4
                         ,1000),povprecja,label="N=25")
    plt.xlabel(r"$\beta_r$")
    plt.ylabel("t")
    plt.legend(loc="best")
    plt.title(r"Povprečen čas izumrtja v odvisnosti od $\beta_r$, $\Delta_t = 0.01, \beta_s = 5$")
    plt.savefig("prva/rojstvarojstva.pdf")
    plt.show()
if 0:
    color = plt.get_cmap("hsv")
    casi = [casizumrtja(250,5,0.01,4) for i in range(1000)]
    y = np.histogram(casi,bins="auto",normed=True)
    plt.plot([(y[1][i]+y[1][i+1])/2 for i in range(len(y[0]))],y[0],label=r"$N=250, \Delta t = 0.01$", color = color(0))
    casi = [casizumrtja(250,5,1,4) for i in range(1000)]
    y = np.histogram(casi,bins="auto",normed=True)
    plt.plot([(y[1][i]+y[1][i+1])/2 for i in range(len(y[0]))],y[0],label=r"$N=250, \Delta t = 1$", color = color(0.3))
    casi = [casizumrtja(25,5,0.01,4) for i in range(1000)]
    y = np.histogram(casi,bins="auto",normed=True)
    plt.plot([(y[1][i]+y[1][i+1])/2 for i in range(len(y[0]))],y[0],label=r"$N=25, \Delta t = 0.01$", color = color(0.6))
    casi = [casizumrtja(25,5,1,4) for i in range(1000)]
    y = np.histogram(casi,bins="auto",normed=True)
    plt.plot([(y[1][i]+y[1][i+1])/2 for i in range(len(y[0]))],y[0],label=r"$N=25, \Delta t = 1$", color = color(0.9))
    plt.xlabel("t")
    plt.legend(loc="best")
    plt.title("Porazdelitev časov izumrtja")
    plt.savefig("prva/rojstvaporazdelitev.pdf")
    plt.show()
    