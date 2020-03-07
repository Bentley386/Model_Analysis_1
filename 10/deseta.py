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
qwertz = [0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-8,-7,-6,-5,-4,-3,-2,-1,0]
def timer(code):
    start_time = timeit.default_timer()
    exec(code)
    elapsed = timeit.default_timer() - start_time
    return elapsed
def energijaverizice(alfa,hji):
    energija = 0
    s = np.asarray(hji)
    energija += np.sum(alfa*s)
    for i in range(0,18):
        energija+= 0.5*(s[i+1]-s[i])**2
    return energija
def energijaising(spini,h=0):
    energija = 0
    if h!=0:
        energija -= np.sum(spini*h)
    sosedi = np.zeros_like(spini)
    sosedi += np.roll(spini,1,0) +np.roll(spini,-1,0) + np.roll(spini,1,1) + np.roll(spini,-1,1)
    energija -= np.sum(spini*sosedi)
    return energija
def korak(hji,lokacija,predznak,nakljucen,alfa,T,energija):
    if hji[lokacija]+predznak<-19:
        return (hji,energija)
    deltaE = 1 - predznak*(hji[lokacija+1] - 2*hji[lokacija] + hji[lokacija-1] - alfa)
    boltzmann = np.exp(-deltaE/T)
    if nakljucen <= boltzmann:
        nov = hji
        nov[lokacija] = nov[lokacija] + predznak
        return (nov,energija+deltaE)
    else:
        return (hji,energija)          
def prva(alfa=1,T=1,N=1000000,zacetni=qwertz):
    hji0 = zacetni
    energy = [energijaverizice(alfa,hji0)]
    lokacije = np.random.randint(1,18,N)
    predznaki = np.random.choice((-1,1),N)
    nakljucna = np.random.rand(N)    
    for i in range(N):
        temp = korak(hji0,lokacije[i],predznaki[i],nakljucna[i],alfa,T,energy[i])
        hji0 = temp[0]
        energy.append(temp[1])
    return (temp[0],energy)

def korakising(spini,i,j,nakljucen,energija,T,d,h=0): #i je y j je x
    if h!=0:
        deltaE = 4*spini[i][j]*(spini[(i+1)%d][j]+spini[(i-1)%d][j]+spini[i][(j+1)%d] + spini[i][(j-1)%d]) + 2*h*spini[i][j]
    else:
        deltaE = 4*spini[i][j]*(spini[(i+1)%d][j]+spini[(i-1)%d][j]+spini[i][(j+1)%d] + spini[i][(j-1)%d])
    if nakljucen <= np.exp(-deltaE/T):
        np.put(spini,(i*d+j),spini[i][j]*(-1))
        return energija + deltaE
    else:
        return energija
def razlika(spini,grozd,d):
    meje = (np.amin(grozd[:,0]),np.amax(grozd[:,0]),np.amin(grozd[:,1]),np.amax(grozd[:,1]))
    nova = []
    delta = 0
    for i in range(len(spini)):
        vrstica = []
        for j in range(len(spini)):
            if (i,j) not in grozd:
                vrstica.append(spini[i][j])
            else:
                vrstica.append(-spini[i][j])
        nova.append(vrstica)
    for i in (range((meje[0]-1)%d,(meje[1]+2)%d)):
        for j in range((meje[2]-1)%d,(meje[3]+2)%d):
            delta += -(nova[i][j]*nova[(i+1)%d][j] + nova[i][j]*nova[i][(j+1)%d] - spini[i][j]*spini[(i+1)%d][j]-spini[i][j]*spini[i][(j+1)%d])
    return (nova,delta)
def wolff(spini,i,j,nakljucen,energija,T,d,h=0):
    grozd = [(i,j)]
    def rekurzija(k,l):
        if spini[(k+1)%d][l%d] == spini[k%d][l%d] and ((k+1)%d,l%d) not in grozd:
            if np.random.rand() <= np.exp(-2/T):
                return
            grozd.append(((k+1)%d,l%d))
            rekurzija((k+1)%d,l%d)
        if spini[(k-1)%d][l%d] == spini[k%d][l%d] and ((k-1)%d,l%d) not in grozd:
            if np.random.rand() <= np.exp(-2/T):
                return
            grozd.append(((k-1)%d,l%d))
            rekurzija((k-1)%d,l%d)
        if spini[k%d][(l+1)%d] == spini[k%d][l%d] and (k%d,(l+1)%d) not in grozd:
            if np.random.rand() <= np.exp(-2/T):
                return
            grozd.append((k%d,(l+1)%d))
            rekurzija(k%d,(l+1)%d)
        if spini[k%d][(l-1)%d] == spini[k%d][l%d] and (k%d,(l-1)%d) not in grozd:
            if np.random.rand() <= np.exp(-2/T):
                return
            grozd.append((k%d,(l-1)%d))
            rekurzija(k%d,(l-1)%d)
        return
    rekurzija(i,j)
    return razlika(spini,np.asarray(grozd),d)
    
def drugawolff(zacetni,T=1,h=0,N=100,n=100,cas=False):
    if cas:
        startni = timeit.default_timer()
        casi = [startni]
    spini = zacetni
    energy = [energijaising(spini,h)]
    lokacije = np.random.randint(0,N,(n,2))
    nakljucna = np.random.rand(n)            
    for i in range(n):
        rezultat=wolff(spini,lokacije[i][0],lokacije[i][1],nakljucna[i],energy[-1],T,N,h)
        energy.append(energy[-1]+rezultat[1])
        spini = rezultat[0]
        if cas:
            casi.append(timeit.default_timer())
    if cas:
        return (spini,energy,np.asarray(casi)-startni)
    return (spini,energy)
    
def druga(zacetni,T=1,h=0,N=100,n=100,cas=False):
    if cas:
        startni = timeit.default_timer()
        casi = [startni]
    spini = zacetni
    energy = [energijaising(spini,h)]
    lokacije = np.random.randint(0,N,(n,2))
    nakljucna = np.random.rand(n)            
    for i in range(n):
        energy.append(korakising(spini,lokacije[i][0],lokacije[i][1],nakljucna[i],energy[-1],T,N,h))
        if cas:
            casi.append(timeit.default_timer())
    if cas:
        return (spini,energy,np.asarray(casi)-startni)
    return (np.asarray(spini),energy)
def razdalja(lokacije):
    d = 0
    for i in range(9):
        d+= (lokacije[i+1][0]-lokacije[i][0])**2 + (lokacije[i+1][1]-lokacije[i][1])**2
    return d**0.5
def salesman(T,zacetna):
    lokacije = [zacetna]
    for i in range(5000):
        erm = korak2(lokacije[-1],T)
        lokacije.append(erm[0])
        #print(dolzina
    return (lokacije,erm[1])    
def korak2(lokacije,T):
    dolzina1 = razdalja(lokacije)
    poteza = np.random.randint(0,10,2)
    nova = []
    for i in range(10):
        if i==poteza[0]:
            nova.append(lokacije[poteza[1]])
        elif i==poteza[1]:
            nova.append(lokacije[poteza[0]])
        else:
            nova.append(lokacije[i])
    #np.put(nova,poteza[0],lokacije[poteza[1]])
    #np.put(nova,poteza[1],lokacije[poteza[0]])
    nova = np.asarray(nova)
    dolzina2 = razdalja(nova)
    nakljucen = np.random.rand()
    #print(lokacije)
    #print("er")
    #print(nova)
    if nakljucen <= np.exp(-(dolzina2-dolzina1)/T):

        return [nova,dolzina2]  
    else:
        return [lokacije,dolzina1]

zacetna = np.random.rand(10,2)*10
temperature = np.linspace(1,0.05,50)
#print(zacetna)
#korak2(zacetna,1)
def animiraj(t):
    global ax1
    global zacetna
    temp = salesman(temperature[t],zacetna)
    zacetna = temp[0][-1][:]
    ax1.clear()
    ax1.plot(zacetna[:,0],zacetna[:,1],"o-")
    plt.suptitle("T={}, d={}".format(str(round(temperature[t],2)),str(round(temp[1],2))))
    
fig, ax1 = plt.subplots()
#fig,ax1 = plt.subplots()
ani = animation.FuncAnimation(fig,animiraj,range(0,50),interval=250)    
#plt.show()
ani.save("merchant.mp4") 
if 0:
    #plt.imshow(energije,cmap="hot",interpolation="nearest")
    fig, ax1 = plt.subplots()
    colormap = plt.get_cmap("copper")
    plt.title("T=1,h=0, št. potez = 5M")
    ax1.set_xlabel("N")
    ax1.set_ylabel("E")
    idk,energije,casi = druga(np.random.choice((1,-1),size=(100,100)),n=5000000,N=100,cas=True)
    ax1.plot(range(0,5000001),energije,color=colormap(0.1),label="N=100")  
    energije = druga(np.random.choice((1,-1),size=(200,200)),n=5000000,N=200)[1]
    ax1.plot(range(0,5000001),energije,color=colormap(0.4),label="N=200")  
    energije = druga(np.random.choice((1,-1),size=(300,300)),n=5000000,N=300)[1]
    ax1.plot(range(0,5000001),energije,color=colormap(0.8),label="N=300")
    ax2 = ax1.twinx()
    ax2.plot(range(0,5000001),casi,"--",color="r")  
    ax2.set_ylabel("t")
    ax1.legend(loc="best")
    plt.savefig("druga/wolff.pdf")
    plt.show()        
if 0:
    hji = [0,0.1,1,10]
    barve = np.linspace(0.3,0.9,4)
    colormap = plt.get_cmap("BuGn")
    for i in range(4):
        energije = []
        zacetna = np.random.choice((1,-1),size=(100,100))
        for temperature in np.linspace(5,1,20):
            temp = druga(zacetna,n=1000000,N=100,T=temperature,h=hji[i])
            zacetna = temp[0]
            temp = np.asarray(temp[1][900000::10000])
            povp = np.sum(temp)/temp.size
            povp2 = np.sum(temp**2)/temp.size
            print(povp2-povp**2)
            energije.append((povp2-povp**2)/(10000*temperature**2))
        plt.plot(np.linspace(5,0.1,20),energije,label="h={}".format(str(hji[i])),color=colormap(barve[i]))
    
    #fig, ax1 = plt.subplots()
    #colormap = plt.get_cmap("copper")
    #plt.title("T=1,h=0, št. potez = 5M")
    #ax1.set_xlabel("N")
    #ax1.set_ylabel("E")
    #idk,energije,casi = druga(np.random.choice((1,-1),size=(100,100)),n=1000000,N=100,cas=True)
    #plt.imshow(idk,cmap="hot",interpolation="nearest")
    #plt.plot(range(0,1000001),energije,color=colormap(0.1),label="N=100")  
    #energije = druga(np.random.choice((1,-1),size=(200,200)),n=5000000,N=200)[1]
    #ax1.plot(range(0,5000001),energije,color=colormap(0.4),label="N=200")  
    #energije = druga(np.random.choice((1,-1),size=(300,300)),n=5000000,N=300)[1]
    #ax1.plot(range(0,5000001),energije,color=colormap(0.8),label="N=300")
    #ax2 = ax1.twinx()
    #ax2.plot(range(0,5000001),casi,"--",color="r")  
    #ax2.set_ylabel("t")
    plt.xlabel("T")
    plt.ylabel(r"C")
    plt.title(r"$<C>(T)$")
    plt.axvline(x=2.27,color="r",ls="--")
    plt.legend(loc="best")
    plt.savefig("druga/kapac2.pdf")
    plt.show()      
    
    
    
    
    
    
if 0:
    plt.title(r"Odvisnost energije od T")
    plt.xlabel("T")
    plt.ylabel("E")
    colormap = plt.get_cmap("cool")
    alfe = [0.1,0.5,1,2,10]
    barve = np.linspace(0.05,0.95,len(alfe))
    for t in range(len(alfe)):
        temperature = np.linspace(0.01,100)
        enrgije = []
        for i in temperature:
            energije = prva(alfa=alfe[t],T = i,N=100000)[1]
            energije=np.asarray(energije[:50000:1000])
            enrgije.append(np.sum(energije)/energije.size)
        plt.plot(temperature,enrgije,label=r"$\alpha={}$".format(str(alfe[t])),color=colormap(barve[t]))
    plt.legend(loc="best")
    plt.savefig("prva/EodT.pdf")
    plt.show()    
if 0:
    plt.title(r"Odvisnost energije od N, $\alpha = 1$")
    plt.xlabel("N")
    plt.ylabel("E")
    colormap = plt.get_cmap("summer")
    tt = [10,100,50]
    barve = np.linspace(0,1,len(tt))
    for t in range(len(tt)):
        energije = prva(alfa=10,T=tt[t],N=1000000)[1]
        energije=energije[::10]
        x = np.linspace(0,1000000,len(energije))
        plt.plot(x,energije,label="T={}".format(str(tt[t])),color=colormap(barve[t]))
    plt.legend(loc="best")
    #plt.savefig("prva/energijaN2.png")
    plt.show()
   
    
    
    
    
"""    
temperature = np.linspace(50,0.1,10)
trenutni = qwertz
razporedi = []
for i in range(10):
    razpored = prva(T=temperature[i],N=100000,zacetni=trenutni)
    trenutni = razpored[0]
    razporedi.append((trenutni[:],razpored[1][-1]))
def animate2(t):
    global axi
    ax1.clear()
    ax1.plot(range(0,19),razporedi[t][0],"o--")
    ax1.set_xlim([0,18])
    ax1.set_ylim([-19,0])
    ax1.set_title(r"$\alpha = 1$")
    ax1.text(7,-1,"T={}".format(str(round(temperature[t],2))))
    ax1.text(7,-3,"E={}".format(str(round(razporedi[t][1],2))))    
"""    
"""
temperature = np.linspace(10,0.1,20)

razporedi = []
hji = [0,0.1,1,10]
for i in range(20):
        razpored = druga(trenutni,T=temperature[i],N=100,n=100000,h=0)[0]
        trenutni = razpored[:]
        razporedi.append(trenutni[:])
"""        
#plt.imshow(razporedi[0][5][0],cmap="hot",interpolation="nearest")
#fig, (ax1,ax2) = plt.subplots(2)
#razpored = druga(trenutni,T=5,N=100,n=100000,h=0)
#ax1.imshow(razporedi[0],cmap="hot",interpolation="nearest")
#trenutni = razpored[0][:]
#razpored = druga(trenutni,T=1,N=100,n=100000,h=0)
#ax2.imshow(razporedi[5],cmap="hot",interpolation="nearest")

#temperature = np.linspace(10,0.1,30)
#trenutni = [np.random.choice((1,-1),size=(100,100)),np.random.choice((1,-1),size=(100,100)),np.random.choice((1,-1),size=(100,100)),np.random.choice((1,-1),size=(100,100))]
def animate(t):
    global axi
    global trenutni
    #plt.title("test")
    razpored = druga(trenutni[0],T=temperature[t],N=100,n=100000,h=0)
    trenutni[0] = razpored[0][:]
    axi[0][0].clear()
    axi[0][0].imshow(razpored[0][:],cmap="hot")
    #axi[0][0].set_ylim([-19,0])
    axi[0][0].set_title(r"$h=0$")
    #axi[0][0].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[0][0].text(-23,-2,"E={}".format(str(round(razpored[1][-1],2))))
    axi[0][1].clear()
    razpored = druga(trenutni[1],T=temperature[t],N=100,n=100000,h=0.1)
    trenutni[1] = razpored[0][:]
    axi[0][1].imshow(trenutni[1],cmap="hot")
    axi[0][1].set_title(r"$h=0.1$")
    #axi[0][1].set_ylim([-19,0])
    #axi[0][1].text(-23,-2,"T={}".format(str(round(temperature[t],2))))
    axi[0][1].text(-23,-2,"E={}".format(str(round(razpored[1][-1],2))))
    axi[1][0].clear()
    razpored = druga(trenutni[2],T=temperature[t],N=100,n=100000,h=1)
    trenutni[2] = razpored[0][:]
    axi[1][0].imshow(trenutni[2],cmap="hot")
    #axi[1][0].set_ylim([-19,0])
    #axi[1][0].set_xlim([0,18])
    axi[1][0].set_title(r"$h=1$")
    #axi[1][0].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[1][0].text(-23,-2,"E={}".format(str(round(razpored[1][-1],2))))
    axi[1][1].clear()
    razpored = druga(trenutni[3],T=temperature[t],N=100,n=100000,h=10)
    trenutni[3] = razpored[0][:]
    axi[1][1].imshow(trenutni[3],cmap="hot")
    #axi[1][1].set_xlim([0,18])
    axi[1][1].set_title(r"$h=10$")
    #axi[1][1].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[1][1].text(-23,-2,"E={}".format(str(round(razpored[1][-1],2))))
    plt.suptitle("T={}".format(str(round(temperature[t],2))))
    plt.tight_layout()
    #ax1.text(7,-1,"T={}".format(str(round(temperature[t],2))))
    #ax1.text(7,-3,"E={}".format(str(round(razporedi[t][1],2)))) 
"""   
fig, axi = plt.subplots(2, 2, sharex='col', sharey='row')
#fig,ax1 = plt.subplots()
ani = animation.FuncAnimation(fig,animate,range(0,30),interval=300)    
#plt.show()
ani.save("ising.mp4") 
"""    
"""temperature = np.linspace(10,0.1,100)
trenutni = qwertz
razporedi = [[],[],[],[]]
alfe = [0.1,0.5,1,2]
for j in range(4):
    for i in range(100):
        razpored = prva(T=temperature[i],N=100000,zacetni=trenutni,alfa=alfe[j])
        trenutni = razpored[0]
        razporedi[j].append((trenutni[:],razpored[1][-1]))

def animate(t):
    global axi
    #plt.title("test")
    axi[0][0].clear()
    axi[0][0].plot(range(0,19),razporedi[0][t][0],"o--",color="blue")
    axi[0][0].set_ylim([-19,0])
    axi[0][0].set_title(r"$\alpha=0.1$")
    #axi[0][0].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[0][0].text(7,-3,"E={}".format(str(round(razporedi[0][t][1],2))))
    axi[0][1].clear()
    axi[0][1].plot(range(0,19),razporedi[1][t][0],"o--",color="red")
    axi[0][1].set_title(r"$\alpha=0.5$")
    axi[0][1].set_ylim([-19,0])
    #axi[0][1].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[0][1].text(7,-3,"E={}".format(str(round(razporedi[1][t][1],2))))
    axi[1][0].clear()
    axi[1][0].plot(range(0,19),razporedi[2][t][0],"o--",color="magenta")
    axi[1][0].set_ylim([-19,0])
    axi[1][0].set_xlim([0,18])
    axi[1][0].set_title(r"$\alpha=1$")
    #axi[1][0].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[1][0].text(7,-3,"E={}".format(str(round(razporedi[2][t][1],2))))
    axi[1][1].clear()
    axi[1][1].plot(range(0,19),razporedi[3][t][0],"o--",color="orange")
    axi[1][1].set_xlim([0,18])
    axi[1][1].set_title(r"$\alpha=2$")
    #axi[1][1].text(7,-1,"T={}".format(str(round(temperature[t],2))))
    axi[1][1].text(7,-3,"E={}".format(str(round(razporedi[3][t][1],2))))
    plt.suptitle("T={}".format(str(round(temperature[t],2))))
    plt.tight_layout()
    #ax1.text(7,-1,"T={}".format(str(round(temperature[t],2))))
    #ax1.text(7,-3,"E={}".format(str(round(razporedi[t][1],2)))) 
"""    
    
