# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 21:35:07 2018

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 10:46:11 2017

@author: Admin
"""
import timeit
from PIL import Image
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
from matplotlib.patches import Ellipse
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rc("text",usetex=True)
matplotlib.rcParams["text.latex.unicode"] = True
plt.close("all")
pi = np.pi
#signal = np.loadtxt("val3.dat")

def gama(x,k,theta):
    return x**(k-1)*np.exp(-x/theta) /(theta**k * spec.gamma(k))
#x = np.linspace(0,20,1000)
#plt.plot(x,gama(x,2,2))
def prenosna(x,tau=30,N=512):
    return 1/tau * np.exp(-x/tau)
    if x <=N/2:
        return 1/(tau)*np.exp(-np.abs(x)/tau)
    else:
        return 1/(tau)*np.exp(-np.abs(x-N)/tau)
def prenosna2(x,tau=16):
    return 1/(1+(tau*x)**2)
def okno(ime,x):
    N = x.size
    def turki(i):
        if i<(N-1)/4:
            return 0.5*(1+np.cos(pi*(4*i/(N-1) - 1)))
        elif i >= (N-1)/4 and i< (N-1)*(0.75):
            return 1
        else:
            return 0.5*(1+np.cos(pi*(4*i/(N-1) -3)))
    if ime == "tri":
        return 1 - np.absolute((x-0.5*(N-1))/(0.5*(N-1)))
    elif ime == "welch":
        return 1 - ((x-0.5*(N-1))/(0.5*(N-1)))**2
    elif ime == "hann":
        return np.sin(pi*x/(N-1))**2
    elif ime == "srs":
        return 1 - 1.93*np.cos((2*pi*x)/(N-1)) + 1.29*np.cos((4*pi*x)/(N-1)) - 0.388*np.cos((6*pi*x)/(N-1)) + 0.028* np.cos((8*pi*x)/(N-1))
    elif ime == "turkey":
        return np.vectorize(turki)(x)
    elif ime=="poisson":
        return np.exp(-np.absolute(x-0.5*(N-1))/(0.5*N*8.69/60))

    raise NameError("Napačno okno")

def spekterstari(x, ime=0):
    N = x.size
    if ime==0:
        fourier = np.fft.fft(x)
        moc = []
        for i in range(int(N/2 + 1)):
            if i==0:
                moc.append(np.absolute(fourier[0])**2)
                continue
            elif i==N/2:
                moc.append(np.absolute(fourier[i])**2)
                continue
            moc.append(0.5*(np.absolute(fourier[i])**2+np.absolute(fourier[N-i])**2))
        return [range(int(N/2+1)),moc]
    oknce = np.asarray(range(N))
    signal = x*okno(ime,oknce)
    fourier = np.fft.fft(signal)
    moc = []
    for i in range(int(N/2 + 1)):
        if i==0:
            moc.append(np.absolute(fourier[0])**2)
            continue
        elif i==N/2:
            moc.append(np.absolute(fourier[i])**2)
            continue
        moc.append(0.5*(np.absolute(fourier[i])**2+np.absolute(fourier[N-i])**2))
    return [range(int(N/2+1)),moc]


def gaussovka(x,mu,sigma):
    return 1/np.sqrt(2*pi*sigma**2) * np.exp(-(x-mu)**2/(2*sigma**2))

def avtokorelacija(x,i):
    N = x.size
    y = np.concatenate((x,np.zeros(N)))[i:i+N]
    return 1/(N-i)*np.sum(x*y)        
def aji(x,p):
    rji = np.asarray([avtokorelacija(x,i) for i in range(p)])
    matrika = scipy.linalg.toeplitz(rji,rji)
    return solve(matrika,-np.append(rji[1:],avtokorelacija(x,p)))
def spekter2(a,w):
    arej = np.exp(-1j * 2*pi*w/100 * np.asarray(range(1,a.size+1)))
    return np.abs(1+np.sum(a*arej))**(-2)
def avtokonvolucija(x,i):
    N = x.size
    y = np.concatenate((np.zeros(N),x))[N+1-i:N+1+N-i]
    return 1/(N)*np.sum(x*y)    
def napaka(a,x):
    N = x.size
    E = np.sum(x**2)/N
    eji = np.asarray([avtokonvolucija(x,i) for i in range(1,a.size+1)])
    return E + np.sum(a*eji)
def napovej(x,a):
    signal = x
    N = x.size
    p = a.size
    pass
    for i in range(N):
        novi = - np.sum(a*signal[signal.size:signal.size-p-1:-1])
        signal = np.append(signal,novi)
    return signal
def preslikaj(x):
    novenicle = []
    for i in range(x.size):
        if np.abs(x[i]) > 1:
            x[i] = x[i]/np.abs(x[i])
            novenicle.append(x[i])
    return np.array(novenicle)
def napaka(data, p, normalizacija = True):
    coef = getCoef(data, p)
    a = np.correlate(data, data, mode = 'valid')
    if normalizacija:
        a /= data.size
        a += np.sum(coef * np.correlate(np.concatenate((data,np.zeros(p))), data, mode = "valid")[1:] / (data.size-np.arange(1, p+1)))
    else:
        a += np.sum(coef * np.correlate(np.concatenate((data,np.zeros(p))), data, mode = "valid")[1:])
    return a
def zasumi(x,mu,sigma):
    if sigma==0:
        return x
    gaussi = np.random.normal(mu,sigma,x.size)
    return x+gaussi
if 1:
    #ODSTOPANJE OD PRAVE VREDNOSTI
    podatki = np.loadtxt("luna.dat",usecols=1)
    print(podatki.shape)
    #print(podatki.size)
    #negativne = np.where(podatki[:,1]<0)
    #podatki = np.delete(podatki,negativne,0)
    #koef = np.polyfit(podatki[:,0],podatki[:,1],1)
    #podatki[:,1] = podatki[:,1] - koef[0] * podatki[:,0] - koef[1]
    #podatki = podatki[:,1]
    podatki2 = zasumi(podatki,5,0.1)
    a = aji(podatki2,27)
    nicle = np.roots(np.append(a[::-1],1)[::-1])            
    preslikaj(nicle)
    a = np.poly(nicle)[1:]
    polovica = podatki2[:int(podatki2.size/2)]
    napovedan = napovej(polovica,a)    
    if 0:
        sigme = [0,0.01,0.05,0.1,0.15,0.2,0.5,1]
        for j in sigme:
            odstopanja = []
            for k in range(100):
                odstopanjatemp = []
                podatki2 = zasumi(podatki,0,j)
                pji = range(2,21)
                for i in pji:
                    a = aji(podatki2,i)
                    nicle = np.roots(np.append(a[::-1],1)[::-1])            
                    preslikaj(nicle)
                    a = np.poly(nicle)[1:]
                    polovica = podatki2[:int(podatki2.size/2)]
                    napovedan = napovej(polovica,a)
                    odstopanje= np.sum(np.abs(napovedan[int(podatki2.size/2):] - podatki2[int(podatki2.size/2):]))/polovica.size
                    odstopanjatemp.append(odstopanje)
                if k==0:
                    odstopanja = np.array(odstopanjatemp)
                else:
                    odstopanja = odstopanja + np.array(odstopanjatemp)
            
            plt.plot(pji,odstopanja/100,"o-",label=r"$\sigma = {0}$".format(str(j)))
    plt.plot(range(podatki2.size),podatki2,label="Pravi")
    plt.plot(range(podatki2.size),napovedan,label="Napovedan")
    plt.title(r"Linearna napoved za luno - RA, zašumljeno z $\mu=5, \sigma=0.1$")
    plt.ylabel("RA [h]")
    plt.xlabel("Cas [dan]")
    #plt.yscale("log")
    plt.legend(loc="best")
    plt.savefig("druga/sumi9.pdf")
#podatki = np.loadtxt("luna.dat",usecols=2)
#plt.xlabel("Čas [dan]")
#plt.ylabel("Dec[stopinje]")
#plt.title("Lunarne efemeride, 1995-2000")
#plt.plot(range(podatki.size),podatki)
#plt.savefig("druga/lunaraw2.pdf")
if 0:
    podatki = np.loadtxt("luna.dat",usecols=2)
    #negativne = np.where(podatki[:,1]<0)
    #podatki = np.delete(podatki,negativne,0)
    #koef = np.polyfit(podatki[:,0],podatki[:,1],1)
    #podatki[:,1] = podatki[:,1] - koef[0] * podatki[:,0] - koef[1]
    #podatki2=podatki[:,0]
    #podatki = podatki[:,1]
    a = aji(podatki,25)
    #LINEARNA NAPOVED
    if 1:
        #PREVERI POLE
        x = np.linspace(-1,1,1000)
        nicle = np.roots(np.append(a[::-1],1)[::-1])            
        #plt.plot(np.real(nicle),np.imag(nicle),"mo",label="p=15")
        preslikaj(nicle)
        #plt.plot(np.real(nicle),np.imag(nicle),"ro")        
        #plt.plot(x,np.sqrt(1-x**2),"b")
        #plt.plot(x,-np.sqrt(1-x**2),"b")
        #plt.savefig("druga/val33poli.pdf")
        a = np.poly(nicle)[1:]
    #raise NameError("Zaenkrat samo preverjamo pole")
    print(podatki.size)
    x = range(podatki.size)
    polovica = podatki[:int(podatki.size/2)]
    napovedan = napovej(polovica,a)
    plt.xlabel("Čas [dan]")
    plt.ylabel("Dec [stopinje]")
    plt.plot(range(podatki.size),podatki,label="Pravi")
    plt.plot(range(podatki.size),napovedan,label="Napovedan")
    plt.title("Linearna napoved za luno, p=25")
    plt.legend(loc="best")
    plt.savefig("druga/luna2.pdf")
if 0:
    #LOCLJIVOST NE POZABI SPREMENIT NJA V SPEKTR2
    x = np.linspace(0,100,10000)
    y = np.sin(2*x) + np.sin(2.2*x) + np.sin(2.4*x)
    #omege = np.linspace(0,300/2/pi,10000)
    #a = aji(y,10)
    #spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    spektr = spekterstari(y,"hann")
    plt.yscale("log")
    plt.title(r"MEM - Spekter $sin(2t) + sin(2.1 t) + sin(2.2 t)$")
    plt.ylabel(r"$|P(\omega)|^2$")
    plt.xlabel(r"$\omega$")
    plt.plot(np.array(spektr[0])*2*pi/10000,spektr[1])
    #a = aji(y,15)
    #spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    #plt.plot(omege*2*pi/100,spekter,label="p=15")
    #a = aji(y,20)
    #spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    #plt.plot(omege*2*pi/100,spekter,label="p=20")
    plt.xlim(0,3)
    plt.legend(loc="best")
    #plt.savefig("prva/locljivost5.pdf")

if 0:
    #POLI mogli bi bit pri k=46.5, k=93 in k=139
    podatki = np.loadtxt("co2.dat")
    negativne = np.where(podatki[:,1]<0)
    podatki = np.delete(podatki,negativne,0)
    koef = np.polyfit(podatki[:,0],podatki[:,1],1)
    podatki[:,1] = podatki[:,1] - koef[0] * podatki[:,0] - koef[1]
    a = aji(podatki[:,1],6)
    x = np.linspace(-1,1,1000)
    nicle = np.roots(np.append(a[::-1],1)[::-1])
    plt.plot(x,np.sqrt(1-x**2),"b")
    plt.plot(x,-np.sqrt(1-x**2),"b")
    plt.plot(np.real(nicle),np.imag(nicle),"go",label="p=6")
    a = aji(podatki[:,1],10)
    nicle = np.roots(np.append(a[::-1],1)[::-1])
    #plt.plot(np.real(nicle[1:]),np.imag(nicle[1:]),"ro",label="p=10")
    a = aji(podatki[:,1],4)
    nicle = np.roots(np.append(a[::-1],1)[::-1])
    plt.plot(np.real(nicle),np.imag(nicle),"bo",label="p=4")
    a = aji(podatki[:,1],15)
    nicle = np.roots(np.append(a[::-1],1)[::-1])
    plt.plot(np.real(nicle),np.imag(nicle),"mo",label="p=15")
    N = podatki[:,1].size
    kji = np.asarray([46.7,93,139])
    nicle = np.exp(-2j*pi*kji/N)
    plt.plot(np.real(nicle),np.imag(nicle),"ro",label="Vrhi")
    plt.axes().set_aspect("equal")
    plt.legend(loc="upper left")
    plt.savefig("prva/poli.pdf")
if 0:
    #SPEKTER VAL2 VAL3 Z MEM
    #podatki = np.loadtxt("val2.dat")
    if 1:
        podatki = np.loadtxt("co2.dat")
        negativne = np.where(podatki[:,1]<0)
        podatki = np.delete(podatki,negativne,0)
        print(podatki[:,0].size)
        koef = np.polyfit(podatki[:,0],podatki[:,1],1)
        podatki[:,1] = podatki[:,1] - koef[0] * podatki[:,0] - koef[1]
    #pji = range(1,21)
    #ajcki = [aji(podatki,p) for p in pji]
    #napakce = [napaka(ajcek,podatki) for ajcek in ajcki]
    #plt.plot(pji,napakce)
    omege = np.linspace(0,303,10000)
    a = aji(podatki[:,1],15)
    spektrstari = spekterstari(podatki[:,1],"hann")
    plt.yscale("log")
    plt.plot(spektrstari[0],spektrstari[1],"--",label="FFT")
    spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    #plt.plot(omege,spekter,label="p=3")
    #a = aji(podatki[:,1],5)
    omege = np.linspace(0,303,10000)
    spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    plt.plot(omege,spekter,label="N=10000")
    #a = aji(podatki[:,1],10)
    omege = np.linspace(0,303,1000)
    spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    plt.plot(omege,spekter,label="N=1000")
    #a = aji(podatki[:,1],15)
    omege = np.linspace(0,303,100)
    spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    plt.plot(omege,spekter,label="N=100")
    #a = aji(podatki[:,1],20)
    omege = np.linspace(0,303,50)
    spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    plt.plot(omege,spekter,label="N=50")
    plt.title(r"Spekter $CO_2$, p=15")
    plt.legend(loc="best")
    plt.savefig("prva/co2vecN.pdf")
if 0:
    t = np.linspace(0,1,512)
    y = np.sin(50*t)
    a = aji(y,10)
    omege = np.linspace(0,256,100)
    spekter = np.vectorize(spekter2,excluded={0})(a,omege)
    plt.plot(omege,spekter)    
    

if 0:
    #ALISAING
    #plt.title("val2, Hannovo okno")
    plt.ylabel("P")
    plt.xlabel("f")
    plt.yscale("log")
    signal = np.loadtxt("val3.dat")
    haha = spekter(signal)
    dodatna = haha[1][128:][::-1]
    plt.plot(haha[0],haha[1],label="N=512")
    haha = spekter(signal[::2])
    plt.plot(np.asarray(haha[0]),haha[1],label="N=256")  
    plt.plot(haha[0],dodatna,"--",label="Preslikana N=512")
    plt.title("val3, Potujitev")
    #haha = spekter(signal[:128])
    #plt.plot(np.asarray(haha[0])*4,haha[1],label="N=128")
    #haha = spekter(signal[:64])
    #plt.plot(np.asarray(haha[0])*8,haha[1],label="N=64")
    plt.legend(loc="best")
    plt.savefig("prva/alias.png")    
if 0:
    #VEC N
    #plt.title("val2, Hannovo okno")
    plt.ylabel("P")
    plt.xlabel("f")
    plt.yscale("log")
    signal = np.loadtxt("val3.dat")
    haha = spekter(signal)
    plt.plot(haha[0],haha[1],label="N=512")
    haha = spekter(signal[:256])
    plt.plot(np.asarray(haha[0]),haha[1],label="N=256")  
    #haha = spekter(signal[:128])
    #plt.plot(np.asarray(haha[0])*4,haha[1],label="N=128")
    #haha = spekter(signal[:64])
    #plt.plot(np.asarray(haha[0])*8,haha[1],label="N=64")
    plt.legend(loc="best")
    plt.savefig("prva/zdravobresar.png")
if 0:
    #OKNJENI SPEKTRI
    signal = np.loadtxt("val2.dat")
    plt.title("Spekter val2, N=512")
    plt.xlabel("f")
    plt.ylabel("P")
    plt.yscale("log")
    #plt.xlim(45,55)
    #plt.ylim(1000,1000000)
    haha = spekter(signal)
    plt.plot(haha[0],haha[1],label="Brez okna")
    #haha = spekter(signal,"tri")
    #plt.plot(haha[0],haha[1],label="Trikotno")
    #haha = spekter(signal,"welch")
    #plt.plot(haha[0],haha[1],label="Welch")
    #haha = spekter(signal,"hann")
    #plt.plot(haha[0],haha[1],label="Hann")
    #haha = spekter(signal,"turkey")
    #plt.plot(haha[0],haha[1],label="Turkey")
    haha = spekter(signal,"poisson")
    plt.plot(haha[0],haha[1],label="Poisson")
    #haha = spekter(signal,"srs")
    #plt.plot(haha[0],haha[1],color="m",label="SRS")
    plt.legend(loc="best")
    plt.savefig("prva/val2okna5.pdf")
if 0:
    #OKNA
    fig, ax = plt.subplots()
    x = np.asarray(range(512))
    ax.plot(x,okno("tri",x),label="Trikotno")
    ax.plot(x,okno("welch",x),label="Welch")
    ax.plot(x,okno("hann",x),label="Hann")
    ax2 = ax.twinx()
    ax2.plot(x,okno("srs",x),"--",label="SRS",color="m")
    ax.plot(x,okno("turkey",x),label="Turkey")
    ax.plot(x,okno("poisson",x),label="Poisson")
    ax.set_xlim(0,512)
    ax.legend(loc="best")
    ax2.legend(loc="best")
    plt.savefig("prva/okna2.pdf")
if 0:
    #PLOTTANJE SIGNALOV
    signal = np.loadtxt("val3.dat")
    povp = np.sum(signal)/signal.size
    print(povp)
    signal = signal - povp
    signal = [np.linspace(0,1,signal.size),signal]
    plt.ylabel("s(t)")
    plt.xlabel("t")
    plt.title("Val3")
    plt.plot(signal[0],signal[1])
    plt.savefig("prva/val2test.pdf")
























