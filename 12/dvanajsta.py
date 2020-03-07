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
def shitwiener(x,N):
    def erm(ermm):
        if ermm<0:
            return 0
        else:
            return ermm
    spektr = np.abs(np.fft.fft(x))**2
    sumcek = np.ones(x.size)*np.sum(spektr[100:412])/312
    signal = spektr - sumcek
    fi =signal/(spektr)
    return np.fft.ifft((np.fft.fft(x)/np.fft.fft(np.vectorize(prenosna)(np.asarray(range(512))))*fi))
def wiener(x,N):
    def erm(ermm):
        if ermm<0:
            return 0
        else:
            return ermm
    spektr = np.abs(np.fft.fft(x))**2
    sumcek = np.ones(x.size)*np.sum(spektr[100:412])/312
    signal = spektr - sumcek
    signal = np.vectorize(erm)(signal)
    aux = [1 for i in range(N)] + [0 for i in range(512-2*N)] + [1 for i in range(N)]
    fi =signal/(spektr)*np.asarray(aux)
    return np.fft.ifft((np.fft.fft(x)/np.fft.fft(np.vectorize(prenosna)(np.asarray(range(512))))*fi))
def spekter(x, ime=0):
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
def preberisliko(strink):
    slika = []
    slikca = []
    with open(strink,"r") as f:
        f.readline()
        (sirina, visina) = f.readline().split()
        sirina = int(sirina)
        visina = int(visina)
        f.readline()
        for line in f:
            for stevilo in line.split():
                slika.append(int(stevilo.strip()))
    for i in range(visina):
        vrstica = []
        for j in range(sirina):
            vrstica.append(slika[visina*i+j])
        slikca.append(vrstica)
    return slikca
def zapisisliko(stirnk):
    with open("hahahaha.txt","w") as f:
        f.write("test")
        f.write("nova?")
def readpgm(name):
    with open(name) as f:
        lines = f.readlines()

    # Ignores commented lines
    for l in list(lines):
        if l[0] == '#':
            lines.remove(l)

    # Makes sure it is ASCII format (P2)
    assert lines[0].strip() == 'P2' 

    # Converts data to a list of integers
    data = []
    for line in lines[1:]:
        data.extend([int(c) for c in line.split()])

    return (np.array(data[3:]),(data[1],data[0]),data[2])
def transformiraj(x,N):
    def uredi(xx):
        if np.real(xx)>255:
            return 255
        elif np.real(xx)<0:
            return 0
        else:
            return np.real(xx)
    y = np.asarray((np.fft.ifft(np.fft.fft(x)/np.fft.fft(np.vectorize(prenosna)(np.asarray(range(N)),N=N)))).T)
    return np.vectorize(uredi)(y)
def wienerlincoln(x,N,n=0):
    def erm(ermm):
        if ermm<0:
            return 0
        else:
            return ermm
    spektr = np.abs(np.fft.fft(x))**2
    sumcek = np.ones(x.size)*np.sum(spektr[100:200])/100
    signal = spektr - sumcek
    signal = np.vectorize(erm)(signal)
    aux = [1 for i in range(N)] + [0 for i in range(n-2*N)] + [1 for i in range(N)]
    fi =signal/(spektr)*np.asarray(aux)
    return np.fft.ifft((np.fft.fft(x)/np.fft.fft(np.vectorize(prenosna)(np.asarray(range(n))))*fi))
def transformiraj2(x,n,N=60):
    def uredi(xx):
        if np.real(xx)>255:
            return 255
        elif np.real(xx)<0:
            return 0
        else:
            return np.real(xx)
    y = wienerlincoln(x,N,n)
    return np.vectorize(uredi)(y)
def gaussovka(x,mu,sigma):
    return 1/np.sqrt(2*pi*sigma**2) * np.exp(-(x-mu)**2/(2*sigma**2))
def skonvuliraj(x,sigma): 
    N = x.size
    #k = N/2/(sigma)
    iksi = np.array(range(N))
    def pravagaussovka(x):
        if x<=N/2:
            return gaussovka(N/2+x,N/2,sigma)
        elif x>N/2:
            return gaussovka(x-N/2,N/2,sigma)
    def pravagama(x):
        if x<=N/2:
            return gama(N/2+x,k,sigma)
        elif x>N/2:
            return gama(x-N/2,k,sigma)
    #gaussi = np.fft.fft(np.vectorize(pravagama)(iksi))        
    gaussi = np.fft.fft(np.vectorize(pravagaussovka)(iksi))
    #gaussi = np.fft.fft(gaussovka(iksi,mu,sigma))
    return np.real(np.fft.ifft(np.fft.fft(x)*gaussi))
def zgladi(x,sigma):    
    transformiranka = np.apply_along_axis(skonvuliraj,0,x,sigma)
    transformiranka2 = np.apply_along_axis(skonvuliraj,1,transformiranka,sigma)
    return transformiranka2
    
"""    
slika = np.asarray([[[(i+j)%256,(i-j)%256,(i/j)%256]for j in range(1,1921)]for i in range(1080)])
img = Image.fromarray(slika.astype('uint8'))
img.save("champa.jpg","JPEG")
"""

im = Image.open("zid.jpg")
data = np.asarray(im)
print(data.shape)
r = data[:,:,0]
rr = zgladi(r,30)
g = data[:,:,1]
gg = zgladi(g,10)
b = data[:,:,2]
bb = zgladi(b,50)
slika = np.dstack((rr,g,b))
img= Image.fromarray(slika.astype('uint8'))
#img = Image.fromarray(slika)
img.save("bluranzid.jpg","JPEG")

    
    
if 0:
    data = readpgm("tretjaslike/lincoln_L30_N40.pgm")
    slika = np.reshape(data[0],data[1])
    N = data[1][1]
    n = data[1][0]
    def vmesna(N):
        def hahaha(x):
            return transformiraj2(x,n,N)
        return hahaha
    def vmesna2(x):
        return transformiraj(x,n)
    fig , (ax1,ax2,ax3) = plt.subplots(1,3)
    ax1.axis("off")
    ax2.axis("off")
    ax3.axis("off")
    if 0:
        ax1.set_title("N40, Original")
        ax1.imshow(slika,cmap="gray")
        slika1 = np.apply_along_axis(vmesna2,0,slika)
        #plt.suptitle("N40")
        ax2.set_title("N40, Brez filtra")
        ax2.imshow(slika1,cmap="gray")
        slika1 = np.apply_along_axis(vmesna(50),0,slika)    
        ax3.set_title("N40, N=50")
        ax3.imshow(slika1,cmap="gray")    
    if 1:
        slika1 = np.apply_along_axis(vmesna(70),0,slika)    
        ax1.set_title("N40, N=30")
        ax1.imshow(slika1,cmap="gray")    
        slika1 = np.apply_along_axis(vmesna(50),0,slika)    
        ax2.set_title("N40, N=20")
        ax2.imshow(slika1,cmap="gray") 
        slika1 = np.apply_along_axis(vmesna(30),0,slika)    
        ax3.set_title("N40, N=10")
        ax3.imshow(slika1,cmap="gray")
    plt.savefig("tretja/lincoln42.png",bbox_inches="tight")
if 0:
    #DEKONVOLUCIJA S ŠUMOM            
    cm = plt.get_cmap("brg")
    #signal0 = np.loadtxt("signal0.dat")
    signal1 = np.loadtxt("signal1.dat")
    signal2 = np.loadtxt("signal2.dat")
    signal3 = np.loadtxt("signal3.dat")
    iksi = np.asarray(range(512))
    #plt.plot(iksi[:256],np.fft.fft(signal0)[:256],label="signal0",color=cm(0))
    #plt.plot(iksi[:256],np.fft.fft(signal1)[:256],label="signal1",color=cm(0.3))
    #plt.plot(iksi[:256],np.fft.fft(signal2)[:256],label="signal2",color=cm(0.6))
    #plt.plot(iksi[:256],np.fft.fft(signal3)[:256],label="signal3",color=cm(0.9))

    #fourier = np.fft.fft(signal)
    #drugo = np.fft.fft(np.vectorize(prenosna)(iksi))
    #plt.plot(iksi,np.fft.ifft(fourier/drugo))
    #plt.plot(spekter(signal0)[0],spekter(signal0)[1],label="signal0",color = cm(0))
    #plt.plot(iksi,wiener(signal3,100),label="signal3",color=cm(0.6))
    #plt.plot(iksi,wiener(signal2,1),label="signal2",color=cm(0.9))
    plt.plot(iksi,wiener(signal3,50),label="N=50",color = cm(0.3))
    plt.plot(iksi,wiener(signal3,30),label="N=30",color = cm(0.6))
    plt.plot(iksi,wiener(signal3,15),label="N=15",color = cm(0.9))

    #plt.yscale("log")
    plt.ylabel(r"$u(t)$")
    plt.xlabel("t")
    plt.title(r"signal3.dat")
    plt.legend(loc="best")
    plt.savefig("druga/9.pdf")
    

if 0:
    #DEKONVOLUCIJA            
    signal = np.loadtxt("signal0.dat")
    iksi = np.asarray(range(512))
    fourier = np.fft.fft(signal)
    #drugo = np.fft.fft(np.vectorize(prenosna)(iksi))
    plt.plot(spekter(signal)[0],spekter(signal)[1])
    plt.ylabel(r"$|C(f)|^2$")
    plt.yscale("log")
    plt.xlabel("f")
    #plt.title("Rekonstruiran nezašumljen signal")
    plt.savefig("druga/2.pdf")





if 0:
    #SAMO R
    iksi = np.asarray(range(512))
    plt.plot(iksi,np.vectorize(prenosna)(iksi))
    plt.title("r(t)")
    plt.xlabel("t")
    plt.ylabel("r")
    plt.savefig("druga/r.pdf")

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
























