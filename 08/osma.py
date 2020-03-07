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

def timer(code):
    start_time = timeit.default_timer()
    exec(code)
    elapsed = timeit.default_timer() - start_time
    return elapsed
    
def hikvadrat(izmerki, model, parametri, napake):
    rezid = (izmerki - model(izmerki,parametri))**2
    return rezid/(napake**2)
def stetje(mejefi,mejetheta,x):
    counter = 0
    for i in x:
        if i[0] < mejefi[1] and i[0] > mejefi[0] and i[1] < mejetheta[1] and i[1] > mejetheta[0]:
            counter+=1
    return counter
def boxmuller(size=1):
    u1, u2 = np.random.rand(2)
    z0 = np.sqrt(-2*np.log(u1))*np.cos(2*np.pi*u2)
    z1 = np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)
    return (z0,z1)
def muller(size=1):
    return [boxmuller()[0] for i in range(size)]
def mullern(N):
    rezultati = []
    for i in range(int(N/2)):
        rezultati.append(boxmuller())
    return rezultati
def analgauss(x):
    return 1/(np.sqrt(2*np.pi)) * np.exp(-x*x /2)
def gaussverjetnost(a,b):
    return 0.5*(spec.erf(b/np.sqrt(2))-spec.erf(a/np.sqrt(2)))
def gausskum(x,size=1):
    return 0.5*(1+spec.erf(x/np.sqrt(2)))
def konvolucijska(size=1):
    if size==1:
        return np.sum(np.random.rand(6)) - np.sum(np.random.rand(6))
    else:
        return [konvolucijska() for i in range(size)]
def naredibine(nn):
    n = nn+1
    razredi = [(-100,-3)]
    x = np.linspace(-3,3,n-2)
    for i in range(0,x.size-1):
        razredi.append((x[i],x[i+1]))
    razredi.append((3,100))
    return razredi
def spremeni(x):
    spremenjeno = []
    for i in x:
        stevilke = i.split(":")
        #print(stevilke)
        cas = float(stevilke[2])+float(stevilke[1])*60 + float(stevilke[0])*60*24
        spremenjeno.append(cas)
    return spremenjeno
def generirajfi(N):
    return 2*np.pi*np.random.rand(N)
def generirajtheta(N):
    return np.arccos(2*np.random.rand(N)-1)
def generirajtheta2(N):
    def funk(x,p):
        return np.sin(x/2)**4 * (np.cos(x)+2) - p
    x0 = np.ones(N)*np.pi/2
    pji = np.random.rand(N)
    resitev = []
    for i in range(N):
        resitev.append(fsolve(funk,x0[i],args=pji[i])[0])
    return np.asarray(resitev)
"""
dji = []    
for i in ["01","02","03","04","05","06"]:
    prva = []
    with open("mod_times/podatki/mod_tm11_1{0}.dat".format(str(i)),"r") as f:
        for line in f:
            prva.append(line)
    prva = spremeni(prva)
    druga = []
    with open("mod_times/podatki/mod_tm13_1{0}.dat".format(str(i)),"r") as f:
        for line in f:
            druga.append(line)
    druga = spremeni(druga)
    dji.append(scipy.stats.ks_2samp(prva,druga)[0])
plt.plot([1,2,3,4,5,6],dji,"o")
plt.xlabel("Zaporedna številka naloge")
plt.ylabel("D")
plt.title("Primerjava časov oddaje nalog v letih 11 in 13")
plt.savefig("tretja/1113.pdf")
plt.show()
"""
#stetje praznih celic - np.where
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_aspect(1)
# Make data.
#print(np.sum(scipy.special.sph_harm(0,1,fiji,np.zeros(10000)))/10000)
#thete = generirajtheta2(10000)
#plt.plot(x[:,0],x[:,1],"o")
#plt.title(r"Dipolno sevanje, N=100")
#plt.xlabel(r"Številka intervala $\phi$")
#plt.ylabel(r"Številka intervala $\cos\theta$")
#x = np.sin(thete)*np.cos(fiji)
#y = np.sin(thete)*np.sin(fiji)
#z = np.cos(thete)
#ax.scatter(x, y, z,".",s=1)
#plt.title("Dipolno sevanje, N=10000")
#plt.savefig("druga/dipol3.pdf")
#plt.hist(fiji)
#plt.title(r"Dipolno prostorsko porazdeljeni $\theta, N=100000$")
#plt.savefig("druga/graf3.pdf")
#plt.show()    

#VARIACIJE V MOMENTIH
varcos = []
varcos2 = []
varcos3= []
vary1 = []
vary2 = []
vary3 = []
for n in np.linspace(100,10000,15):
    N = int(n)
    povp = []
    for i in range(10):
        thete = generirajtheta(N)
        povp.append(np.sum(np.cos(thete))/N)
    povp = np.asarray(povp)
    varcos.append(np.sum((povp - np.sum(povp)/10)**2)/10)
    povp = []
    for i in range(10):
        thete = generirajtheta(N)
        povp.append(np.sum(np.cos(thete)**2)/N)
    povp = np.asarray(povp)
    varcos2.append(np.sum((povp - np.sum(povp)/10)**2)/10)    
    povp = []
    for i in range(10):
        thete = generirajtheta(N)
        povp.append(np.sum(np.cos(thete)**3)/N)
    povp = np.asarray(povp)
    varcos3.append(np.sum((povp - np.sum(povp)/10)**2)/10)
    povp = []
    for i in range(10):
        thete = generirajtheta(N)
        povp.append(np.sum(spec.sph_harm(0,1,thete,np.zeros(N)))/N)
    povp = np.asarray(povp)
    vary1.append(np.sum((povp - np.sum(povp)/10)**2)/10)
    povp = []
    for i in range(10):
        thete = generirajtheta(N)
        povp.append(np.sum(spec.sph_harm(0,2,thete,np.zeros(N)))/N)
    povp = np.asarray(povp)
    vary2.append(np.sum((povp - np.sum(povp)/10)**2)/10)
    povp = []
    for i in range(10):
        thete = generirajtheta(N)
        povp.append(np.sum(spec.sph_harm(0,3,thete,np.zeros(N)))/N)
    povp = np.asarray(povp)
    vary3.append(np.sum((povp - np.sum(povp)/10)**2)/10)
N = np.linspace(100,10000,15)
#plt.plot(N,varcos,"--.",label=r"$\cos\theta$")
#plt.plot(N,varcos2,"--.",label=r"$\cos^2\theta$")
#plt.plot(N,varcos3,"--.",label=r"$\cos^3\theta$")
plt.plot(N,vary1,"--.",label=r"$Y_{10}(\theta)$")
plt.plot(N,vary2,"--.",label=r"$Y_{20}(\theta)$")
plt.plot(N,vary3,"--.",label=r"$Y_{30}(\theta)$")

plt.legend(loc="best")
plt.title("Variacije v momentih za enakomerno porazdelitev")    
plt.xlabel("N")
#plt.savefig("Druga/variacije1.pdf")
plt.show()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
"""
fiji = np.linspace(0,2*np.pi,10)
intfi = []
for i in range(9):
    intfi.append((fiji[i],fiji[i+1]))
thete = np.linspace(-1,1,10)
intth = []
for i in range(9):
    intth.append((thete[i],thete[i+1]))
"""    
    
    
"""
x = []
for i in range(20):
    N = 10000
    fiji = generirajfi(N)
    thete = np.cos(generirajtheta(N))
    pojavitve = []
    fith = []
    pari = np.asarray([[fiji[i],thete[i]] for i in range(N)])
    for i in range(9):
        for j in range(9):
            koliko = stetje(intfi[i],intth[j],pari)
            pojavitve.append(koliko)
            fith.append((i,j))
    x.append(fith[pojavitve.index(min(pojavitve))])
x = np.asarray(x)
"""    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#for i in [10,50,100,500]:
    #y = []
    #for n in np.linspace(500,500000,15):
        #x = np.asarray([boxmuller() for i in range(int(n/2))])
        #x = np.ndarray.flatten(x)
        #bini = np.asarray(naredibine(i))
        #verjetnosti = np.asarray([gaussverjetnost(i[0],i[1]) for i in bini])
        #x = np.histogram(x,np.append(bini[:,0],100))[0]
        #izmerki = [scipy.stats.chisquare(x,verjetnosti*n)[0] for i in range(10)]
        #y.append(sum(izmerki)/(10*(i-1)))
    #plt.plot(np.linspace(500,500000,15),y,"--.",label="n={0}".format(str(i)))
#izmerki = [scipy.stats.kstest(muller,gausskum,N=500)[1] for i in range(100)]
#plt.hist(izmerki)
#plt.title("Kolmogorov-Smirnov test za konvolucijski generator")
#plt.ylabel(r"D")
#plt.xlabel("N")
"""
N = 10000
casimuller = [timer("mullern({0})".format(int(i))) for i in N]
plt.plot(N,casimuller,"--o",label="Box-Muller")
casimuller = [timer("konvolucijska({0})".format(int(i))) for i in N]
plt.plot(N,casimuller,"--o",label="Konvolucijski")
plt.legend(loc="best")
plt.title("Primerjava hitrosti izvajanje generatorjev")
plt.xlabel("N")
plt.ylabel(r"Čas [ms]")

plt.savefig("prva/hitrosti.pdf")
plt.show()
"""















"""
bini = naredibine(100)
bini = np.asarray(bini)
#razredi = [(-100,-3),(-3,-2.5),(-2.5,-2),(-2,-1.5),(-1.5,-1),(-1,-0.5),(-0.5,0),(0,0.5),(0.5,1),(1,1.5),(1.5,2),(2,2.5),(2.5,3),(3,100)]
verjetnosti = np.asarray([gaussverjetnost(i[0],i[1]) for i in bini])
#bini = [-100,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,100]
iksi = []
for i in range(1000):
    x = np.asarray([boxmuller() for i in range(5000)])
    #x = np.random.normal(0,1,1000000)
    x = np.ndarray.flatten(x)
    #scipy.stats.probplot(x,plot=plt)
    #print(scipy.stats.normaltest(x))
    x = np.histogram(x,np.append(bini[:,0],100))[0]
    hi2 = scipy.stats.chisquare(x,verjetnosti*10000)[1]
    iksi.append(hi2)
#x = [scipy.stats.kstest(konvolucijska,gausskum,N=10000)[1] for i in range(100)]
plt.hist(iksi)
razredi = np.asarray([(0,0.2),(0.2,0.4),(0.4,0.6),(0.6,0.8),(0.8,1)])
print(scipy.stats.chisquare(np.histogram(iksi,np.append(razredi[:,0],1))[0],np.ones(5)*0.2)[0])
plt.ylabel("P vrednost")
plt.title(r"$\chi^2, N = 10000$, Box-Muller")
plt.savefig("prva/histogram.pdf")
plt.show()
#print(scipy.stats.kstest(konvolucijska,gausskum,N=n)[1])    
#print(hi2[1])   
"""    
    
    
"""
RISANJE HISTOGRAMOV
x = [boxmuller() for i in range(500)]
x = np.ndarray.flatten(np.asarray(x))    
plt.hist(x,normed=True,bins=50,label=r"$N=500$")
x = [boxmuller() for i in range(1000)]
x = np.ndarray.flatten(np.asarray(x))    
plt.hist(x,normed=True,bins=50,label=r"$N=1000$")
x = [boxmuller() for i in range(10000)]
x = np.ndarray.flatten(np.asarray(x))    
plt.hist(x,normed=True,bins=50,label=r"$N=10000$")
x = [boxmuller() for i in range(100000)]
x = np.ndarray.flatten(np.asarray(x))    
plt.hist(x,normed=True,bins=50,label=r"$N=100000$")
x = [boxmuller() for i in range(1000000)]
x = np.ndarray.flatten(np.asarray(x))    
plt.hist(x,normed=True,bins=50,label=r"$N=1000000$")
x = np.linspace(-5,5,1000)
plt.plot(x,np.vectorize(analgauss)(x),label="Normalna porazdelitev")
plt.title("Box-Mullerjev generator")
plt.xlim(-4,4)
plt.legend(loc="best")
plt.savefig("prva/muller.pdf")
plt.show()
"""    
    
    
    
    
    
    
    
    
    
    
"""
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


