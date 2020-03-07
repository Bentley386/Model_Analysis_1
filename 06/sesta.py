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
from scipy.stats import linregress
from scipy.linalg import svd
plt.rc("text",usetex=True)
matplotlib.rcParams["text.latex.unicode"] = True
plt.close("all")
#izmerki = np.loadtxt("farmakoloski.dat")
izmerki = np.loadtxt("thtg-xfp-thfp.dat") #target, x fp, thetafp
#izmerki = np.loadtxt("CdL3_linfit.norm")
#print(izmerki)
#x doza y odziv
def hikvadrat(izmerki, model, parametri, napake):
    rezid = (izmerki - model(izmerki,parametri))**2
    return rezid/(napake**2)
def narediA(izmerki,napake):
    matrika = []
    for i,j in zip(izmerki[:,0], napake):
        matrika.append([1/j,i/j])
    return np.matrix(matrika)
def narediAiksi(izmerki,napake,funkcije):
    matrika = []
    for i,j in zip(izmerki,napake):
        vrstica = [f(i)/j for f in funkcije]
        matrika.append(vrstica)
    return np.matrix(matrika)
def narediAmesani(izmerki1,izmerki2,napake,funkcije):
    matrika = []
    napaka = napake[0]
    for i,j in zip(izmerki1,izmerki2):
        vrstica = [f(i,j)/napaka for f in funkcije]
        matrika.append(vrstica)
    return np.matrix(matrika)
def ita(dekompozicija,b,i):
    return np.sum(dekompozicija[0][:,i]*b)/dekompozicija[1][i] * dekompozicija[2][i]
def naredipotencno(p):
    def potencna(x):
        return x**p
    return potencna
def narediceb(p):
    def ceb(x):
        return spec.eval_chebyt(p,x/200)
    return ceb
def narediceb2(p):
    def ceb2(x):
        return spec.eval_chebyt(p,x/95)
    return ceb2    
def naredileg(p):
    def leg(x):
        return spec.eval_legendre(p,x/200)
    return leg
def naredilag(p):
    def lag(x):
        return spec.eval_laguerre(p,x+200)
    return lag
def narediher(p):
    def her(x):
        return spec.eval_hermite(p,x)
    return her
def naredimesano(p1,p2):
    def mesana(x,y):
        return (x)**p1 * (y)**p2
    return mesana
def naredimesanopoli(p1,p2,p):
    def mesana(x,y):
        return (x)**p1 * (y)**p2
    def mesanaceb(x,y):
        return spec.eval_chebyt(p1,x/200)*spec.eval_chebyt(p2,y/200)
    def mesanaleg(x,y):
        return spec.eval_legendre(p1,x/200)*spec.eval_legendre(p2,y/200)
    def mesanalag(x,y):
        return spec.eval_laguerre(p1,x+200)*spec.eval_laguerre(p2,y+200)
    def mesanaher(x,y):
        return spec.eval_hermite(p1,x)*spec.eval_hermite(p2,y)
    if p==0:
        return mesana
    elif p==1:
        return mesanaceb
    elif p==2:
        return mesanaleg
    elif p==3:
        return mesanalag
    elif p==4:
        return mesanaher
        
def naredimesano2(p1,p2):
    def mesana2(x,y):
        return (x**p1 * y)**p2
    return mesana2
def naredimesano3(p1,p2):
    def mesana3(x,y):
        return(y**p1 * x)**p2
    return mesana3
def evaluiraj(x,y,funkcije1,funkcije2,parametri):
    vrednost = 0
    for i,j in zip(funkcije1,parametri[0]):
        vrednost += j*i(x)
    for i,j in zip(funkcije2,parametri[1]):
        vrednost +=j*i(y)
    return vrednost
def evaluiraj2(x,y,funkcije1,funkcije2,funkcije3,parametri): #še mešani členi
    vrednost = 0
    for i,j in zip(funkcije1,parametri[0]):
        vrednost += j*i(x)
    for i,j in zip(funkcije2,parametri[1]):
        vrednost +=j*i(y)
    for i,j in zip(funkcije3,parametri[2]):
            vrednost +=j*i(x,y)
    return vrednost
#napake1=np.ones(izmerki[:,0].size)*0.06 #iksi
napake1 = np.ones(izmerki[:,0].size)*0.06
N = izmerki[:,0].size    
def naredifunkcije(spodnja,zgornja,polinom):
    if polinom == 0:
        return [naredipotencno(i) for i in range(spodnja,zgornja)]
    elif polinom == 1:
        return [narediceb(i) for i in range(spodnja,zgornja)]
    elif polinom == 12:
        return [narediceb2(i) for i in range(spodnja,zgornja)]        
    elif polinom == 2:
        return [naredileg(i) for i in range(spodnja,zgornja)]
    elif polinom == 3:
        return [naredilag(i) for i in range(spodnja,zgornja)]
    elif polinom == 4:
        return [narediher(i) for i in range(spodnja,zgornja)]
def reduciranhipotence(px,py,izmerki):
    funkcijex = naredifunkcije(0,px+1,4)
    #funkcijex = naredifunkcije(0,px+1,polinom1) #normalen, ceb, leg , lag ,her
    funkcijey = naredifunkcije(1,py+1,4)
    #funkcijey = naredifunkcije(1,py+1,polinom2)
    #funkcijemesane =[naredimesano2(i,j) for i in range(1,pmesan1+1) for j in range(1,pmesan2+1)]
    #funkcijemesane = funkcijemesane + [naredimesano(1,1)] #dodam normalen xy člen
    funkcijemesane = [naredimesanopoli(i,i,4) for i in range(1,4)]
    funkcijemesane = funkcijemesane + [naredimesanopoli(2*i,i,4) for i in range(1,3)]
    funkcijemesane = funkcijemesane + [naredimesanopoli(3*i,i,4) for i in range(1,2)]
    funkcijemesane = funkcijemesane + [naredimesanopoli(i,2*i,4) for i in range(1,3)]
    funkcijemesane = funkcijemesane + [naredimesanopoli(i,3*i,4) for i in range(1,2)]
    funkcijemesane = funkcijemesane + [naredimesanopoli(-i,-i,4) for i in range(1,4)]                                       
    Ax = narediAiksi(izmerki[:,1],napake1,funkcijex) # xfp
    Ay = narediAiksi(izmerki[:,2],napake1,funkcijey) # theta fp
    Amesan = narediAmesani(izmerki[:,1],izmerki[:,2],napake1,funkcijemesane) # theta fp
    A = np.hstack((Ax,Ay,Amesan))
    b = izmerki[:,0]/napake1
    dekompozicija = svd(A,full_matrices=False)    
    m = px+py+1+len(funkcijemesane) #+1 na koncu zaradi xy
    suma = np.zeros(m)
    print(len(funkcijemesane))
    print("ja")
    for i in range(m):
        suma+= ita(dekompozicija,b,i)
    parametri=[suma[:(px+1)],suma[(px+1):(px+1+py)],suma[(px+1+py):]]
    ipsiloni = [evaluiraj2(izmerki[:,1][i],izmerki[:,2][i],funkcijex,funkcijey,funkcijemesane,parametri) for i in range(izmerki[:,0].size)]
    ipsiloni = np.asarray(ipsiloni)
    hi2 = (izmerki[:,0]-ipsiloni)**2/(napake1**2)
    return [np.sum(hi2)/(N-m),parametri]
def varianca(dekompozicija,j):
    V = dekompozicija[2].T
    s = dekompozicija[1]
    suma = 0
    for i in range(len(s)):
        suma+= V[i][j]**2 /(s[i]**2)
    return suma
def kovarianca(dekompozicija,j,k):
    V = dekompozicija[2].T
    s = dekompozicija[1]
    suma = 0
    for i in range(len(s)):
        suma+= V[i][j]*V[i][k]/(s[i]**2)
    return suma

a= reduciranhipotence(4,4,izmerki)  
print(a[0])
"""
hikvadrati = [reduciranhipotence(4,4,3,i,izmerki) for i in range(1,2)]
plt.plot(range(1,2),hikvadrati,"o",)
#plt.legend(loc="best")
plt.xlabel("s")
#plt.yscale("log")
plt.ylabel(r"$\chi^2$")
plt.title(r"$\chi^2$ še z mešanimi členi. $t=2$")
#plt.savefig("druga/poblize2.pdf")
plt.show()    
"""    
"""
#napake2=np.ones(izmerki[:,0].size)*0.18/np.pi #thete
print(np.amin(izmerki[:,2]))
maks = 7
hikvadrati = [reduciranhipotence(i,i,0,izmerki,0,0) for i in range(2,maks)]
plt.plot(range(2,maks),hikvadrati,"o",label=r"Navadne potence")
hikvadrati = [reduciranhipotence(i,i,0,izmerki,1,1) for i in range(2,maks)]
plt.plot(range(2,maks),hikvadrati,"o",label=r"Čebišev")
hikvadrati = [reduciranhipotence(i,i,0,izmerki,2,2) for i in range(2,maks)]
plt.plot(range(2,maks),hikvadrati,"o",label=r"Legendre")
hikvadrati = [reduciranhipotence(i,i,0,izmerki,3,3) for i in range(2,maks)]
plt.plot(range(2,maks),hikvadrati,"o",label=r"Laguerre")
hikvadrati = [reduciranhipotence(i,i,0,izmerki,4,4) for i in range(2,maks)]
plt.plot(range(2,maks),hikvadrati,"o",label=r"Hermite")
plt.xlabel("s")
plt.title(r"$\chi^2$ z več vrstami polinomov. $q=1$")
plt.ylabel(r"$\chi^2$")
#plt.yscale("log")
plt.legend(loc="best")
plt.savefig("druga/specialni2.pdf")
plt.show()
"""
























"""
#krovna plast, sredica, zveplo, kisik    
funkcije = naredifunkcije(1,2,0)
funkcijemesane = [[]]
Ax = narediAiksi(izmerki[:,3],napake1,funkcije)
Ay = narediAiksi(izmerki[:,4],napake1,funkcije)
A = np.hstack((Ax,Ay))
b = izmerki[:,2]/napake1
dekompozicija = svd(A,full_matrices=False)
m = 2
suma = np.zeros(m)
for i in range(m):
    suma+= ita(dekompozicija,b,i)
parametri = [[suma[0]],[suma[1]],[]]
ipsiloni = [evaluiraj2(izmerki[:,3][i],izmerki[:,4][i],funkcije,funkcije,funkcijemesane,parametri) for i in range(izmerki[:,0].size)]
ipsiloni = np.asarray(ipsiloni)
hi2 = (izmerki[:,1]-ipsiloni)**2/(napake1**2)
hi2 = np.sum(hi2)/(N-m)
vara = varianca(dekompozicija,0)
varb = varianca(dekompozicija,1)
covar = kovarianca(dekompozicija,0,1)
iksi = izmerki[:,0]
zveplo = izmerki[:,3]
kisik = izmerki[:,4]
krovna = izmerki[:,1]    
sredica = izmerki[:,2]
plt.title("Sredica lističa")
plt.xlabel("Energija")
plt.ylabel("Absorbcija")
tekst = r"$\\ \chi^2 = {0}$\\$a = {1}\\b={2}\\ \sigma^2(a)={3}\\ \sigma^2(b)={4}\\ cov(a,b)={5}\\ \rho = {6}$".format(round(hi2,5),round(parametri[0][0],5),round(parametri[1][0],5),round(vara,5),round(varb,5),round(covar,5),round(covar/(vara**0.5 * varb**0.5),5))
#plt.plot(iksi,zveplo,"--",label="Vezan na žveplo")
#plt.plot(iksi,kisik,"--",label="Vezan na kisik")  
#plt.plot(iksi,sredica,label="Krovna plast lističa")
plt.plot(iksi,sredica,".")
plt.plot(iksi,ipsiloni)
plt.text(3560,0.3,tekst)
#plt.legend(loc="best")
#plt.savefig("tretja/sredica2.pdf")
plt.show()  
"""

"""
fig = plt.figure()
ax = fig.gca(projection="3d")
ax.scatter(izmerki[:,1],izmerki[:,2],izmerki[:,0],s=1)
#plt.plot(izmerki[:,2],izmerki[:,0],"o")
plt.xlabel(r"$x_{fp}$")
plt.ylabel(r"$\theta_{fp}$")
#ax.set_zlabel(r"$\theta_{tg}$")
plt.savefig("druga/scatter.pdf")
plt.show()
"""
def prva(izmerki,napake): #resi y=y0 + ax
    iksi = izmerki[:,0]
    ipsiloni = izmerki[:,1]
    A = napake**(-2)
    Ax = iksi*A
    Axx = Ax*iksi
    Axy = Ax*ipsiloni
    Ay = ipsiloni*A
    Ayy = Ay*ipsiloni
    y0=(np.sum(Axx)*np.sum(Ay) - np.sum(Ax)*np.sum(Axy))/(np.sum(Axx)*np.sum(A)-np.sum(Ax)**2)
    a=(np.sum(A)*np.sum(Axy) - np.sum(Ax)*np.sum(Ay))/(np.sum(Axx)*np.sum(A)-np.sum(Ax)**2)
    vary0 = np.sum(Axx)/(np.sum(Axx)*np.sum(A)-np.sum(Ax)**2)
    vara = np.sum(A)/(np.sum(Axx)*np.sum(A)-np.sum(Ax)**2)
    cov = -np.sum(Ax)/(np.sum(Axx)*np.sum(A)-np.sum(Ax)**2)
    return [y0,a,vary0,vara,cov]
"""
#napake = np.ones(8)*3*izmerki[:,1]**(-2)
napake = np.ones(8)
A = narediA(izmerki**(-1),napake)
dekompozicija = svd(A,full_matrices=False)
print(dekompozicija[2])
b = izmerki[:,1]**(-1)/(napake)
suma = np.zeros(2)
for i in range(2):
    suma+= ita(dekompozicija,b,i)

napake = np.ones(8)*3*izmerki[:,1]**(-2)
#rezultati = prva(izmerki**(-1),napake)
rezultati = np.polyfit(izmerki[:,0]**(-1),izmerki[:,1]**(-1),1,w=napake**(-1))
#suma = [rezultati[1],rezultati[0]]
#suma=[1/suma[0],suma[1]/suma[0]] #y0, a suma je prej bla 1/y0,drugo
hi = np.sum((izmerki[:,1]**(-1)-(rezultati[1]+rezultati[0]*izmerki[:,0]**(-1)))**2/(napake**2))
print(hi)
print(suma)
x = np.linspace(0,1000,5000)
plt.plot(izmerki[:,0],izmerki[:,1],"ro")
#print(suma[0]*x/(x+suma[1]))
plt.plot(x,suma[0]*x/(x+suma[1]))
#plt.text(0.5,0.5,r"$\chi^2 = 26.08$")
plt.xlabel("1/x")
plt.ylabel("1/y")
#plt.savefig("prva/netransf.pdf")
plt.show()
"""





