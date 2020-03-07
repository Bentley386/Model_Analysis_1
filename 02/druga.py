# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 20:58:56 2016

@author: Admin
"""
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


data =np.genfromtxt("tabela-zivil.dat",usecols=range(1,10),dtype=None) #na 100g 
l = len(data)
#print(data)
data = [[data[i][j] for j in range(9)] for i in range(l)]
data = np.around(data,decimals=3)
c = [data[i][0] for i in range(len(data)) if i not in [15,16,17,18,29,32,33,34,35,36]]
#print(len(c))
cc = [data[i][1] for i in range(len(data)) if i not in [15,16,17,18,29,32,33,34,35,36]]
data1 = [[-data[i][j] for j in range(1,9)] for i in range(l) if i not in [15,16,17,18,29,32,33,34,35,36]] #min kalorij
data1 = np.around(data1,decimals=3)
aub = np.matrix(data1).T
data2 = [[-data[i][j] for j in (0,2,3,4,5,6,7,8)] for i in range(l) if i not in [15,16,17,18,29,32,33,34,35,36]] #min mascob
data2 = np.around(data2,decimals=3)
data3 = [[-data[i][j] for j in range(9)] for i in range(l) if i not in [15,16,38,47,48]] #min cene
data3 = np.around(data3,decimals=3)
aubbb = np.matrix(data3).T
aubb = np.matrix(data2).T
#print(aub)
extrarow = [1]*l
aub = np.vstack([aub,[1]*(l-10)])
aub = np.vstack([aub,-aub[-2]])
aubb = np.vstack([aubb,[1]*(l-10)])
aubb = np.vstack([aubb,-aubb[-2]])
aubbb = np.vstack([aubbb,[1]*(l-5)])
aubbb = np.vstack([aubbb,-aubbb[-2]])
b = np.matrix([-70,-310,-50,-1000,-18,-60,-3500,-500,20,2400]).T
bb = np.matrix([-2000,-310,-50,-1000,-18,-60,-3500,-500,20,2400]).T 
bbb = np.matrix([-2000,-70,-310,-50,-1000,-18,-60,-3500,-500,20,2400]).T 
#bb = np.matrix([-70,-310,-50,-1000,-18,-60,-3500,-500,20,2400]).T               


print(l-4)

"""
mascobe gljikovi hidrati[g]	proteini[g]	Ca[mg]	Fe[mg]	Vitamin C[mg]	Kalij[mg] Natrij[mg]
ce je priporoˇcen minimalni dnevni
vnos 70 g maˇsˇcob, 310 g ogljikovih hidratov, 50 g proteinov, 1000 mg kalcija ter 18 mg ˇzeleza.
"""
"""
cene=[1.09/5,0.14,0.19,3.59/4.2,2.29/4,0.35,0.08,1.89/2.5,0.25,1.99/4,1.59/2.5,
1.19/1.1,2.49,0.3,0.89/80*100,0.59/4,0.1,1.19,0.05,3.96/3,0.23,2.29/2.5,
1.48/7.5,2.19/2.5,3.59/3.5,4.99/7.5,2.29/50,0.1,0.55,2.51/10,4.19/3.2,
2.59/10,2.22/10,0.08,1.09/4.5,0.59/4,0.1,0.35,10.5/10,1.49/1.5,6.2/10,
0.7,2.14/3.5,6.99/10,1.79/4,0.11,0.35/5,5.46/7.5,4.59/7.5]    
#cene = [cene[i] for i in range(len(cene)) if i not in [15,16,17,18,29,32,33,34,35,36]]
#se 60 mg kalcij    
"""
cene =[7.9,198/16,379/32,628/28,1162/16,26.3,26.8,2499/32,
38.9,391/16,499/8,238/15,20.3,3.1,39.3,259/24,12.4,235/16,211/16,269/16,
289/15,180/16,1399/12.2,33.7,16.5,87/16,82.2,319/16,19.1,11.3,1329/16, 10.5,
10.3,5.5,11.4,549/(3*16),783/16,121,29.9,10.9,15.6,55.5,6.4,
5.2]
cene = [i*3.5*0.85/100 for i in cene] #v eur/100g    
print(cene)
#print(len(cene))
#y = linprog(c,aub,b,bounds=(0,2)).x
yy = linprog(cene,aubbb,bbb,bounds=(0,2)).x
#yx= [i for i in range(len(y)) if y[i]!=0]
yyx= [i for i in range(len(yy)) if yy[i]!=0]
#print(yx)
print(yyx)
#y = [y[i] for i in range(len(y)) if y[i]!=0]
yy = [yy[i] for i in range(len(yy)) if yy[i]!=0]
#print(yy[-1])
#x = np.arange(len(y))
xx = np.arange(len(yy))
#bar1=plt.bar(x,y,color="orange")
bar2=plt.bar(xx,yy,color="green")
plt.ylabel("Količina [100g]")
plt.yticks([i/4 for i in range(11)])
tiki = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]
tiki2=["ovseni kosmiči","edamec","skuša","riž","čokolada","banana","korenje","fižol","pomaranča"]
#tiki3=["bel kruh","fižol","solata","puran","marmelada"]
plt.xticks(tiki,tiki2,rotation=45,ha="right")
plt.tick_params(axis="x",length=0)
cena = 0
counter = 0
for rect in bar2: #NE POZABIT BAR1
    cenca = cene[yyx[counter]]
    cena+=cenca
    counter += 1
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width()/2.0, height+0.2, str(round(cenca,2)), ha='center', va='bottom')
    plt.text(rect.get_x() + rect.get_width()/2.0, height, str(round(height,2)), ha='center', va='bottom')
    
#plt.text(6.5,2.2,"Skupna cena: \n"+str(round(cena,2))+"€")
    #plt.text(bar2[-1].get_x() + bar2[-1].get_width()/2.0, bar2[-1].get_height(), str(round(bar2[-1].get_height(),3)), ha='center', va='bottom')    
plt.grid(True,axis="y")
#plt.legend(loc=(0.01,0.75))
plt.title("Minimizacija cene. Skupna cena:" + str(round(cena,2))+"€") 
plt.tight_layout()
plt.savefig("usa.pdf")    

    
    