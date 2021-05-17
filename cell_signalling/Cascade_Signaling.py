
# **Study of MAPK dynamics following Fig. 1 from B.N. Khodolenko, J.F. Hancock, W. Kolch, Signalling ballet in space and time, Nat Rev Mol Cell Biol 11, pp.414-426 (2010).**

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint

avo = 6.022*10**(23)       ## 1/mol
V = 1000*10**(-15)         ## L

## parameters changed into #-of-molecules/s or #-of-molecules 
F = 5
kf = 100         ## nM
k1 = 0.025        ## 1/s
k11 = 0.025     ## nM/s
k7 = 3.0        ## nM/s
rasgtp = 10     ## nM

km1 = 100
km2 = 200
k2 = 0.25
k3 = 2.5
km3 = 50
km4 = 100
k4 = 3.75
k5 = 2.5
km5 = 250
km6 = 250
k6 = 0.5
km7 = 250
km8 = 80
k8 = 3.75
k9 = 0.125
km9 = 250
km10 = 250
k10 = 0.125
km11 = 120
km12 = 20
k12 = 5
km13 = 300


#initial conditions
y0 =  [0] * 9
y0[0] = y0[3] = y0[6] = 300.0
y0[8] = 45
y0[6] = 300.0 - y0[8]

print(y0)

#Derivatives
def dy_dt(y,t):
    v1= ((k1*rasgtp*y[0]/km1)/(1+y[0]/km1+y[1]*km2))*((1+F*y[8]/kf)/(1+y[8]/kf))
    v2= ((k2*rasgtp*y[1]/km2)/(1+y[0]/km1+y[1]*km2))*((1+F*y[8]/kf)/(1+y[8]/kf))
    v3= (k3*y[2]/km3)/(1 + (y[3]/km3) + (y[1]/km4))
    v4= (k4*y[1]/km4)/(1 + (y[2]/km3) + y[1]/km4)

    v5= (k5*y[2]*y[3]/km5)/(1 + (y[3]/km5) + (y[4]/km6))
    v6= (k6*y[2]*y[4]/km6)/(1 + (y[3]/km5) + (y[4]/km6))
    v7= (k7*y[5]/km7)/(1 + (y[5]/km7) + (y[4]/km8))
    v8= (k8*y[4]/km8)/(1 + (y[5]/km7) + (y[4]/km8))

    v9= (k10*y[5]*y[6]/km9)/(1 + (y[6]/km9) + (y[7]/km10))
    v10= (k9*y[5]*y[7]/km10)/(1 + (y[6]/km9) + (y[7]/km10))
    v11= (k11*y[8]/km11)/(1 + (y[8]/km11) + (y[7]/km12))
    v12= (k12*y[7]/km12)/(1 + (y[8]/km11) + (y[7]/km12))
    
    dy = []
    dy.append(v4-v1)
    dy.append(v1+v3-v2-v4)
    dy.append(v2-v3)
    dy.append(v8-v5)
    dy.append(v5+v7-v6-v8)
    dy.append(v6-v7)
    dy.append(v12-v9)
    dy.append(v9+v11-v10-v9)
    dy.append(v10-v11)
    return dy    

# raf,praf,ppraf,mek,pmek,ppmek,erk,perk,pperk


#Integration of ODEs
t = np.linspace(0, 120*60,10000)
sol = odeint(func=dy_dt, y0=y0, t=t)
print(sol.shape)
print(type(sol))
sol = pd.DataFrame(sol, columns=["Raf", "pRaf","ppRaf","MEK","pMEK","ppMEK","ERK","pERK","ppERK"])

sol.index = t
sol.index = sol.index/60

#Plotting
fig, ax = plt.subplots(figsize=(8,6))
sol.plot(subplots=True,layout=(3,3),figsize=(10,6),ax=ax, xlabel="Time (min)", ylabel = "Concentration (nM)")
fig.subplots_adjust(top=0.8)
fig.align_ylabels()
plt.tight_layout()
plt.savefig("ODE")
