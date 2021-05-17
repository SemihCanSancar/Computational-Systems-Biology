import numpy as np
import math
import random
import matplotlib.pyplot as plt

def gillespie(tsim):

    ## constant

    avo = 6.022*10**(23)       ## 1/mol
    V = 1000*10**(-15)         ## L

    ## parameters changed into #-of-molecules/s or #-of-molecules 
    F = 0.01
    kf = 25         ## nM
    kf = (kf*10**(-9))*(avo*V)
    k1 = 1.0        ## 1/s
    k11 = 3.75      ## nM/s
    k11 = (k11*10**(-9))*(avo*V)
    k7 = 3.0        ## nM/s
    k7 = (k7*10**(-9))*(avo*V)
    rasgtp = 10     ## nM
    rasgtp = (rasgtp*10**(-9))*(avo*V)

    km1 = 100
    km1 = (km1*10**(-9))*(avo*V)
    km2 = 200
    km2 = (km2*10**(-9))*(avo*V)
    k2 = 0.25
    k3 = 2.5
    k3 = (k3*10**(-9))*(avo*V) 
    km3 = 50
    km3 = (km3*10**(-9))*(avo*V)
    km4 = 100
    km4 = (km4*10**(-9))*(avo*V)
    k4 = 3.75
    k4 = (k4*10**(-9))*(avo*V)
    k5 = 2.5
    km5 = 250
    km5 = (km5*10**(-9))*(avo*V)
    km6 = 250
    km6 = (km6*10**(-9))*(avo*V)
    k6 = 0.5
    km7 = 250
    km7 = (km7*10**(-9))*(avo*V)
    km8 = 80
    km8 = (km8*10**(-9))*(avo*V)
    k8 = 3.75
    k8 = (k8*10**(-9))*(avo*V)
    k9 = 0.125
    km9 = 250
    km9 = (km9*10**(-9))*(avo*V)
    km10 = 250
    km10 = (km10*10**(-9))*(avo*V)
    k10 = 0.125
    km11 = 120
    km11 = (km11*10**(-9))*(avo*V)
    km12 = 20
    km12 = (km12*10**(-9))*(avo*V)
    k12 = 5
    k12 = (k12*10**(-9))*(avo*V)
    km13 = 300
    km13 = (km13*10**(-9))*(avo*V)

    ## lists for concentrations and time

    raf = [(300*10**(-9))*avo*V]
    mek = [(300*10**(-9))*avo*V]
    erk = [(300*10**(-9))*avo*V]
    praf = [0]
    ppraf = [0]
    pmek = [0]
    ppmek = [0]
    perk = [0]
    pperk = [0]

    time = [0]

    while time[-1]<tsim:

        ## propensities computing
        a1 = ((k1*rasgtp*raf[-1]/km1)/(1+(raf[-1]/km1)+(praf[-1]/km2)))*((1+F*pperk[-1]/kf)/(1+pperk[-1]/kf))
        a2 = ((k2*rasgtp*praf[-1]/km2)/(1+ (raf[-1]/km1) + (praf[-1]/km2))) * ((1+F*pperk[-1]/kf)/(1+pperk[-1]/kf))
        a3 = (k3*ppraf[-1]/km3)/(1 + (ppraf[-1]/km3) + (praf[-1]/km4))
        a4 = (k4*praf[-1]/km3)/(1 + (ppraf[-1]/km3) + (praf[-1]/km4))
        a5 = (k5*ppraf[-1]*mek[-1]/km5)/(1 + (mek[-1]/km5) + (pmek[-1]/km6))
        a6 = (k6*ppraf[-1]*pmek[-1]/km6)/(1 + (mek[-1]/km5) + (pmek[-1]/km6))
        a7 = (k7*ppmek[-1]/km7)/(1 + (ppmek[-1]/km7) + (pmek[-1]/km8))
        a8 = (k8*pmek[-1]/km8)/(1+ (ppmek[-1]/km7) + (pmek[-1]/km8))
        a9 = (k9*ppmek[-1]*erk[-1]/km9)/(1+ (erk[-1]/km9) + (perk[-1]/km10))
        a10 = (k10*ppmek[-1]*perk[-1]/km10)/(1+ (erk[-1]/km9) + (perk[-1]/km10))
        a11 = (k11*pperk[-1]/km11)/(1+ (pperk[-1]/km11) + (perk[-1]/km12))
        a12 = (k12*perk[-1]/km12)/(1 + (pperk[-1]/km11) + (perk[-1]/km12))

        a0 = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12

        ## compute tau and add to simulation time

        tau = - 1/(a0)*np.log(random.uniform(0,1))
        
        time.append(time[-1] + tau)
        #print(time[-1])

        ## compute mu and decide which reaction happens

        mu = a0*random.uniform(0,1)
        

        if mu >= 0 and mu < a1:
            raf.append(raf[-1]-1)
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1]+1)
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1])
            pperk.append(pperk[-1])
        
        elif mu >= a1 and mu < a1+a2:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1]-1)
            ppraf.append(ppraf[-1]+1)
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1])
            pperk.append(pperk[-1])
        
        elif mu >= a1+a2 and mu < a1+a2+a3:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1]+1)
            ppraf.append(ppraf[-1]-1)
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1])
            pperk.append(pperk[-1])
        
        elif mu >= a1+a2+a3 and mu < a1+a2+a3+a4:
            raf.append(raf[-1]+1)
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1]-1)
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1])
            pperk.append(pperk[-1])

        elif mu>= a1+a2+a3+a4 and mu < a1+a2+a3+a4+a5:
            raf.append(raf[-1])
            mek.append(mek[-1]-1)
            erk.append(erk[-1])
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1]+1)
            ppmek.append(ppmek[-1])
            perk.append(perk[-1])
            pperk.append(pperk[-1])
        
        elif mu >= a1+a2+a3+a4+a5 and mu < a1+a2+a3+a4+a5+a6:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1]-1)
            ppmek.append(ppmek[-1]+1)
            perk.append(perk[-1])
            pperk.append(pperk[-1])

        elif mu >= a1+a2+a3+a4+a5+a6 and mu < a1+a2+a3+a4+a5+a6+a7:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1]+1)
            ppmek.append(ppmek[-1]-1)
            perk.append(perk[-1])
            pperk.append(pperk[-1])

        elif mu >= a1+a2+a3+a4+a5+a6+a7 and mu < a1+a2+a3+a4+a5+a6+a7+a8:
            raf.append(raf[-1])
            mek.append(mek[-1]+1)
            erk.append(erk[-1])
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1]-1)
            ppmek.append(ppmek[-1])
            perk.append(perk[-1])
            pperk.append(pperk[-1])

        elif mu >= a1+a2+a3+a4+a5+a6+a7+a8 and mu < a1+a2+a3+a4+a5+a6+a7+a8+a9:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1]-1)
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1]+1)
            pperk.append(pperk[-1])

        elif mu >= a1+a2+a3+a4+a5+a6+a7+a8+a9 and mu < a1+a2+a3+a4+a5+a6+a7+a8+a9+a10:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1]-1)
            pperk.append(pperk[-1]+1)

        elif mu >= a1+a2+a3+a4+a5+a6+a7+a8+a9+a10 and mu < a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1])
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1]+1)
            pperk.append(pperk[-1]-1)

        elif mu >= a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11 and mu < a0:
            raf.append(raf[-1])
            mek.append(mek[-1])
            erk.append(erk[-1]+1)
            praf.append(praf[-1])
            ppraf.append(ppraf[-1])
            pmek.append(pmek[-1])
            ppmek.append(ppmek[-1])
            perk.append(perk[-1]-1)
            pperk.append(pperk[-1])

    ## recalc. #-of-molecules to nM    
    raf = (np.array(raf)*10**(9))/(V*avo)
    praf = (np.array(praf)*10**(9))/(V*avo)
    ppraf = (np.array(ppraf)*10**(9))/(V*avo)
    mek = (np.array(mek)*10**(9))/(V*avo)
    pmek = (np.array(pmek)*10**(9))/(V*avo)
    ppmek = (np.array(ppmek)*10**(9))/(V*avo)
    erk = (np.array(erk)*10**(9))/(V*avo)
    perk = (np.array(perk)*10**(9))/(V*avo)
    pperk = (np.array(pperk)*10**(9))/(V*avo)



    return raf,praf,ppraf,mek,pmek,ppmek,erk,perk,pperk, time

simulation_time = 120*60 ## in seconds

#x1,x2,x3,x4,x5,x6,x7,x8,x9, time = gillespie(simulation_time)
x1,x2,x3,x4,x5,x6,x7,x8,x9, time = gillespie(7200)

print(len(x1))
print(len(time))

## from seconds to mins
time = np.array(time)/60

fig = plt.figure(1)
plt.plot(time,x1,label="raf")
plt.plot(time,x2,label="praf")
plt.plot(time,x3,label="ppraf")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend()
plt.title("case C: RAF ")
fig.savefig("caseC_raf")

fig = plt.figure(2)
plt.plot(time,x4,label="mek")
plt.plot(time,x5,label="pmek")
plt.plot(time,x6,label="ppmek")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend()
plt.title("case C: MEK ")
fig.savefig("caseC_mek")

fig = plt.figure(3)
plt.plot(time,x7,label="erk")
plt.plot(time,x8,label="perk")
plt.plot(time,x9,label="pperk")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend()
plt.title("case C: ERK ")
fig.savefig("caseC_erk")

"""
##saving the arrays
np.save("raf_C",x1)
np.save("praf_C",x2)
np.save("ppraf_C",x3)
np.save("mek_C",x4)
np.save("pmek_C", x5)
np.save("ppmek_C", x6)
np.save("erk_C", x7)
np.save("perk_C", x8)
np.save("pperk_C", x9)
np.save("time_C", time)
"""