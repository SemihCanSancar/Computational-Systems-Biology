import numpy as np
import math
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint

## 1a, Deterministic Dynamics

def diff_eqns(concs,t):

    ## parameters
    alph_m = 100 * 10**(-6) ## M/min
    delta_m = 1 # 1/min

    alph_p = 10 ## 1/min
    delta_p = 0.05 ## 1/min

    ## variables
    m = concs[0]
    p = concs[1]

    ## diff. eqn

    dmdt = alph_m - delta_m*m

    dpdt = alph_p*m - delta_p*p

    return ([dmdt,dpdt])

time = np.linspace(0,100,num=300001)
conc0 = [0,0]      ##position 0: mRNA; position 1: protein 

y = odeint(diff_eqns,conc0,time)

## 1b, Stochastic Dynamics

def gillespie(concs,tsim):

    ## preliminaries
    alpha_m = 100              ## nM/min
    delta_m = 1               ## 1/min

    alpha_p = 10        ## 1/min
    delta_p = 0.05      ## 1/min

    ## arrays for concentration and time

    m_array = [concs[0]]
    p_array = [concs[1]]

    time = [0]

    while time[-1]<tsim:

        ## propensities computing

        a1 = alpha_m
        a2 = delta_m*m_array[-1]
        a3 = alpha_p*m_array[-1]
        a4 = delta_p*p_array[-1]

        a0 = a1+a2+a3+a4

        ## compute tau and add to simulation time

        tau = - 1/(a0)*np.log(random.uniform(0,1))
        
        time.append(time[-1] + 1)

        ## compute mu and decide which reaction happens

        mu = a0*random.uniform(0,1)
        #print("Mu: " + str(mu))
        #print("a0: " + str(a0))

        if mu >= 0 and mu < a1:
            m_array.append(m_array[-1] + 1)
            p_array.append(p_array[-1])
        
        elif mu >= a1 and mu < a1+a2:
            m_array.append(m_array[-1] - 1)
            p_array.append(p_array[-1])
        
        elif mu >= a1+a2 and mu < a1+a2+a3:
            p_array.append(p_array[-1] + 1)
            m_array.append(m_array[-1])
        
        elif mu >= a1+a2+a3 and mu < a0:
            p_array.append(p_array[-1] - 1)
            m_array.append(m_array[-1])
        

    return m_array, p_array, time

m,p, time = gillespie(conc0,3*10**5)
m = np.array(m)
p = np.array(p)

## Plotting the results

plt.figure(1)
plt.plot(time, y[:,0],label="Deterministic")
plt.plot(time, m/10**(6),label =r"Stochastic $1\mu m^3$")
plt.plot(time, m/10**(7),label =r"Stochastic $10\mu m^3$")
plt.legend()
plt.xlabel("Time t")
plt.ylabel("mRNA [nM]")
plt.title("Transcript conc over time")
plt.savefig("mRNA")

plt.figure(2)
plt.plot(time, y[:,1],label="Deterministic")
plt.plot(time, p/10**(6),label =r"Stochastic 1 $\mu m^3$")
plt.plot(time, p/10**(7),label =r"Stochastic 10 $\mu m^3$")
plt.legend()
plt.xlabel("Time t")
plt.ylabel("Protein [nM]")
plt.title("Protein conc over time")
plt.savefig("protein")

plt.figure(3)
plt.plot(y[:,1],y[:,0])
plt.title("Time evolution protein conc. vs mRNA conc.")
plt.xlabel("Protein conc [nM]")
plt.ylabel("mRNA conc [nM]")
plt.savefig("phase plane")