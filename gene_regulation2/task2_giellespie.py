import numpy as np
import math
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import norm


## 1b, Stochastic Dynamics

def gillespie(concs,tsim,V):

    ## constant

    avo = 6.022*10**(23)       ## 1/mol

    ## preliminaries
    alpha_m = 100*100*10**(-9)    ## M/min
    delta_m = 1               ## 1/min

    alpha_p = 10/100        ## 1/min
    delta_p = 0.05      ## 1/min

    ## arrays for concentration and time

    m_array = [concs[0]]
    p_array = [concs[1]]

    time = [0]

    while time[-1]<tsim:

        ## propensities computing

        a1 = alpha_m*V*avo
        a2 = delta_m*m_array[-1]
        a3 = alpha_p*m_array[-1]
        a4 = delta_p*p_array[-1]

        a0 = a1+a2+a3+a4

        ## compute tau and add to simulation time

        tau = - 1/(a0)*np.log(random.uniform(0,1))
        
        time.append(time[-1] + tau)

        ## compute mu and decide which reaction happens

        mu = a0*random.uniform(0,1)
        

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
        
    m_array = np.array(m_array)/(V*avo)

    p_array = np.array(p_array)/(V*avo)


    return m_array[-100000:], p_array[-100000:], time[-100000:]
    #return m_array, p_array, time

V1 = 1*10**(-15)    ## L
V2 = 10*10**(-15)   ## L 

conc0 = [0,0]

m,p, time = gillespie(conc0,100,V2)
"""
print(m.shape)


plt.plot(p)
plt.savefig("bild2")

"""