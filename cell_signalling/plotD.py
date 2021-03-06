import numpy as np
import matplotlib.pyplot as plt

x1 = np.load("raf_D.npy")
x2 = np.load("praf_D.npy")
x3 = np.load("ppraf_D.npy")

x4 = np.load("mek_D.npy")
x5 = np.load("pmek_D.npy")
x6 = np.load("ppmek_D.npy")

x7 = np.load("erk_D.npy")
x8 = np.load("perk_D.npy")
x9 = np.load("pperk_D.npy")

time = np.load("time_D.npy")

plt.figure(1)
plt.plot(time,x1,label="raf")
plt.plot(time,x2,label="praf")
plt.plot(time,x3,label="ppraf")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend(loc="best")
plt.title("case D: RAF ")
plt.savefig("caseD_raf")

plt.figure(2)
plt.plot(time,x4,label="mek")
plt.plot(time,x5,label="pmek")
plt.plot(time,x6,label="ppmek")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend(loc="best")
plt.title("case D: MEK ")
plt.savefig("caseD_mek")

plt.figure(3)
plt.plot(time,x7,label="erk")
plt.plot(time,x8,label="perk")
plt.plot(time,x9,label="pperk")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend(loc="best")
plt.title("case D: ERK ")
plt.savefig("caseD_erk")

plt.figure(4)
plt.plot(time,x3,label="ppraf")
plt.plot(time,x6,label="ppmek")
plt.plot(time,x9,label="pperk")
plt.xlabel("time [min]")
plt.ylabel("concentration [nM]")
plt.legend(loc="best")
plt.title("case D: Double phosphated ")
plt.savefig("caseD_pp")
