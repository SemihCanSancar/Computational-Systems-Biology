import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.stats import norm


mrna_a1 = np.load("mrna_V1_partA.npy",allow_pickle=True).astype("float64")
mrna_a2 = np.load("mrna_V2_partA.npy",allow_pickle=True).astype("float64")
mrna_b1 = np.load("mrna_V1_partB.npy",allow_pickle=True).astype("float64")
mrna_b2 = np.load("mrna_V2_partB.npy",allow_pickle=True).astype("float64")

protein_a1 = np.load("protein_V1_partA.npy",allow_pickle=True).astype("float64")
protein_a2 = np.load("protein_V2_partA.npy",allow_pickle=True).astype("float64")
protein_b1 = np.load("protein_V1_partB.npy",allow_pickle=True).astype("float64")
protein_b2 = np.load("protein_V2_partB.npy",allow_pickle=True).astype("float64")


"""
## stationary state plot for each case 
x = np.linspace(0,250,100000)

plt.figure(1)
plt.plot(x,mrna_a1[:100000],color="b",label=r"V = 1 $\mu m^3$")
plt.plot(x,mrna_a2[:100000],color="k",label=r"V = 10 $\mu m^3$")
plt.title("mRNA concentration default case in stationary state")
plt.xlabel("simulation time")
plt.ylabel("concentration [nM]")
plt.legend()
plt.savefig("mrna_A")

plt.figure(2)
plt.plot(x,mrna_b1[:100000],color="b",label=r"V = 1 $\mu m^3$")
plt.plot(x,mrna_b2[:100000],color="k",label=r"V = 10 $\mu m^3$")
plt.title("mRNA concentration case I in stationary state")
plt.xlabel("simulation time")
plt.ylabel("concentration [nM]")
plt.legend()
plt.savefig("mrna_B")

plt.figure(3)
plt.plot(x,protein_a1[:100000],color="b",label=r"V = 1 $\mu m^3$")
plt.plot(x,protein_a2[:100000],color="k",label=r"V = 10 $\mu m^3$")
plt.title("protein concentration default case in stationary state")
plt.xlabel("simulation time")
plt.ylabel("concentration [nM]")
plt.legend()
plt.savefig("protein_A")

plt.figure(4)
plt.plot(x,protein_b1[:100000],color="b",label=r"V = 1 $\mu m^3$")
plt.plot(x,protein_b2[:100000],color="k",label=r"V = 10 $\mu m^3$")
plt.title("protein concentration case I in stationary state")
plt.xlabel("simulation time")
plt.ylabel("concentration [nM]")
plt.legend()
plt.savefig("protein_B")

"""
# Fit a normal distribution to the data:
mu, std = norm.fit(mrna_a1)

plt.figure(1)

# Plot the histogram.
plt.hist(mrna_a1, bins=35, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"mRNA in V = 1 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_mrna_a11")

# Fit a normal distribution to the data:
mu, std = norm.fit(mrna_a2)

plt.figure(2)

# Plot the histogram.
plt.hist(mrna_a2, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"mRNA in V = 10 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_mrna_a2")

# Fit a normal distribution to the data:
mu, std = norm.fit(mrna_b1)

plt.figure(3)

# Plot the histogram.
plt.hist(mrna_b1, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"mRNA in V = 1 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_mrna_b1")

# Fit a normal distribution to the data:
mu, std = norm.fit(mrna_b2)

plt.figure(4)

# Plot the histogram.
plt.hist(mrna_b2, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"mRNA in V = 10 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_mrna_b2")

plt.figure(5)
# Fit a normal distribution to the data:
mu, std = norm.fit(protein_a1)

# Plot the histogram.
plt.hist(protein_a1, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"protein in V = 1 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_protein_a1")

plt.figure(6)
# Fit a normal distribution to the data:
mu, std = norm.fit(protein_a2)

# Plot the histogram.
plt.hist(protein_a2, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"protein in V = 10 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_protein_a2")

plt.figure(7)
# Fit a normal distribution to the data:
mu, std = norm.fit(protein_b1)

# Plot the histogram.
plt.hist(protein_b1, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"protein in V = 1 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_protein_b1")

plt.figure(8)
# Fit a normal distribution to the data:
mu, std = norm.fit(protein_b2)

# Plot the histogram.
plt.hist(protein_b2, bins=100, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = r"protein in V = 10 $\mu$M; mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.xlabel("concentration [nM]")
plt.ylabel("Density")

plt.savefig("erg_protein_b2")
