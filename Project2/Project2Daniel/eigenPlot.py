import numpy as np
import matplotlib.pyplot as plt


rho = np.linspace(0,10,100)

u = np.loadtxt("eigenvec.txt")


plt.plot(rho,u)
plt.show()
