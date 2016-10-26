import numpy as np
import matplotlib.pyplot as plt


folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project3/build-Project3-Desktop_Qt_5_7_0_GCC_64bit-Release/"
data = np.loadtxt(folder + "errorEnergyMomentum.txt")

print data.shape
errorMomentum = (data[:,0])
errorEnergy = (data[:,1])

#x = np.linspace(0,np.log10(errorEnergy.size),errorEnergy.size)

plt.plot(errorEnergy)
plt.show()
plt.plot(errorMomentum)
plt.show()
