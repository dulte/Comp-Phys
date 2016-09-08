import numpy as np
import matplotlib.pyplot as plt


u = np.fromfile("numericalSolution.bin")
error = np.fromfile("logerror.bin")

x = np.linspace(0,1,1002)
print x.size, u.size

log_h = [-(1+i/4.) for i in range(error.size)]

exact = 1 - (1-np.exp(-10))*x - np.exp(-10*x)



plt.plot(x,u)
plt.plot(x,exact,"r--")
plt.legend(["Numerical", "Exact"])

plt.show()

plt.xlabel("$log_{10}(h)$")
plt.ylabel("$log_{10}(error)$")
plt.plot(log_h,error)
plt.show()
