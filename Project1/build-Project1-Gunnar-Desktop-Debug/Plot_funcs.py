import numpy as np
import matplotlib.pyplot as plt
def exact_func(x):
	return 1-(1-np.exp(-10))*x-np.exp(-10*x)
infile=np.loadtxt("num.txt")
np.insert(infile, 0, 0)
np.append(infile, 0)
print infile
N=len(infile)
print N
infile_2=np.fromfile("numericalSolution.bin")
x=np.linspace(0,1,N)
x_2=np.linspace(0,1,len(infile_2))
exact=exact_func(x)
plt.plot(x, infile)
plt.hold("on")
plt.plot(x, exact)
plt.hold("on")
plt.plot(x_2, infile_2)
plt.legend(['LU', 'Exact', 'Gauss'])
plt.show()
