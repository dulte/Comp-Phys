import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


exact = np.loadtxt('../build-Project1Halvard-Desktop-Debug/ExactSolution.txt')
computed = np.loadtxt('../build-Project1Halvard-Desktop-Debug/ComputedSolution.txt')


x1 = np.linspace(0,1,len(exact))
x2 = np.linspace(0,1,len(computed))
plt.plot(x1, exact)
plt.plot(x2, computed)

plt.title('N = %d'% len(computed))
plt.show()
