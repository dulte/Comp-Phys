import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


import sys
'''
try:
    filename = sys.argv[1]
except IndexError:
    print "BAD USAGE. Should take filename in build dir"
    sys.exit(1)
'''


directory = '../build-Project1Halvard-Desktop-Debug/%s' 
arma_name = "ArmaSolution.txt"
exact_name = "ComputedSolution.txt"
comp_name = "ExactSolution.txt"
lu_decomp = np.loadtxt(directory % arma_name)
exact = np.loadtxt(directory % exact_name)
computed = np.loadtxt(directory % comp_name)


x1 = np.linspace(0,1,len(exact))
x2 = np.linspace(0,1,len(computed))
x3 = np.linspace(0,1,len(lu_decomp))
plt.plot(x1, exact, '--')
plt.plot(x2, computed)
plt.plot(x3, lu_decomp)
plt.legend(["exact","computed", "lu"])

plt.title('N computed = %d, N lu decomposition = %d'% (len(computed),
    len(lu_decomp)))
plt.show()
