import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy.special import gamma
from scipy.optimize import curve_fit


def powerLaw(x,a,b):

    return a*x**(-1-b)






#folderSen = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/FinancialEngineeringSen/build-FinancialEngineeringSen-Desktop_Qt_5_7_0_GCC_64bit-Release/"
#folderSen = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/data/lamb9/"

#folderSen = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/data/alpha2/gamma2/"
binSen = np.loadtxt(folderSen + "binParameters.txt")
binParaData = np.loadtxt(folderSen + "binParameters.txt")
dataSen = np.fromfile(folderSen + "data.bin")
binDataSen = np.fromfile(folderSen + "bins.bin")
xSen = np.linspace(0,binSen[-4],binDataSen.size)

xSen /= binParaData[2]

lamb = binParaData[-3]


plotBin = binDataSen/binParaData[0]/binParaData[1]*binParaData[2]/binParaData[3]



#lamb = 0.0

xMax = 2./np.sqrt(lamb+0.1) -.1
#xSen = np.linspace(xMax,xMax+1,1000)


n = 1 + 3.0*lamb/(1 - lamb)
a = n**n / gamma(n)
P = a*xSen**(n-1)*np.exp(-n*xSen)


guess = 1e9,1.
#powerLaw = xSen**(-1-2)
popt, pcov = curve_fit(powerLaw, xSen[xSen>=xMax], plotBin[xSen>=xMax],guess)#,bounds=([0.,1.], [10000.,2.]))
print popt

b = findCurve(powerLaw,xSen,P)

alpha = binParaData[-2]
gamma = binParaData[-1]

errSq = np.sum((powerLaw(xSen,popt[0],popt[1])[xSen>=xMax] - plotBin[xSen>=xMax])**2)
print "Error: ", errSq


plt.title(r"Log of tail and $m^{-1-\beta}$ for $\lambda = %g$, $\alpha = %g$,$\gamma = %g$ " %(lamb,alpha,gamma))
plt.xlabel(r"$m/\langle m \rangle$")
plt.ylabel(r"$log(F(m/\langle m \rangle))$")
plt.semilogy(xSen[xSen>=xMax],powerLaw(xSen,popt[0],popt[1])[xSen>=xMax])
#plt.semilogy(xSen,powerLaw(xSen,1,b))
#plt.semilogy(xSen,powerLaw)
plt.semilogy(xSen[xSen>=xMax],P[xSen>=xMax])
plt.semilogy(xSen[xSen>=xMax],plotBin[xSen>=xMax])

plt.legend(["Power law  $%.1f \cdot m^{-1-%.2f}$" %(popt[0],popt[1]),"Tail, Analytical","Tail, Numerical"])
plt.show()
