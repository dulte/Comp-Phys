import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy.special import gamma

#folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/build-FinancialEngineering-Desktop_Qt_5_7_0_GCC_64bit-Release/"
#folderSlanina = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/FinancialEngineeringSlanina/build-FinancialEngineeringSlanina-Desktop_Qt_5_7_0_GCC_64bit-Release/"
#folderSen = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/FinancialEngineeringSen/build-FinancialEngineeringSen-Desktop_Qt_5_7_0_GCC_64bit-Release/"






binSen = np.loadtxt(folderSen + "binParameters.txt")
binParaData = np.loadtxt(folderSen + "binParameters.txt")
dataSen = np.fromfile(folderSen + "data.bin")
binDataSen = np.fromfile(folderSen + "bins.bin")
xSen = np.linspace(0,binSen[-4],binDataSen.size)

xSen /= binParaData[2]

lamb = binParaData[-3]


n = 1 + 3.0*lamb/(1 - lamb)
a = n**n / gamma(n)
P = a*xSen**(n-1)*np.exp(-n*xSen)

powerLaw = xSen**(-1-0.9)/binParaData[0]/binParaData[1]*binParaData[2]/binParaData[3]

plotBin = binDataSen/binParaData[0]/binParaData[1]*binParaData[2]/binParaData[3]

alpha = binParaData[-2]
gamma = binParaData[-1]


#plt.plot(xSen,powerLaw)
plt.title(r"Distribution of Wealth $F(m/\langle m \rangle)$ for $\lambda$ = %g, $\alpha = %g$,$\gamma = %g$ " %(lamb,alpha,gamma))
plt.xlabel(r"$m/\langle m \rangle$")
plt.ylabel(r"$F(m/\langle m \rangle)$")
plt.plot(xSen,plotBin)
#plt.plot(xSen,P)
plt.legend(["Numerical", "Analytical"])
plt.show()

plt.title(r"Log of distribution of Wealth $F(m/\langle m \rangle)$ for $\lambda$ = %g , $\alpha = %g$,$\gamma = %g$ " %(lamb,alpha,gamma))
plt.xlabel(r"$m/\langle m \rangle$")
plt.ylabel(r"$log(F(m/\langle m \rangle))$")
plt.loglog(xSen,(plotBin))
#plt.plot(xSen,np.log(P))
plt.legend(["Numerical", "Analytical"])
plt.show()

try:
    errSen = np.loadtxt(folderSen + "distError.txt")
    plt.plot(errSen[:,0],errSen[:,1])
    plt.title(r"Relative difference for $\lambda$ = %g, $\alpha = %g$,$\gamma = %g$ " %(lamb,alpha,gamma))
    plt.xlabel("Monte Carlo Steps")
    plt.ylabel("Relative difference")
    plt.show()
except:
    pass


try:

    dataVar = np.loadtxt(folderSen + "variance.txt")
    x = dataVar[:,0]

    lamb = binSen[-3]#binParaData[-1]



    var = np.ones(x.size)*((1+3.0*lamb/(1-lamb))**(-1))


    plt.plot(x,dataVar[:,1])
    #plt.plot(x,var)
    plt.title(r"Variance $\sigma ^2$ for $\lambda = %g$, $\alpha = %g$,$\gamma = %g$ " %(lamb,alpha,gamma))
    plt.xlabel("Transactions")
    plt.ylabel("$\sigma ^2$")
    #plt.legend(["Numerical","Analytical"])
    plt.show("$\sigma ^2$")
    # plt.plot(x,dataVar[:,-1])
    # plt.show()
except:
    pass
