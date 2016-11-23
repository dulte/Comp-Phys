import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy.special import gamma

folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/build-FinancialEngineering-Desktop_Qt_5_7_0_GCC_64bit-Release/"
folderSlanina = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/FinancialEngineeringSlanina/build-FinancialEngineeringSlanina-Desktop_Qt_5_7_0_GCC_64bit-Release/"
folderSen = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project5/Kode/FinancialEngineeringSen/build-FinancialEngineeringSen-Desktop_Qt_5_7_0_GCC_64bit-Release/"




# data = np.fromfile(folder + "data.bin")
# print data.shape
# meanData = np.sum(data,axis = 1)/data.shape[1]
# print meanData.shape
#


#
# b = meanData.size/10.
#
#
# plt.hist(meanData, bins = b,normed=False)
#
# plt.show()


binSen = np.loadtxt(folderSen + "binParameters.txt")
print binSen
dataSen = np.fromfile(folderSen + "data.bin")
binDataSen = np.fromfile(folderSen + "bins.bin")
xSen = np.linspace(0,binSen[-4],binDataSen.size)
print binDataSen

plt.plot(xSen,binDataSen)
plt.show()


try:
    binParaData = np.loadtxt(folder + "binParameters.txt")
    dataVar = np.loadtxt(folderSen + "variance.txt")
    x = dataVar[:,0]

    lamb = binSen[-3]#binParaData[-1]
    print binSen


    var = np.ones(x.size)*((1+3.0*lamb/(1-lamb))**(-1))


    plt.plot(x,dataVar[:,1])
    plt.plot(x,var)
    plt.title("Variance $\sigma ^2$")
    plt.xlabel("Monte Carlo Cycle")
    plt.show()
except:
    pass








# binSlanina = np.loadtxt(folderSlanina + "binParameters.txt")
# dataSlanina = np.fromfile(folderSlanina + "data.bin")
# binDataSlanina = np.fromfile(folderSlanina + "bins.bin")
# xSlanina = np.linspace(0,binSlanina[-2],binDataSlanina.size)
# print binDataSlanina
# plt.plot(xSlanina,binDataSlanina)
# plt.show()
#
# #plt.hist(np.log(dataSlanina), cumulative = 1)
# plt.hist(np.log(dataSlanina), cumulative = 1)
# plt.gca().set_xscale("log")
# #plt.gca().set_yscale("log")
# #plt.plot(dataSlanina)
# plt.show()
#
#
# binParaData = np.loadtxt(folder + "binParameters.txt")
# binData = np.fromfile(folder + "bins.bin")
# x = np.linspace(0,binParaData[-2],binData.size)
#
# x /= binParaData[2]
#
# lamb = binParaData[-1]
#
#
# n = 1 + 3.0*lamb/(1 - lamb)
# a = n**n / gamma(n)
# P = a*x**(n-1)*np.exp(-n*x)
#
#
# print np.sum(binData)/binParaData[1]/binParaData[0]
# beta = 1./binParaData[2]
#
# #analytical = 5*beta*np.exp(-beta*x)
#
#
#
# plotBin = binData/binParaData[0]/binParaData[1]*binParaData[2]/binParaData[3]
#
#
# plt.plot(x,plotBin)
# plt.plot(x,P)
# plt.legend(["Numerical","Analytical"])
# plt.show()
#
# plt.plot(x,np.log(P))
# plt.plot((x),np.log(plotBin))
# plt.show()
