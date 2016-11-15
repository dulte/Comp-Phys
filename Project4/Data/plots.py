import numpy as np
import collections as cl
import matplotlib.pyplot as plt
import seaborn

class Lattice_2_by_2:

	def __init__(self, J, k_B,T):
		self.J=J
		self.k_B=k_B
		self.T=T
		self.beta=1.0/(k_B*self.T)
		self.beta2=1.0/(k_B*self.T*self.T)
		self.arg=8.0*J*self.beta

	def compute_partition_function(self):
		partition_func=4.0*np.cosh(self.arg)+12
		return partition_func

	def compute_energy(self):
		first_fac=-32.0*self.J*np.sinh(self.arg)
		second_fac=4.0*np.cosh(self.arg)+12
		energy=first_fac/second_fac
		return energy

	def compute_energy_squared(self):
		first_fac=256.0*self.J*self.J*np.cosh(self.arg)
		second_fac=4.0*np.cosh(self.arg)+12.0
		energy_squared=first_fac/second_fac

		return energy_squared

	def compute_magnetization(self):
		return 0

	def compute_abs_magnetization(self):
		abs_magnetization_num=2.0*np.exp(self.arg)+4
		abs_magnetization_denom=np.cosh(self.arg)+3
		return abs_magnetization_num/abs_magnetization_denom

	def compute_magnetization_squared(self):
		first_fac=8.0*np.exp(self.arg)+4.0
		second_fac=np.cosh(self.arg)+3

		magnetization_squared=first_fac/second_fac
		return magnetization_squared	

	def compute_CV(self):
		first_fac=128.0*(self.J**2)*(np.exp(self.arg)+np.exp(-self.arg))
		second_fac=-16.0*self.J*np.exp(self.arg)+16.0*self.J*np.exp(-self.arg)
		denominator=2.0*np.exp(self.arg)+2.0*np.exp(-self.arg)+12.0
	
		C_V=self.beta2*((first_fac/denominator)-(second_fac/denominator)**2)
		return C_V

	def compute_CV_hyp(self):
		first_fac=256.0*self.J*self.J*np.cosh(self.arg)
		second_fac=-32.0*self.J*np.sinh(self.arg)
		denominator=4.0*np.cosh(self.arg)+12.0

		C_V=self.beta2*((first_fac/denominator)-(second_fac/denominator)**2)
		return C_V

	def compute_MagSus(self):
		first_fac=32.0*np.exp(self.arg)+32
		denominator=2.0*np.exp(self.arg)+2.0*np.exp(-self.arg)+12.0

		MagSus=self.beta*(first_fac/denominator)
		return MagSus

	def compute_MagSus_hyp(self):
		first_fac=8.0*np.exp(self.arg)+8.0
		denominator=np.cosh(self.arg)+3.0

		MagSus=self.beta*(first_fac/denominator)
		return MagSus

folder='/home/gunnnar/Documents/FYS3150/Comp-Phys/Project4/Data/'

data = np.loadtxt(folder +"N60Fabs.txt")
lattice=Lattice_2_by_2(1,1, data[:,0])

L = 60

plt.title(r"Heat capacity per spin, $C_v$ for L=%i" %L)
plt.xlabel("Temperature [K]")
plt.ylabel("$C_v$ [J/K]")
plt.plot(data[:,0],data[:,3])
plt.plot(data[np.argmax(data[:,3]),0],np.max(data[:,3]),"ro")
plt.text(data[np.argmax(data[:,3]),0]-0.01,np.max(data[:,3])+0.1,"$T_c = %g$" %data[np.argmax(data[:,3]),0] )
plt.show()
plt.title(r"Magnetic susceptibility per spin, $\chi$ for L=%i" %L)
plt.xlabel("Temperature [K]")
plt.ylabel("Magnetic susceptibility [$Js^2$]")
plt.plot(data[:,0],data[:,-1])
plt.plot(data[np.argmax(data[:,-1]),0],np.max(data[:,-1]),"ro")
plt.text(data[np.argmax(data[:,-1]),0]-0.05,np.max(data[:,-1])-2,"$T_c = %g$" %data[np.argmax(data[:,-1]),0] )
plt.show()
plt.plot(data[:,0],data[:,1])
plt.xlabel("Temperature [K]")
plt.ylabel("Mean energy [J]")
plt.title(r"Mean Energy, $\langle E \rangle$ per spin for L=%i" %L)
mean_energy=lattice.compute_energy()
plt.show()
plt.plot(data[:,0],data[:,-2]/float(L*L))
plt.ylabel("Magnetic Moment [J/T]")
plt.xlabel("Temperature [K]")
plt.title(r"Mean absolute Magnetic Moment $\langle |M| \rangle$ per spin for L=%i" %L)
mean_magnetic_moment=lattice.compute_abs_magnetization()
plt.show()







#folder2 = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project4/Kode/build-Project4c-Desktop_Qt_5_7_0_GCC_64bit-Release/"
#dataOneTemp = np.loadtxt(folder2 + "dataOneTemp.txt")

# plt.title("Energy")
# plt.xlabel("Monte Carlo Step")
# plt.ylabel("Energy [J]")
# plt.plot(dataOneTemp[:,0],dataOneTemp[:,1])
# plt.show()
# plt.title("Mean Energy")
# plt.xlabel("Monte Carlo Step")
# plt.ylabel("Energy [J]")
# plt.plot((dataOneTemp[:,0]),dataOneTemp[:,2])
# plt.show()
# plt.title("Mean Magnetic moment")
# plt.xlabel("Monte Carlo Step")
# plt.ylabel("Magnetic Moment [J/T]")
# plt.plot(dataOneTemp[:,0],dataOneTemp[:,-3])
# plt.show()
# plt.title("Accepted Flips")
# plt.xlabel(["Monte Carlo Step"])
# plt.plot((dataOneTemp[:,0]),dataOneTemp[:,-2]/dataOneTemp[:,0])
# plt.show()
#
# plt.title("Variance Of Enery")
# plt.xlabel("Monte Carlo Step")
# plt.plot((dataOneTemp[:,0]),dataOneTemp[:,3])
# plt.show()


# #
# uniqueE,countsE = np.unique(dataOneTemp[:,1], return_counts = True)
# plt.plot(uniqueE,countsE)
# print uniqueE,countsE
# #plt.hist(uniqueE,countsE)
# plt.show()
