import numpy as np
import matplotlib.pyplot as plt
import seaborn

error = np.fromfile("logerror.bin")
log_h = np.fromfile("logh.bin")

plt.title('Log-log plot of the relative error')
plt.xlabel("$log_{10}(h)$")
plt.ylabel("$log_{10}(error)$")
plt.plot(log_h,error)
plt.savefig('Project1_log_error.pdf')
plt.show()
