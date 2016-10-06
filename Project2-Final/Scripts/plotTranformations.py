import numpy as np
import matplotlib.pyplot as plt
import seaborn


col = np.array([50,100,150,200,250])
transformations = np.array([3438,14284,32110,56412,87163])


plt.plot(col,transformations)
plt.xlabel("Dimensionality")
plt.ylabel("Number of transformations")
plt.show()
