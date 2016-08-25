# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:14:08 2016

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt

err = np.fromfile("output/errors.bin")

h = [(np.log10(10**(-i))) for i in range(20)]
print err, h

plt.plot(h,err)