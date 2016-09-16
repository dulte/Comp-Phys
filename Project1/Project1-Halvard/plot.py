import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob


import sys
if len(sys.argv) >1:
    directory = sys.argv[1]
else:
    directory = './data'



'''
try:
    filename = sys.argv[1]
except IndexError:
    print "BAD USAGE. Should take filename in build dir"
    sys.exit(1)
'''

class Plotter:
    def __init__(self):
        self.ax = plt.subplot(111)
        self.legend = []
        self.title = "N = ("

    def load_and_plot(self, filename):
        ax = self.ax
        filetype =filename.split('.')[-1] 
        if  filetype == 'txt':
            data = np.loadtxt(filename)
        elif filetype == 'bin':
            data = np.load(filename)
        else:
            print "BAD USAGE. ", filetype, " not recognized. Should take\
            bin or txt"

        x = np.linspace(0,1,len(data))
        name = filename.split('.')[-2].split('/')[-1]
        print filename
        self.legend.append(name + ", N = %d"%len(data))
        self.title += '%d ,' % len(data)
        ax.plot(x,data)

        return ax

    def show(self):
        plt.legend(self.legend)
        self.title += '\b)'
        plt.title(self.title)
        plt.show()

if __name__ == "__main__":
    files = [name for name in glob.glob(directory+"/*")]

    Plotting = Plotter()
    for filename in files:
        Plotting.load_and_plot(filename)

    Plotting.show()
