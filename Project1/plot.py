"""
File: plot.py
Author: Halvard Sutterud, Daniel, Gunnar
Email: halvard.sutterud@gmail.com, + div
Github: https://github.com/halvarsu, /dulte, /gunnarius
Description: Plots all files in a directory into file in figures with name
of directory. Writes new files for each run, unless "clean"

USAGE: python plot.py [-h] [-dir DIRECTORY] [-c]
-c cleans the files with same name before saving file
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import argparse

'''
try:
    filename = sys.argv[1]
except IndexError:
    print "BAD USAGE. Should take filename in build dir"
    sys.exit(1)
'''

class Plotter:
    def __init__(self, directory):
        self.ax = plt.subplot(111)
        self.legend = []
        self.directory = directory

    def build_plot(self, filename):
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
        print 'plotting file ', filename
        self.legend.append(name + ", N = %d"%len(data))
        ax.plot(x,data)

        return ax

    def plot(self, figname, clean = False, save = True, maketex = True):
        plt.legend(self.legend)
        plt.title(self.directory)
        plt.show()
        fig_dir = "./figure/"
        fig_filename  = fig_dir + directory 

        if clean:
            import os
            old_files = glob.glob(fig_filename + '*.png')
            for old_file in old_files:
                print "Removing file     : ", old_file
                os.remove(old_file)

        if save:
            num = len(glob.glob(fig_filename + '*')) + 1
            fig_filename += str(num)
            print "Saving fig to file: ", fig_filename + '.png'
            plt.savefig(fig_filename)

            
            '''
            if maketex:
                tex_dir = "./tex_test/"
                tex_filename  = tex_dir + directory 
                num = len(glob.glob(tex_filename + '*')) + 1

                plt.savefig(tex_filename + str(num))
                num = len(glob.glob(tex_file + '*')) + 1
                tex_file = open(tex_name+str(num), 'w')
                tex_file.write(tex_figure_text)
                '''
           

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", "--directory", help="Directory to plot",
             default = 'data')
    parser.add_argument("-c", "--clean", help="Clean the figures with same\
            name", action = "store_true", default = False)
    parser.add_argument("-ds", "--dont_save", help="Dont save", action =
            "store_true", default = False)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    directory = args.directory
    files = [name for name in glob.glob(directory+"/*")]

    Plotting = Plotter(directory)
    for filename in files:
        Plotting.build_plot(filename)
    Plotting.plot(directory, save = not args.dont_save, clean = args.clean)
