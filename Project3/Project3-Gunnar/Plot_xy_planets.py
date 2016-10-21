import numpy as np
import matplotlib.pyplot as plt
import time, glob, os



with open("planetData.txt","r") as infile:
    bodies = []
    for line in infile:
        words = line.split()
        bodies.append(words[0])

    #bodies=['Sun','Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']


no_frames=1

#folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project3/build-Project3-Desktop_Qt_5_7_0_GCC_64bit-Release/"
infile=open(folder + 'positionsEarthSunEuler.xyz')
lines=np.array(infile.readlines())
no_bodies= int(lines[0])
position=np.zeros(shape=(no_bodies, 2))
lines_true= np.logical_and(lines != str(no_bodies)+'\n', lines != 'New Time step\n')
lines=lines[lines_true]
position=np.zeros(shape=(int(len(lines)/float(no_bodies)),no_bodies, 2))


for i in range(0,np.shape(lines)[0], no_bodies):
	for k in range(no_bodies):
		position[int(i/float(no_bodies)),k,0]=float(((lines[i+k]).split())[0])
		position[int(i/float(no_bodies)),k,1]=float(((lines[i+k]).split())[1])

legend=[]

for j in range(no_bodies):
	legend.append(bodies[j])



old_files = glob.glob('solar_system*.png')
for file in old_files:
	os.remove(file)


frames=np.linspace(0, int(np.shape(lines)[0]/float(no_bodies))-1, no_frames)
colors=plt.cm.rainbow(np.linspace(0,1,no_bodies))

frame_numbering=np.linspace(0, no_frames, no_frames)


for k in range(len(frames)):
    for j in range(no_bodies):
        plt.plot(position[:, j, 0], position[:, j,1], c=colors[j])
        plt.xlim([-40, 60])
        plt.ylim([-40, 50])
        plt.legend(legend, loc='upper left', prop={'size':8},handles="line")
        plt.xlabel('x [AU]')
        plt.ylabel('y [AU]')
        plt.title('2D Projection of Solar System movements')
        plt.grid('on')
        plt.axis('equal')
        plt.plot(position[int(frames[k]), j, 0], position[int(frames[k]),j,1], marker ='*')
    print "Done with:", k
    plt.show()
    #plt.savefig('solar_system%06d.png' %frame_numbering[k])
    plt.clf()
