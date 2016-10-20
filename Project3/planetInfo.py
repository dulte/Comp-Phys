import telnetlib
import numpy as np
import sys



class planetInfo:

    def __init__(self):
        pass


    def printInfoToFile(self,startdate,enddate,IDs=[],filename = ""):



        if len(IDs) == 0 and len(filename) >0:
            IDs = np.loadtxt(filename)
        elif len(IDs) == 0 and len(filename) == 0:
            print "No IDs or filename given"
            sys.exit(1)
        #print IDs
        print "Getting info for %g planet(s)" %len(IDs)

        data = np.zeros((len(IDs),2,3))
        names = []
        IDs = np.array(IDs)
        for i in range(IDs.size):
            name,data[i,:,:] = self.getInfo(int(IDs[i]),startdate,enddate)
            names.append(name)
            print "Found info for %s" %name

        np.savetxt("planetPositions.txt",data[:,0,:].reshape(len(names),3))
        np.savetxt("planetVelocities.txt",data[:,1,:].reshape(len(names),3))

        with open("planetConfig_" + str(len(names)) + ".txt","w") as outFile:

            outFile.write(str(len(names)) + "\n")
            for name in names:
                outFile.write(name + "\n")



    # def getIDsFromList(self,filename):
    #     return np.loadtxt(filename)


    def getInfo(self,ID,startdate,enddate):

        self.t = telnetlib.Telnet()
        self.t.open('horizons.jpl.nasa.gov', 6775)

        expect = ( ( r'Horizons>', str(ID) + '\n' ),
           ( r'< Scroll.*:', 'q' ),
           ( r'Select.*E.phemeris.*:', 'E\n'),
           ( r'Observe.*:', 'v\n' ),
           ( r'Coordinate center.*:', '500@0\n' ),
           ( r'Confirm selected station.*>', 'y\n'),
           ( r'Reference plane.*:', 'eclip\n'),
           ( r'Starting .* :', str(startdate)+'\n' ),
           ( r'Ending .* :', str(enddate)+'\n' ),
           ( r'Output interval.*:', '1d\n' ),
           ( r'Accept default output.* :', 'y\n' ),
           ( r'Scroll . Page: .*%', ' '),
           ( r'Select\.\.\. .A.gain.* :', 'x\n' ))




        with open('temp.txt', 'w') as fp:
            while True:
                try:
                    answer = self.t.expect(list(i[0] for i in expect), 10)
                except EOFError:
                    break


                fp.write(answer[2])
                fp.flush()
                self.t.write(expect[answer[0]][1])

        data = np.zeros((2,3))
        name = ""
        exponant = 0
        mass = 0
        foundMass = False
        with open('temp.txt', 'r') as fp:
            for line in fp:
                words = line.split()
                if len(words) == 0:
                    continue

                if words[0] == "Target":
                    name = words[3]
                    continue
                if (("Mass" in words) or ("Mass," in words)) and (not foundMass):# or ("mass" in words) or ("mass," in words):
                    for i in range(len(words)):
                        if (words[i] == "Mass") or (words[i] == "Mass,"):

                            try:
                                if len(words[i+1]) == 6:
                                    exponant = float(words[i+1][4:])
                                else:

                                    exponant = float(words[i+1][3:])

                                if (words[i+2] == "g") or (words[i+2] == "g)"):
                                    exponant -= 3
                            except:
                                if len(words[i+2]) == 6:
                                    exponant = float(words[i+2][4:])
                                else:

                                    exponant = float(words[i+2][3:])

                                if (words[i+2] == "g") or (words[i+2] == "g)"):
                                    exponant -= 3
                            for j in range(len(words)-i):
                                try:
                                    mass = float(words[j+i])*10**exponant
                                    break
                                except:
                                    try:
                                        nwords = words[j+i].split("+")

                                        mass = float(nwords[0])*10**exponant
                                        break
                                    except:
                                        continue
                            foundMass = True
                            break

                if words[0] == "$$SOE":
                    fp.next()
                    data[0,:] = np.array(fp.next().split())
                    data[1,:] = np.array(fp.next().split())
                    data[1,:] *= 365.24
                    break

        print "Mass: ",mass
        return name,data


info = planetInfo()

start = "2016-Oct-18"
end = "2016-Oct-19"

IDs = [399,801]
info.printInfoToFile(start,end,filename = "targetids.txt")
