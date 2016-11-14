for i in range (MCSteps):
   for k in range(N*N):
      x=int(random.uniform(0,1)*N)
      y=int(random.uniform(0,1)*N)
      dEnergy=2*grid(x,y)*(grid(x+1,y)+grid(x-1,y)+grid(x,y+1)+grid(x,y-1))
      boltzmann=boltzmann_factor(dEnergy) //Precomputed factors
      if random.uniform(0,1)=<boltzmann: //Metropolis test
          lattice(x,y)*=-1
	  energy +=dEnergy
	  magnetic_moment=2*lattice(x,y)

		
