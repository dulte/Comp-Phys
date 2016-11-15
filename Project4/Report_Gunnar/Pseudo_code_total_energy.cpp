for i in range(Nspins):
   for j in range(Nspins):
      energy -= spinMatrix(i,j)*(spinMatrix(i-1,j)+ spinMatrix(i,j-1));	
