num=zeros(shape=bins)
m0=m0 //Initial amount of money
N=N //Initial number of players
m=zeros(N)
for i in range(no_of_experiments):
	m.fill(m0)
	for k in range(no_of_transactions):
		i=random.uniformint(0, N)
		j=random.uniformint(0, N)
		eps=random.uniform(0,1)
		money_i=m[i]
		money_j=m[k]
		m[i]=eps*(money_i+money_k)
		m[k]=(1-eps)*(money_i+money_k)
	for bin in bins:
		for agent in m:
			if agent in bin:
				num[bin]+=1
