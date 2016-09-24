for(i = 2; i< N+1):
    f[i] += f[i-1]*a_inv[i-1];

u[N] = f[N]*a_inv[N];

for(i = N-1; i > 0):
    u[i] = (f[i] + u[i+1])*a_inv[i];
