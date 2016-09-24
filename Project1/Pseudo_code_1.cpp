for(i = 2; i < N+1):
    d = c[i-1]/a[i-1];
    a[i] -= b[i-1]*d;
    f[i] -= f[i-1]*d;

u[N] = f[N]/a[N];

for(i = N-1; i >0)_
    u[i] = (f[i] - b[i]*u[i+1])/a[i]
