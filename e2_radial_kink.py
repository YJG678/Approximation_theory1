import numpy as np, time
from diskpoly import *
# (iii) radial kink along the circle t = 1/2
Kmax=200
Q = disk_rule(0.0, 130, 32, t_breaks=[0,0.5,1.0])
t = Q["x"]**2+Q["y"]**2
g = np.abs(t-0.5); m0=np.sum(Q["w"]*g)
w = Q["w"]*(1+0.8*(g-m0))
print("w range", (1+0.8*(g-m0)).min(), (1+0.8*(g-m0)).max())
B = DiskBasis(Q,4); Bi = DiskBasis(Q,Kmax,mmax=3)
Ks = list(range(3,Kmax+1))
t0=time.time(); d = delta_all_K(B,3,w,Ks,Binput=Bi); print("time",time.time()-t0)
d2=d**2; ref=d2[-1]
for K in [3,4,5,6,8,10,15,20,30,40,60,80,100,150]:
    print(K, "%.10e"%d2[K-3], "gap %.3e"%(ref-d2[K-3]))
np.save("e2_radial_kink.npy", np.vstack([Ks,d2]))
