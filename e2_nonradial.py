import numpy as np, time, sys
from diskpoly import *
Kmax=int(sys.argv[1]) if len(sys.argv)>1 else 56   # K_ref = 56 in the paper
sect=[0,np.pi/2,np.pi,1.5*np.pi,2*np.pi]
Q = disk_rule(0.0, Kmax+16, Kmax+16, sectors=sect, radial_in_r=True)
x,y=Q["x"],Q["y"]; print("nodes",x.size)
t0=time.time(); B = DiskBasis(Q,Kmax); print("basis",time.time()-t0, B.val.shape)
z=x+1j*y
g=np.abs(x)+np.abs(y); m0=np.sum(Q["w"]*g)
dens = {"smooth": 1+0.1*np.real(z**5), "kink": 1+0.5*(g-m0)}
print("kink w range", dens["kink"].min(), dens["kink"].max())
Ks=list(range(3,Kmax+1))
out={}
for name,dd in dens.items():
    for N in (2,3):
        t0=time.time(); d=delta_all_K(B,N,Q["w"]*dd,Ks); out[(name,N)]=d**2
        print(name,N,"time %.1f"%(time.time()-t0), " ".join("K%d:%.8e"%(K,d[K-3]**2) for K in (5,10,20,30,40,Kmax)))
np.savez("e2_nonradial_K%d.npz"%Kmax, Ks=Ks, **{f"{a}_{b}":v for (a,b),v in out.items()})
