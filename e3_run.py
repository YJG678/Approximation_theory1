import numpy as np, time, json
from numpy.polynomial import legendre as Lg
from diskpoly import *
def hp(t,p):
    c=np.zeros(p+1); c[p]=1; return Lg.legval(2*t-1,c)
def W1_iso(p, eps):
    r = np.linspace(0,1,400001); s=r**2
    c=np.zeros(p+1); c[p]=1; ci=Lg.legint(c, lbnd=-1)
    H = Lg.legval(2*s-1, ci)/2
    return 2*eps*np.trapezoid(np.abs(H), r)       # factor 2: isotropic scaling (radius 2)
def delta3(p, eps):
    Kmax = 2*p+24
    Q = disk_rule(0.0, Kmax+p+12, 32)
    t = Q["x"]**2+Q["y"]**2
    w = Q["w"]*(1+eps*hp(t,p))
    B = DiskBasis(Q,4); Bi = DiskBasis(Q,Kmax,mmax=3)
    d = delta_all_K(B,3,w,[Kmax-8,Kmax],Binput=Bi)
    return d[-1], abs(d[-1]-d[0])/d[-1]
rows=[]
for p in [3,4,5,6,8,10,12,16,20,24,32,40]:
    for eps in [0.1]:
        D,conv = delta3(p,eps); rows.append((p,eps,D,W1_iso(p,eps),conv))
for p in [3,12]:
    for eps in [1e-3,3e-3,1e-2,3e-2,0.3,0.5]:
        D,conv = delta3(p,eps); rows.append((p,eps,D,W1_iso(p,eps),conv))
rows=np.array(rows)
np.save("e3_radial.npy", rows)
for p_, eps_, D_, W_, conv_ in rows:
    print("p=%2d eps=%.3g Delta=%.4e W1=%.4e W1/Delta=%.3f conv=%.1e" % (p_, eps_, D_, W_, W_ / D_, conv_))
