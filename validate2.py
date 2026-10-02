import numpy as np, time
from diskpoly import *
# Table 1 measure (unit disk, similarity invariant)
Q = disk_rule(0.0, 40, 96); B = DiskBasis(Q, 12)
z = Q["x"] + 1j*Q["y"]; w = Q["w"]*(1+0.1*np.real(z**5))
print("table N=2:", delta_all_K(B,2,w,[5,6,7])**2, "  N=3:", delta_all_K(B,3,w,[5,6,7])**2)
print("  paper: 6.944444e-4 6.944444e-4 9.424603e-4 | 6.944444e-4 7.931806e-4 9.424603e-4")
# Gegenbauer compatibility at high degree
for al in (-0.5, 0.0, 2.0):
    Q = disk_rule(al, 70, 150); t=time.time(); B = DiskBasis(Q, 60)
    d = [delta_all_K(B,N,Q["w"],[60])[0] for N in (2,3,6)]
    print("alpha",al,"delta_{N,60}(nu):", ["%.1e"%v for v in d], "%.1fs"%(time.time()-t))
# (4.22) example
Q = disk_rule(0.0, 40, 96); B = DiskBasis(Q, 10)
t = Q["x"]**2+Q["y"]**2; h = 20*t**3-30*t**2+12*t-1
gu = (t**2+4*Q["x"]**2*t, 4*Q["x"]*Q["y"]*t)
C,_,_ = defect(B,3,gu,Q["w"]*(1+0.5*h)); print("(4.22):", vnorm(C,Q["w"])[0], np.sqrt(6)/70*0.5)
# radial restriction check (rho_* radial example, delta_3,K)
Q = disk_rule(0.0, 60, 64); w = Q["w"]*(2/3)*(1+Q["x"]**2+Q["y"]**2)*1.0
w = w/ w.sum()
B = DiskBasis(Q, 30); B3 = DiskBasis(Q, 30, mmax=3); B6 = DiskBasis(Q, 30, mmax=6)
print("rho* full  ", delta_all_K(B,3,w,[5,10,20,30]))
print("rho* m<=3  ", delta_all_K(B,3,w,[5,10,20,30],Binput=B3))
print("rho* m<=6  ", delta_all_K(B,3,w,[5,10,20,30],Binput=B6))
