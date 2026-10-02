import numpy as np, time
from diskpoly import *
sect=[0,np.pi/2,np.pi,1.5*np.pi,2*np.pi]
K=32
Q=disk_rule(0.0,K+16,K+16,sectors=sect,radial_in_r=True)
x,y=Q["x"],Q["y"]; g=np.abs(x)+np.abs(y); w=Q["w"]*(1+0.5*(g-8/(3*np.pi)))
B=DiskBasis(Q,4); Bf=DiskBasis(Q,K)
classes={"A1":lambda m,k: m%4==0 and k=="c","A2":lambda m,k: m%4==0 and k=="s",
         "B1":lambda m,k: m%4==2 and k=="c","B2":lambda m,k: m%4==2 and k=="s",
         "E":lambda m,k: m%2==1 and k=="c"}
Ks=[10,20,32]
for N in (2,3):
    full=delta_all_K(B,N,w,Ks,Binput=Bf)
    per={c:delta_all_K(B,N,w,Ks,Binput=DiskBasis(Q,K,select=f)) for c,f in classes.items()}
    mx=np.max(np.array(list(per.values())),axis=0)
    print(N,"full",full,"max-class",mx,"argmax",[max(per,key=lambda c:per[c][i]) for i in range(len(Ks))])
