"""Analytic validation of the new gpuShadow outmodes (classical schlieren, BOS).

Test field: a LINEAR refractive-index ramp in x, n = 1 + K*field with field = G*x.
For a ray crossing a slab of length L with a transverse index gradient, the deflection is

    eps_x = (1/n0) * (dn/dx) * L  ~=  K*G*L        (n0 ~= 1, small angle)

which is exact for a uniform gradient. That gives three independent checks:

  outmode 1  mean eps_x [mrad]      -> K*G*L*1000
  outmode 5  BOS dx     [mm]        -> Lbg*tan(eps_x)*1000
  outmode 3  schlieren  [0-1]       -> clamp(cutoff + (f2/a)*eps_x, 0, 1)

The schlieren check is run at several `sgain` values to confirm the response is LINEAR in
sensitivity and saturates at 0 and 1, which is the physical knife-edge behaviour.
"""
import numpy as np, struct, subprocess, os, sys

W='/scratch/run/gpuShadow/'
NX,NY,NZ=64,64,256
X0,X1,Y0,Y1,Z0,Z1=0.0,0.02,0.0,0.02,0.0,0.08     # metres; L = 0.08 m
G=50.0                                            # field gradient [1/m]: field = G*x
K=1.0e-4                                          # n = 1 + K*field  -> dn/dx = K*G

def write_grid(path):
    x=np.linspace(X0,X1,NX,dtype=np.float32)
    fld=np.broadcast_to(G*x[None,None,:],(NZ,NY,NX)).astype(np.float32)  # (k,j,i), x fastest
    with open(path,'wb') as f:
        f.write(struct.pack('3i',NX,NY,NZ))
        f.write(struct.pack('6f',X0,X1,Y0,Y1,Z0,Z1))
        fld.tofile(f)

def read_img(path):
    with open(path,'rb') as f:
        rx,ry=struct.unpack('2i',f.read(8))
        return np.fromfile(f,np.float32).reshape(ry,rx)

def run(outmode,extra=()):
    #      grid img  mode    RESX ACC DN   K   absorb nLiq outmode naper apR zf grid1 K1 srcAng knife cutoff sgain Lbg
    cmd=[W+'gpuShadow',W+'val.bin',W+'valo.bin','eikonal','128','90','0.33',str(K),'0','1.0',
         str(outmode),'1','0','0.04',W+'val.bin','0','0']+[str(e) for e in extra]
    r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode!=0: print('FAIL:',r.stdout,r.stderr); sys.exit(1)
    return read_img(W+'valo.bin')

write_grid(W+'val.bin')
L=Z1-Z0
eps_pred=K*G*L                    # radians
print('Test slab: L=%.3f m, dn/dx=%.3e /m  ->  predicted eps_x = %.4f mrad'%(L,K*G,eps_pred*1e3))
print()

ok=True
# ---- 1. raw deflection ----
img=run(1); meas=float(np.median(img))
err=100*(meas-eps_pred*1e3)/(eps_pred*1e3)
print('[1] deflection eps_x   pred %8.4f mrad   meas %8.4f mrad   err %+6.2f%%  %s'
      %(eps_pred*1e3,meas,err,'OK' if abs(err)<5 else 'FAIL'))
ok &= abs(err)<5

# ---- 2. BOS displacement ----
Lbg=0.30
img=run(5,extra=(2,0.5,1000,Lbg)); meas=float(np.median(img))
pred=Lbg*np.tan(eps_pred)*1e3
err=100*(meas-pred)/pred
print('[5] BOS dx (Lbg=%.2f m)  pred %8.4f mm     meas %8.4f mm     err %+6.2f%%  %s'
      %(Lbg,pred,meas,err,'OK' if abs(err)<5 else 'FAIL'))
ok &= abs(err)<5

# ---- 3. schlieren: linearity in sensitivity, and saturation ----
print()
print('[3] classical schlieren, knife=1 (cuts eps_x), cutoff=0.5')
print('    %10s %10s %10s %8s'%('f2/a[1/rad]','predicted','measured','status'))
for sg in [0,100,500,1000,2000,5000,20000]:
    img=run(3,extra=(1,0.5,sg,0)); meas=float(np.median(img))
    pred=min(max(0.5+sg*eps_pred,0.0),1.0)
    good=abs(meas-pred)<0.02
    ok &= good
    print('    %10d %10.4f %10.4f %8s'%(sg,pred,meas,'OK' if good else 'FAIL'))
print()
print('OVERALL:','PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
