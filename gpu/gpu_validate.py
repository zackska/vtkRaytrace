import numpy as np, subprocess, struct
NX,NY,NZ=64,64,160
X0,X1,Y0,Y1,Z0,Z1=0.,0.02,0.,0.02,0.,0.02
kappa=0.5
ys=(Y0+(np.arange(NY)+0.5)/NY*(Y1-Y0)).astype(np.float32)
field=np.broadcast_to(ys[None,:,None],(NZ,NY,NX)).astype(np.float32).copy()  # field = y
with open('grid.bin','wb') as f:
    f.write(struct.pack('3i',NX,NY,NZ)); f.write(struct.pack('6f',X0,X1,Y0,Y1,Z0,Z1)); field.tofile(f)
r=subprocess.run(['./gpuShadow','grid.bin','img.bin','eikonal','160','8','0',str(kappa),'0','1.33','2'],capture_output=True,text=True)
print(r.stdout.strip())
with open('img.bin','rb') as f:
    rx,ry=struct.unpack('2i',f.read(8)); img=np.fromfile(f,np.float32).reshape(ry,rx)
cy,cx=ry//2,rx//2; gpu=float(img[cy,cx])
L=Z1-Z0; y_c=Y0+(cy+0.5)/ry*(Y1-Y0); n0=1+kappa*y_c; ana=kappa/n0*L*1000
print('EIKONAL: GPU deflection_y=%.4f mrad  analytic paraxial=%.4f mrad  rel err %.2e'%(gpu,ana,abs(gpu-ana)/ana))
