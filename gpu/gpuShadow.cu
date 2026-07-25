// gpuShadow.cu — GPU synthetic shadowgraph. Modes: sharp | eikonal | hybrid.
// One CUDA thread per pixel; rays march through 3D texture(s), hardware trilinear sampling.
//  sharp   : threshold field0(alpha) at 0.5 -> Snell/TIR + Beer-Lambert absorption.
//  eikonal : continuous RK4 (Sharma) through n=1+K*field0.
//  hybrid  : gas eikonal bend (field1) + Snell interface (field0=alpha) + absorption (evap).
// Finite DoF: aperture sampling (naper rays/pixel converging at focal plane zf, radius apR).
// outmode: 0 shadowgraph | 1 mean eps_x [mrad] | 2 mean eps_y [mrad]
//          3 CLASSICAL SCHLIEREN (knife edge at the focal plane, finite source image)
//          4 BOS |displacement| [mm] | 5 BOS dx [mm] | 6 BOS dy [mm]
//   Classical schlieren: a source image of height a is cut by a knife so a fraction `cutoff`
//   passes undeflected; a ray deflected by eps is displaced f2*eps at the knife, so
//   T = clamp(cutoff + (f2/a)*eps, 0, 1)  -- the standard Settles contrast dI/I0 = f2*eps/a.
//   BOS: a background pattern a distance Lbg behind the test section appears displaced by
//   Lbg*tan(eps); cross-correlating that displacement field is what a BOS rig measures.
// grid.bin: int32 NX,NY,NZ ; float X0,X1,Y0,Y1,Z0,Z1 ; NX*NY*NZ float32 (k,j,i).
// Args: grid0 img mode RESX ACCdeg DN K absorb nLiq [outmode naper apR zf grid1 K1 srcAng knife cutoff sgain Lbg]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cuda_runtime.h>

struct Params {
  int NX,NY,NZ, RESX,RESY, mode, nstep, outmode, naper, knife;
  float X0,X1,Y0,Y1,Z0,Z1, ACCcos, DN, K, absorb, nLiq, ds, apR, zf, K1, srcAng, cutoff, sgain, Lbg;
};
__device__ __forceinline__ float3 mk(float a,float b,float c){ float3 r; r.x=a;r.y=b;r.z=c; return r; }
__device__ __forceinline__ float3 operator+(float3 a,float3 b){ return mk(a.x+b.x,a.y+b.y,a.z+b.z); }
__device__ __forceinline__ float3 operator-(float3 a,float3 b){ return mk(a.x-b.x,a.y-b.y,a.z-b.z); }
__device__ __forceinline__ float3 operator*(float3 a,float s){ return mk(a.x*s,a.y*s,a.z*s); }
__device__ __forceinline__ float dot3(float3 a,float3 b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
__device__ __forceinline__ float len3(float3 a){ return sqrtf(dot3(a,a)); }
__device__ __forceinline__ float3 norm3(float3 a){ float l=len3(a); return l>1e-20f? a*(1.f/l):a; }
__device__ __forceinline__ float sampleF(cudaTextureObject_t tex,const Params&P,float3 w){
  return tex3D<float>(tex,(w.x-P.X0)/(P.X1-P.X0)*P.NX,(w.y-P.Y0)/(P.Y1-P.Y0)*P.NY,(w.z-P.Z0)/(P.Z1-P.Z0)*P.NZ);
}
__device__ __forceinline__ float3 gradF(cudaTextureObject_t tex,const Params&P,float3 w){
  float hx=(P.X1-P.X0)/P.NX*0.5f, hy=(P.Y1-P.Y0)/P.NY*0.5f, hz=(P.Z1-P.Z0)/P.NZ*0.5f;
  return mk((sampleF(tex,P,mk(w.x+hx,w.y,w.z))-sampleF(tex,P,mk(w.x-hx,w.y,w.z)))/(2*hx),
            (sampleF(tex,P,mk(w.x,w.y+hy,w.z))-sampleF(tex,P,mk(w.x,w.y-hy,w.z)))/(2*hy),
            (sampleF(tex,P,mk(w.x,w.y,w.z+hz))-sampleF(tex,P,mk(w.x,w.y,w.z-hz)))/(2*hz));
}
__device__ __forceinline__ float3 refract(float3 d,float3 n,float eta){
  float ci=-dot3(n,d); if(ci<0){ n=n*-1.f; ci=-ci; }
  float k=1.f-eta*eta*(1.f-ci*ci);
  if(k<0.f) return norm3(d + n*(2.f*ci));
  return norm3(d*eta + n*(eta*ci - sqrtf(k)));
}
struct RayOut{ float3 dir; float absorb; };
// BOS validity: the ray must still be travelling forward (it has not been turned back or
// sideways by TIR at the interface) AND the background must still be visible through it.
// Rays failing either carry no background signal -- in a real rig the cross-correlation
// simply fails there -- so they are excluded rather than allowed to blow up tan(eps).
__device__ __forceinline__ bool bosValid(float3 d,float att,float minCos,float minT){
  return (d.z>minCos)&&(att>minT)&&isfinite(d.x)&&isfinite(d.y)&&isfinite(d.z);
}
__device__ RayOut marchRay(cudaTextureObject_t t0,cudaTextureObject_t t1,const Params&P,float3 pos,float3 dir){
  float absorb=0.f;
  if(P.mode==1){ // eikonal Sharma RK4
    float n0=1.f+P.K*sampleF(t0,P,pos); float3 D=dir*n0;
    for(int s=0;s<P.nstep;s++){ if(pos.z>=P.Z1)break; float dt=P.ds;
      #define AF(pp) ( gradF(t0,P,(pp))*(P.K*(1.f+P.K*sampleF(t0,P,(pp)))) )
      float3 a1=D,a1D=AF(pos); float3 a2=D+a1D*(0.5f*dt),a2D=AF(pos+a1*(0.5f*dt));
      float3 a3=D+a2D*(0.5f*dt),a3D=AF(pos+a2*(0.5f*dt)); float3 a4=D+a3D*dt,a4D=AF(pos+a3*dt);
      pos=pos+(a1+a2*2.f+a3*2.f+a4)*(dt/6.f); D=D+(a1D+a2D*2.f+a3D*2.f+a4D)*(dt/6.f);
      #undef AF
    }
    dir=norm3(D);
  } else { // sharp (mode 0) or hybrid (mode 2)
    bool inside=sampleF(t0,P,pos)>=0.5f;
    for(int s=0;s<P.nstep;s++){
      if(P.mode==2 && !inside){ float ng=1.f+P.K1*sampleF(t1,P,pos);
        float3 gn=gradF(t1,P,pos)*P.K1; float3 gp=gn-dir*dot3(gn,dir); dir=norm3(dir+gp*(P.ds/ng)); }
      pos=pos+dir*P.ds; if(pos.z>=P.Z1)break;
      bool now=sampleF(t0,P,pos)>=0.5f;
      if(now!=inside){ float3 nrm=norm3(gradF(t0,P,pos)); float eta=now?(1.f/P.nLiq):P.nLiq; dir=refract(dir,nrm,eta); inside=now; }
      if(inside) absorb+=P.absorb*P.ds;
    }
  }
  RayOut o; o.dir=dir; o.absorb=absorb; return o;
}
__global__ void traceKernel(cudaTextureObject_t t0,cudaTextureObject_t t1,Params P,float* img){
  int ix=blockIdx.x*blockDim.x+threadIdx.x, iy=blockIdx.y*blockDim.y+threadIdx.y;
  if(ix>=P.RESX||iy>=P.RESY) return;
  float x=P.X0+(ix+0.5f)/P.RESX*(P.X1-P.X0), y=P.Y0+(iy+0.5f)/P.RESY*(P.Y1-P.Y0);
  float sumT=0.f, dfx=0.f, dfy=0.f, schl=0.f, bosx=0.f, bosy=0.f; int nbos=0;
  for(int ai=0;ai<P.naper;ai++){
    float3 start,dir,entry;
    if(P.srcAng>0.f){                        // diffuse extended source: sample incident ray over the source cone
      float sa=ai*2.399963f, sr=P.srcAng*sqrtf((ai+0.5f)/P.naper);
      start=mk(x,y,P.Z0); dir=norm3(mk(tanf(sr)*cosf(sa),tanf(sr)*sinf(sa),1.f)); entry=mk(0,0,1); }
    else if(P.naper<=1||P.apR<=0.f){ start=mk(x,y,P.Z0); dir=mk(0,0,1); entry=dir; }
    else{ float ang=ai*2.399963f, rr=P.apR*sqrtf((ai+0.5f)/P.naper);   // thin-lens aperture disk
          start=mk(x+rr*cosf(ang),y+rr*sinf(ang),P.Z0); dir=norm3(mk(x,y,P.zf)-start); entry=norm3(dir); }
    RayOut o=marchRay(t0,t1,P,start,dir);
    float cosang=dot3(norm3(o.dir),entry);   // net deviation from the undeviated aperture path
    float att=expf(-o.absorb);
    sumT+=(cosang>=P.ACCcos)? att:0.f;
    dfx+=atan2f(o.dir.x,o.dir.z); dfy+=atan2f(o.dir.y,o.dir.z);
    // deflection RELATIVE to this ray's own incident direction (correct for aperture /
    // diffuse sources, where the undeviated ray is not along +z)
    float ex=atan2f(o.dir.x,o.dir.z)-atan2f(entry.x,entry.z);
    float ey=atan2f(o.dir.y,o.dir.z)-atan2f(entry.y,entry.z);
    float eps=(P.knife==1)? ex : ey;                       // knife==1 cuts x-deflection
    float T=P.cutoff+P.sgain*eps; T=fminf(fmaxf(T,0.f),1.f);
    schl+=T*att;
    // exclude rays that cannot reach the background (TIR / opaque liquid): cos(85deg)=0.0872
    if(bosValid(o.dir,att,0.0872f,0.02f)){ bosx+=P.Lbg*tanf(ex); bosy+=P.Lbg*tanf(ey); nbos++; }
  }
  float inv=1.f/P.naper, out;
  switch(P.outmode){
    case 1: out=dfx*inv*1000.f; break;                       // mean eps_x [mrad]
    case 2: out=dfy*inv*1000.f; break;                       // mean eps_y [mrad]
    case 3: out=schl*inv; break;                             // classical schlieren intensity
    // BOS modes report NaN where no ray reached the background (masked in a real rig)
    case 4: out=nbos? sqrtf(bosx*bosx+bosy*bosy)/nbos*1000.f : nanf(""); break;
    case 5: out=nbos? bosx/nbos*1000.f : nanf(""); break;
    case 6: out=nbos? bosy/nbos*1000.f : nanf(""); break;
    default: out=sumT*inv;                                   // shadowgraph
  }
  img[iy*P.RESX+ix]=out;
}
static cudaTextureObject_t upload(const char* path,int&NX,int&NY,int&NZ,float bnd[6]){
  FILE* f=fopen(path,"rb"); if(!f){printf("no grid %s\n",path);exit(1);}
  int dims[3]; if(fread(dims,4,3,f)!=3||fread(bnd,4,6,f)!=6){printf("hdr\n");exit(1);}
  NX=dims[0];NY=dims[1];NZ=dims[2]; size_t N=(size_t)NX*NY*NZ;
  float* h=(float*)malloc(N*4); if(fread(h,4,N,f)!=N){printf("data\n");exit(1);} fclose(f);
  cudaArray_t arr; cudaChannelFormatDesc ch=cudaCreateChannelDesc<float>(); cudaExtent ext=make_cudaExtent(NX,NY,NZ);
  cudaMalloc3DArray(&arr,&ch,ext);
  cudaMemcpy3DParms cp={}; cp.srcPtr=make_cudaPitchedPtr(h,NX*4,NX,NY); cp.dstArray=arr; cp.extent=ext; cp.kind=cudaMemcpyHostToDevice; cudaMemcpy3D(&cp);
  cudaResourceDesc rd={}; rd.resType=cudaResourceTypeArray; rd.res.array.array=arr;
  cudaTextureDesc td={}; td.filterMode=cudaFilterModeLinear; td.addressMode[0]=td.addressMode[1]=td.addressMode[2]=cudaAddressModeClamp; td.normalizedCoords=0;
  cudaTextureObject_t tex; cudaCreateTextureObject(&tex,&rd,&td,0); return tex;
}
int main(int argc,char**argv){
  if(argc<10){ printf("usage: grid0 img mode RESX ACCdeg DN K absorb nLiq [outmode naper apR zf grid1 K1 srcAng knife cutoff sgain Lbg]\n  outmode: 0 shadowgraph, 1 eps_x mrad, 2 eps_y mrad, 3 schlieren, 4 BOS|d| mm, 5 BOS dx, 6 BOS dy\n"); return 1; }
  Params P; float bnd[6]; int NX,NY,NZ;
  cudaTextureObject_t t0=upload(argv[1],NX,NY,NZ,bnd);
  P.mode = strcmp(argv[3],"sharp")==0?0 : strcmp(argv[3],"eikonal")==0?1 : 2;
  int RESX=atoi(argv[4]); P.ACCcos=cosf(atof(argv[5])*3.14159265f/180.f);
  P.DN=atof(argv[6]); P.K=atof(argv[7]); P.absorb=atof(argv[8]); P.nLiq=atof(argv[9]);
  P.outmode=argc>10?atoi(argv[10]):0; P.naper=argc>11?atoi(argv[11]):1;
  P.apR=argc>12?atof(argv[12]):0.f; P.zf=argc>13?atof(argv[13]):0.5f*(bnd[4]+bnd[5]);
  P.K1=argc>15?atof(argv[15]):0.f;
  P.srcAng=argc>16?atof(argv[16])*3.14159265f/180.f:0.f;   // diffuse-source half-angle (deg)
  P.knife =argc>17?atoi(argv[17]):2;      // 1 = knife cuts eps_x (vertical edge), 2 = cuts eps_y
  P.cutoff=argc>18?atof(argv[18]):0.5f;   // fraction of the source image passed at zero deflection
  P.sgain =argc>19?atof(argv[19]):1000.f; // f2/a  [1/rad]  schlieren sensitivity
  P.Lbg   =argc>20?atof(argv[20]):0.f;    // BOS background distance behind the test section
  cudaTextureObject_t t1=t0;
  if(P.mode==2 && argc>14){ int a,b,c; float bb[6]; t1=upload(argv[14],a,b,c,bb); }
  P.NX=NX;P.NY=NY;P.NZ=NZ; P.X0=bnd[0];P.X1=bnd[1];P.Y0=bnd[2];P.Y1=bnd[3];P.Z0=bnd[4];P.Z1=bnd[5];
  P.RESX=RESX; P.RESY=(int)(RESX*(P.Y1-P.Y0)/(P.X1-P.X0)+0.5f);
  P.ds=(P.Z1-P.Z0)/(NZ*2.f); P.nstep=NZ*3; if(P.naper<1)P.naper=1;
  float* dimg; cudaMalloc(&dimg,(size_t)P.RESX*P.RESY*4);
  dim3 blk(16,16), grd((P.RESX+15)/16,(P.RESY+15)/16);
  cudaEvent_t e0,e1; cudaEventCreate(&e0);cudaEventCreate(&e1); cudaEventRecord(e0);
  traceKernel<<<grd,blk>>>(t0,t1,P,dimg); cudaEventRecord(e1); cudaEventSynchronize(e1);
  float ms=0; cudaEventElapsedTime(&ms,e0,e1);
  cudaError_t e=cudaGetLastError(); if(e){printf("CUDA err: %s\n",cudaGetErrorString(e));return 1;}
  float* himg=(float*)malloc((size_t)P.RESX*P.RESY*4); cudaMemcpy(himg,dimg,(size_t)P.RESX*P.RESY*4,cudaMemcpyDeviceToHost);
  FILE* o=fopen(argv[2],"wb"); int od[2]={P.RESX,P.RESY}; fwrite(od,4,2,o); fwrite(himg,4,(size_t)P.RESX*P.RESY,o); fclose(o);
  printf("GPU %s naper=%d: %dx%d in %.2f ms (%.1f Mrays/s)\n",argv[3],P.naper,P.RESX,P.RESY,ms,(P.RESX*(double)P.RESY*P.naper)/ms/1e3);
  return 0;
}
