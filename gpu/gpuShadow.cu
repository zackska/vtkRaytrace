// gpuShadow.cu — GPU synthetic shadowgraph. Modes: sharp | eikonal | hybrid.
// One CUDA thread per pixel; rays march through 3D texture(s), hardware trilinear sampling.
//  sharp   : threshold field0(alpha) at 0.5 -> Snell/TIR + Beer-Lambert absorption.
//  eikonal : continuous RK4 (Sharma) through n=1+K*field0.
//  hybrid  : gas eikonal bend (field1) + Snell interface (field0=alpha) + absorption (evap).
// Finite DoF: aperture sampling (naper rays/pixel converging at focal plane zf, radius apR).
// grid.bin: int32 NX,NY,NZ ; float X0,X1,Y0,Y1,Z0,Z1 ; NX*NY*NZ float32 (k,j,i).
// Args: grid0 img mode RESX ACCdeg DN K absorb nLiq [outmode naper apR zf grid1 K1]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cuda_runtime.h>

struct Params {
  int NX,NY,NZ, RESX,RESY, mode, nstep, outmode, naper;
  float X0,X1,Y0,Y1,Z0,Z1, ACCcos, DN, K, absorb, nLiq, ds, apR, zf, K1;
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
  float sumT=0.f, dfx=0.f, dfy=0.f;
  for(int ai=0;ai<P.naper;ai++){
    float3 start,dir;
    if(P.naper<=1||P.apR<=0.f){ start=mk(x,y,P.Z0); dir=mk(0,0,1); }
    else{ float ang=ai*2.399963f, rr=P.apR*sqrtf((ai+0.5f)/P.naper);   // thin-lens aperture disk
          start=mk(x+rr*cosf(ang),y+rr*sinf(ang),P.Z0); dir=norm3(mk(x,y,P.zf)-start); }
    float3 entry=norm3(dir);                 // collection is relative to the (possibly tilted) chief ray
    RayOut o=marchRay(t0,t1,P,start,dir);
    float cosang=dot3(norm3(o.dir),entry);   // net deviation from the undeviated aperture path
    sumT+=(cosang>=P.ACCcos)? expf(-o.absorb):0.f;
    dfx+=atan2f(o.dir.x,o.dir.z); dfy+=atan2f(o.dir.y,o.dir.z);
  }
  float inv=1.f/P.naper, out;
  if(P.outmode==1) out=dfx*inv*1000.f; else if(P.outmode==2) out=dfy*inv*1000.f; else out=sumT*inv;
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
  if(argc<10){ printf("usage: grid0 img mode RESX ACCdeg DN K absorb nLiq [outmode naper apR zf grid1 K1]\n"); return 1; }
  Params P; float bnd[6]; int NX,NY,NZ;
  cudaTextureObject_t t0=upload(argv[1],NX,NY,NZ,bnd);
  P.mode = strcmp(argv[3],"sharp")==0?0 : strcmp(argv[3],"eikonal")==0?1 : 2;
  int RESX=atoi(argv[4]); P.ACCcos=cosf(atof(argv[5])*3.14159265f/180.f);
  P.DN=atof(argv[6]); P.K=atof(argv[7]); P.absorb=atof(argv[8]); P.nLiq=atof(argv[9]);
  P.outmode=argc>10?atoi(argv[10]):0; P.naper=argc>11?atoi(argv[11]):1;
  P.apR=argc>12?atof(argv[12]):0.f; P.zf=argc>13?atof(argv[13]):0.5f*(bnd[4]+bnd[5]);
  P.K1=argc>15?atof(argv[15]):0.f;
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
