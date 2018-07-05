#include <particle_simulator.hpp>
#include <omp.h>
#include "cuda_pointer.h"
#include "helper_cuda.h"
#include "DF.hpp"
#include "use_gpu.hpp"

__constant__ PS::F64 KernelSupportRadiusMaximum_C;
__constant__ DF eta_c;
__constant__ DF ksrh_c;
__constant__ DF ceff0_c;
__constant__ DF ceff1_c;
/*
  enum{
  N_THREAD_GPU = 32,
  N_WALK_LIMIT = 1000,
  NI_LIMIT = N_WALK_LIMIT*1000,
  NJ_LIMIT = N_WALK_LIMIT*10000,
  NumberOfDimension = 3
  };
*/
class DensityEPI
{
public:
  PS::S32    id;
  PS::F64    mass;
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64    ksr;
  PS::F64    rs;
};

class DensityEPJ
{
public:
  PS::S32    id;
  PS::F64    mass;
  PS::F64vec pos;
  PS::F64vec vel;
};

class Density
{
public:
  PS::F64 dens;
  PS::F64 divv;
  PS::F64 rotv;
  PS::F64 grdh;
  PS::F64 ksr;
  PS::S64 np;
  bool    itr;
};

struct Epi_DensGPU{
  long long int id;
  DF mass;
  DF3 pos;
  DF3 vel;
  DF ksr;
  DF rs;
  int id_walk;
  int nocalc;
};

struct Epj_DensGPU{
  long long int id;
  DF mass;
  DF3 pos;
  DF3 vel;
};

struct DensGPU{
  DF dens;
  DF divv;
  DF rotv;
  DF grdh;
  DF ksr;
  long long int np;
  bool itr;
};

inline __device__ DF kernel0th(const DF q)
{
  //const DF ceff0   = +3.342253804929802286e+00;
  const DF qmin  = ((1. - q > 0.) ? 1. - q : 0.);
  const DF qmin2 = qmin * qmin;
  return ceff0_c * qmin2 * qmin2 * (1. + 4. * q);
}

inline __device__ DF kernel1st(const DF q)
{
  //const DF ceff1 = +1.336901521971920914e+01;
  const DF qmin  = ((1. - q > 0.) ? 1. - q : 0.);
  const DF qmin2 = qmin  * qmin;
  const DF qmin3 = qmin  * qmin2;
  const DF qmin4 = qmin2 * qmin2;
  return ceff1_c * (qmin4 - qmin3 * (1. + 4. * q));
}

inline __device__ void Dens_calc1(const struct Epi_DensGPU ip,
				  struct Epj_DensGPU jp,
				  const DF hi_i,
				  const DF hi3_i,
				  struct DensGPU *dens)
{
  const DF dx_ij = ip.pos.x - jp.pos.x;
  const DF dy_ij = ip.pos.y - jp.pos.y;
  const DF dz_ij = ip.pos.z - jp.pos.z;

  const DF r2_ij = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

  const DF r1_ij = sqrt(r2_ij);
  const DF q_i   = r1_ij * hi_i;

  const DF kw0 = kernel0th(q_i);
  const DF rhj = jp.mass * hi3_i * kw0;
  
  dens->dens += rhj;
  dens->np += ((q_i < 1.) ? 1 : 0);        
}

inline __device__ void Dens_calc2(struct Epi_DensGPU ip,
				  const struct Epj_DensGPU jp,
				  const DF hi_i,
				  const DF hi4_i,
				  DF *grdh_i,
				  DF *divv_i,
				  DF *rotx_i,
				  DF *roty_i,
				  DF *rotz_i)
{
  DF dpx_ij = ip.pos.x - jp.pos.x;
  DF dpy_ij = ip.pos.y - jp.pos.y;
  DF dpz_ij = ip.pos.z - jp.pos.z;
  DF dvx_ij = ip.vel.x - jp.vel.x;
  DF dvy_ij = ip.vel.y - jp.vel.y;
  DF dvz_ij = ip.vel.z - jp.vel.z;

  DF r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
  //DF ri_ij = 1. / sqrt(r2_ij);
  DF ri_ij = rsqrtf(r2_ij);
  ri_ij = ((ip.id != jp.id) ? ri_ij : 0.);
  DF r1_ij = r2_ij * ri_ij;
  DF q_i = r1_ij * hi_i;

  DF kw0 = kernel0th(q_i);
  DF kw1 = kernel1st(q_i);

  DF m_j(jp.mass);
  DF ghj = (DF) NumberOfDimension * kw0;
	   
  ghj += q_i * kw1;
  *grdh_i -= ghj * hi4_i * m_j;

  DF dw_ij  = m_j * hi4_i * kw1 * ri_ij;
  DF dwx_ij = dw_ij * dpx_ij;
  DF dwy_ij = dw_ij * dpy_ij;
  DF dwz_ij = dw_ij * dpz_ij;

  *divv_i -= dvx_ij * dwx_ij;
  *divv_i -= dvy_ij * dwy_ij;
  *divv_i -= dvz_ij * dwz_ij;

  *rotx_i += dvy_ij * dwz_ij;
  *roty_i += dvz_ij * dwx_ij;
  *rotz_i += dvx_ij * dwy_ij;
  *rotx_i -= dvz_ij * dwy_ij;
  *roty_i -= dvx_ij * dwz_ij;
  *rotz_i -= dvy_ij * dwx_ij;  

}

__device__ void DensKernel_1walk(struct Epj_DensGPU *jpsh,
				 struct Epi_DensGPU ip,
				 const int id_walk,
				 const int2 *ij_disp,
				 const Epj_DensGPU *epj,
				 struct DensGPU *dens)
{
  const int tid = threadIdx.x;
  const int j_head = ij_disp[id_walk  ].y;
  const int j_tail = ij_disp[id_walk+1].y;

  for(int repeat=0; repeat<3; repeat++){
    const DF hi_i  = 1. / ip.ksr;
    const DF hi3_i = hi_i * hi_i * hi_i;
    //const DF hi4_i = hi_i * hi3_i;
    dens->dens = 0;
    dens->np = 0;
    //DF rh_i  = 0.;
    //long long int nj_i  = 0.;
    for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
      jpsh[tid] = ((Epj_DensGPU *)(epj + j))[tid];

      if(j_tail-j < N_THREAD_GPU){
	for(int jj=0; jj<j_tail-j; jj++){
	  Dens_calc1(ip, jpsh[jj], hi_i, hi3_i, dens);
	}
      } else {
#pragma unroll
	for(int jj=0; jj<N_THREAD_GPU; jj++){
	  Dens_calc1(ip, jpsh[jj], hi_i, hi3_i, dens);
	}
      }
    }
    /*
      const DF eta     = 1.6;
      const DF ksrh    = 1.936492;
    */
    DF buf_hs = eta_c * ksrh_c * powf(ip.mass / dens->dens, 1. / 3.);
  
    buf_hs = ((buf_hs < KernelSupportRadiusMaximum_C)
	      ? buf_hs : KernelSupportRadiusMaximum_C);

    dens->ksr = buf_hs;
    dens->itr = (buf_hs > ip.rs) ? true : false;
 
    ip.ksr = buf_hs;
  }

  const DF hi_i   = 1. / ip.ksr;
  const DF hi4_i  = hi_i * hi_i * hi_i * hi_i;
  DF grdh_i = 0.;
  DF divv_i = 0.;
  DF rotx_i = 0.;
  DF roty_i = 0.;
  DF rotz_i = 0.;

  for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
    jpsh[tid] = ((Epj_DensGPU *)(epj + j))[tid];
    if(j_tail-j < N_THREAD_GPU){
      for(int jj=0; jj<j_tail-j; jj++){
	Dens_calc2(ip, jpsh[jj], hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
      }
    } else {
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	Dens_calc2(ip, jpsh[jj], hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
      }
    }
  }

  DF dens_i = dens->dens;
  DF deni_i = 1. / dens_i;
  DF omgi_i = 1. / (1. + ip.ksr * deni_i * grdh_i / NumberOfDimension);
  omgi_i = (1. + ip.ksr * deni_i * grdh_i / NumberOfDimension != 0.) ? omgi_i : 1.;
  DF rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
  DF rotv_i = rot2_i * ((rot2_i != 0.) ? 1. / sqrt(rot2_i) : 0.);
  rotv_i *= deni_i * omgi_i;
  divv_i *= deni_i * omgi_i;

  dens->divv = divv_i;
  dens->rotv = rotv_i;
  dens->grdh = omgi_i;
}

__device__ void DensKernel_2walk(struct Epj_DensGPU (*jpsh)[N_THREAD_GPU],
				 struct Epi_DensGPU ip,
				 const int id_walk,
				 const int2 *ij_disp,
				 const Epj_DensGPU *epj,
				 struct DensGPU *dens,
				 const int iwalk0,
				 const int iwalk1)
{
  const int jbeg0 = ij_disp[iwalk0].y;
  const int jbeg1 = ij_disp[iwalk1].y;
  const int jend0 = ij_disp[iwalk0+1].y;
  const int jend1 = ij_disp[iwalk1+1].y;
  const int nj0 = jend0 - jbeg0;
  const int nj1 = jend1 - jbeg1;

  const int nj_longer = nj0 > nj1 ? nj0 : nj1;
  const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
  const int walk_longer = nj0 > nj1 ? 0 : 1;
  const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

  const int mywalk = id_walk==iwalk0 ? 0 : 1;

  const int tid = threadIdx.x;

  for(int repeat=0; repeat<3; repeat++){
    const DF hi_i  = 1. / ip.ksr;
    const DF hi3_i = hi_i * hi_i * hi_i;
    //const DF hi4_i = hi_i * hi3_i;
    dens->dens = 0;
    dens->np = 0;
    //DF rh_i  = 0.;
    //long long int nj_i  = 0.;
    for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
      jpsh[0][tid] = ((Epj_DensGPU *)(epj + jbeg0 + j))[tid];
      jpsh[1][tid] = ((Epj_DensGPU *)(epj + jbeg1 + j))[tid];
      if(nj_shorter-j < N_THREAD_GPU){
	for(int jj=0; jj<nj_shorter-j; jj++){
	  Dens_calc1(ip, jpsh[mywalk][jj], hi_i, hi3_i, dens);
	}
      }else {
#pragma unroll
	for(int jj=0; jj<N_THREAD_GPU; jj++){
	  Dens_calc1(ip, jpsh[mywalk][jj], hi_i, hi3_i, dens);
	}
      }
    }

    for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
      jpsh[0][tid] = ((Epj_DensGPU *)(epj + jbeg_longer + j))[tid];
      int jrem = nj_longer - j;
      if(jrem < N_THREAD_GPU){
	for(int jj=0; jj<jrem; jj++){
	  if(mywalk == walk_longer)
	    Dens_calc1(ip, jpsh[0][jj], hi_i, hi3_i, dens);
	}
      }else {
#pragma unroll
	for(int jj=0; jj<N_THREAD_GPU; jj++){
	  if(mywalk == walk_longer)
	    Dens_calc1(ip, jpsh[0][jj], hi_i, hi3_i, dens); 
	}
      }
    }
    /*
      const DF eta     = 1.6;
      const DF ksrh    = 1.936492;
    */
    DF buf_hs = eta_c * ksrh_c * powf(ip.mass / dens->dens, 1. / 3.);

    buf_hs = ((buf_hs < KernelSupportRadiusMaximum_C)
	      ? buf_hs : KernelSupportRadiusMaximum_C);

    dens->ksr = buf_hs;
    dens->itr = (buf_hs > ip.rs) ? true : false;
 
    ip.ksr = buf_hs;
  }

  const DF hi_i   = 1. / ip.ksr;
  const DF hi4_i  = hi_i * hi_i * hi_i * hi_i;
  DF grdh_i = 0.;
  DF divv_i = 0.;
  DF rotx_i = 0.;
  DF roty_i = 0.;
  DF rotz_i = 0.;

  for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_DensGPU *)(epj + jbeg0 + j))[tid];
    jpsh[1][tid] = ((Epj_DensGPU *)(epj + jbeg1 + j))[tid];
    if(nj_shorter-j < N_THREAD_GPU){
      for(int jj=0; jj<nj_shorter-j; jj++){
	Dens_calc2(ip, jpsh[mywalk][jj], hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
      }
    }else {
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	Dens_calc2(ip, jpsh[mywalk][jj], hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
      }
    }
  }

  for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_DensGPU *)(epj + jbeg_longer + j))[tid];
    int jrem = nj_longer - j;
    if(jrem < N_THREAD_GPU){
      for(int jj=0; jj<jrem; jj++){
	if(mywalk == walk_longer)
	  Dens_calc2(ip, jpsh[0][jj], hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
      }
    }else {
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	if(mywalk == walk_longer)
	  Dens_calc2(ip, jpsh[0][jj], hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
      }
    }
  }

  DF dens_i = dens->dens;
  DF deni_i = 1. / dens_i;
  DF omgi_i = 1. / (1. + ip.ksr * deni_i * grdh_i / NumberOfDimension);
  omgi_i = (1. + ip.ksr * deni_i * grdh_i / NumberOfDimension != 0.) ? omgi_i : 1.;
  DF rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
  DF rotv_i = rot2_i * ((rot2_i != 0.) ? 1. / sqrt(rot2_i) : 0.);
  rotv_i *= deni_i * omgi_i;
  divv_i *= deni_i * omgi_i;

  dens->divv = divv_i;
  dens->rotv = rotv_i;
  dens->grdh = omgi_i;
}

__device__ void DensKernel_multiwalk(struct Epi_DensGPU ip,
				     const int id_walk,
				     const int2 *ij_disp,
				     const Epj_DensGPU *epj,
				     struct DensGPU *dens)
{
  const int j_head = ij_disp[id_walk  ].y;
  const int j_tail = ij_disp[id_walk+1].y;
 
  for(int repeat=0; repeat<3; repeat++){
    const DF hi_i  = 1. / ip.ksr;
    const DF hi3_i = hi_i * hi_i * hi_i;
    //const DF hi4_i = hi_i * hi3_i;
    dens->dens = 0;
    dens->np = 0;
    //DF rh_i  = 0.;
    //long long int nj_i  = 0.;
    
    for(int j=j_head; j<j_tail; j++){
      const struct Epj_DensGPU jp = epj[j];
      Dens_calc1(ip, jp, hi_i, hi3_i, dens);
    }
    /*   
	 const DF eta     = 1.6;
	 const DF ksrh    = 1.936492;
    */
    DF buf_hs = eta_c * ksrh_c * powf(ip.mass / dens->dens, 1. / 3.);
  
    buf_hs = ((buf_hs < KernelSupportRadiusMaximum_C)
	      ? buf_hs : KernelSupportRadiusMaximum_C);
   
    dens->ksr = buf_hs;
    dens->itr = (buf_hs > ip.rs) ? true : false;
 
    ip.ksr = buf_hs;
  }
  const DF hi_i   = 1. / ip.ksr;
  const DF hi4_i  = hi_i * hi_i * hi_i * hi_i;
  DF grdh_i = 0.;
  DF divv_i = 0.;
  DF rotx_i = 0.;
  DF roty_i = 0.;
  DF rotz_i = 0.;
  
  for(int j=j_head; j<j_tail; j++){
    const struct Epj_DensGPU jp = epj[j];
    Dens_calc2(ip, jp, hi_i, hi4_i, &grdh_i, &divv_i, &rotx_i, &roty_i, &rotz_i); 
  }
  
  DF dens_i = dens->dens;
  DF deni_i = 1. / dens_i;
  DF omgi_i = 1. / (1. + ip.ksr * deni_i * grdh_i / NumberOfDimension);
  omgi_i = (1. + ip.ksr * deni_i * grdh_i / NumberOfDimension != 0.) ? omgi_i : 1.;
  DF rot2_i = rotx_i * rotx_i + roty_i * roty_i + rotz_i * rotz_i;
  DF rotv_i = rot2_i * ((rot2_i != 0.) ? 1. / sqrt(rot2_i) : 0.);
  rotv_i *= deni_i * omgi_i;
  divv_i *= deni_i * omgi_i;

  dens->divv = divv_i;
  dens->rotv = rotv_i;
  dens->grdh = omgi_i;
  
}

__global__ void DensKernel(const int2 *ij_disp,
			   const Epi_DensGPU *epi,
			   const Epj_DensGPU *epj,
			   DensGPU *dev_dens)  
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  struct Epi_DensGPU ip = epi[tid];
  const int id_walk = epi[tid].id_walk;
  struct DensGPU dens;

  const int t_head = blockDim.x * blockIdx.x;
  const int t_tail = t_head + N_THREAD_GPU - 1;
  const int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);

  __shared__ struct Epj_DensGPU jpsh[2][N_THREAD_GPU];
    
  if(1 == nwalk_in_block){
    DensKernel_1walk(jpsh[0], ip, id_walk, ij_disp, epj, &dens);
  } else if(2 == nwalk_in_block){
    const int iwalk0 = epi[t_head].id_walk;
    const int iwalk1 = epi[t_tail].id_walk;
    DensKernel_2walk(jpsh, ip, id_walk, ij_disp, epj, &dens, iwalk0, iwalk1);
  } else {
    DensKernel_multiwalk(ip, id_walk, ij_disp, epj, &dens);
  }
  
  //DensKernel_multiwalk(ip, id_walk, ij_disp, epj, &dens);

  dev_dens[tid] = dens;

}

static cudaPointer<Epi_DensGPU>   dev_epi;
static cudaPointer<Epj_DensGPU>   dev_epj;
static cudaPointer<DensGPU> dev_dens;
static cudaPointer<int2>     ij_disp;
static bool init_call = true;

namespace RunParameter {
  extern PS::F64 KernelSupportRadiusMaximum;
};

PS::S32 DispatchKernel_Dens(const PS::S32 tag,
			    const PS::S32 n_walk,
			    const DensityEPI **epi,
			    const PS::S32 *n_epi,
			    const DensityEPJ **epj,
			    const PS::S32 *n_epj)
{
  assert(n_walk <= N_WALK_LIMIT);

  if(init_call){
    dev_epi  .allocate(NI_LIMIT);
    dev_epj  .allocate(NJ_LIMIT);
    dev_dens .allocate(NI_LIMIT);
    ij_disp  .allocate(N_WALK_LIMIT+2);
    init_call = false;
  }

  ij_disp[0].x = 0;
  ij_disp[0].y = 0;
  for(int k=0; k<n_walk; k++){
    ij_disp[k+1].x = ij_disp[k].x + n_epi[k];
    ij_disp[k+1].y = ij_disp[k].y + n_epj[k];
  }
  ij_disp[n_walk+1] = ij_disp[n_walk];

  assert(ij_disp[n_walk].x < NI_LIMIT);
  assert(ij_disp[n_walk].y < NJ_LIMIT);
  ij_disp.htod(n_walk+2);

  int ni_tot_reg = ij_disp[n_walk].x;
  if(ni_tot_reg % N_THREAD_GPU){
    ni_tot_reg /= N_THREAD_GPU;
    ni_tot_reg++;
    ni_tot_reg *= N_THREAD_GPU;
  }

  int ni_tot = 0;
  int nj_tot = 0;
  
  int iw, i, j, num_i, num_j;

#pragma omp parallel for private(i,j, num_i, num_j) reduction(+:ni_tot, nj_tot)
  for(iw=0; iw<n_walk; iw++){
    for(i=0; i<n_epi[iw]; i++){
      num_i = ij_disp[iw].x + i;
      dev_epi[num_i].id = epi[iw][i].id;
      dev_epi[num_i].pos.x = epi[iw][i].pos.x;
      dev_epi[num_i].pos.y = epi[iw][i].pos.y;
      dev_epi[num_i].pos.z = epi[iw][i].pos.z;
      dev_epi[num_i].vel.x = epi[iw][i].vel.x;
      dev_epi[num_i].vel.y = epi[iw][i].vel.y;
      dev_epi[num_i].vel.z = epi[iw][i].vel.z;
      dev_epi[num_i].ksr = epi[iw][i].ksr;
      dev_epi[num_i].mass = epi[iw][i].mass;
      dev_epi[num_i].rs = epi[iw][i].rs;
      dev_epi[num_i].id_walk = iw;
      ni_tot++;
    }
    for(j=0; j<n_epj[iw]; j++){
      num_j = ij_disp[iw].y + j;
      dev_epj[num_j].id = epj[iw][j].id;
      dev_epj[num_j].pos.x = epj[iw][j].pos.x;
      dev_epj[num_j].pos.y = epj[iw][j].pos.y;
      dev_epj[num_j].pos.z = epj[iw][j].pos.z;
      dev_epj[num_j].vel.x = epj[iw][j].vel.x;
      dev_epj[num_j].vel.y = epj[iw][j].vel.y;
      dev_epj[num_j].vel.z = epj[iw][j].vel.z;
      dev_epj[num_j].mass = epj[iw][j].mass;
      nj_tot++;
    }
  }

  for(int i=ni_tot; i<ni_tot_reg; i++){
    dev_epi[i].id_walk = n_walk;
    dev_epi[i].nocalc = 1;
  }

  dev_epi.htod(ni_tot_reg);
  dev_epj.htod(nj_tot);

  const DF eta = 1.6;
  const DF ksrh = 1.936492;
  const DF ceff0   = +3.342253804929802286e+00;
  const DF ceff1 = +1.336901521971920914e+01;

  cudaMemcpyToSymbol(KernelSupportRadiusMaximum_C, &RunParameter::KernelSupportRadiusMaximum, sizeof(PS::F64));
  cudaMemcpyToSymbol(eta_c, &eta, sizeof(DF));
  cudaMemcpyToSymbol(ksrh_c, &ksrh, sizeof(DF));
  cudaMemcpyToSymbol(ceff0_c, &ceff0, sizeof(DF));
  cudaMemcpyToSymbol(ceff1_c, &ceff1, sizeof(DF));

  int nblocks  = ni_tot_reg / N_THREAD_GPU;
  int nthreads = N_THREAD_GPU;
  DensKernel <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_dens);

  return 0;
}

PS::S32 RetrieveKernel_Dens(const PS::S32 tag,
			    const PS::S32 n_walk,
			    const PS::S32 *ni,
			    Density **density)
{
  int ni_tot = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }

  dev_dens.dtoh(ni_tot);

  int n_cnt = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<ni[iw]; i++){
      density[iw][i].dens = dev_dens[n_cnt].dens;
      density[iw][i].divv = dev_dens[n_cnt].divv;
      density[iw][i].rotv = dev_dens[n_cnt].rotv;
      density[iw][i].grdh = dev_dens[n_cnt].grdh;
      density[iw][i].ksr = dev_dens[n_cnt].ksr;
      density[iw][i].np = dev_dens[n_cnt].np;
      density[iw][i].itr = dev_dens[n_cnt].itr;
      n_cnt++;
    }
  }

  return 0;
}
