#include <particle_simulator.hpp>
#include <cuda_runtime.h>
#include <omp.h>
#include "cuda_pointer.h"
#include <helper_cuda.h>
#include "DF.hpp"
#include "use_gpu.hpp"
/*
  enum{
  N_THREAD_GPU = 32,
  N_WALK_LIMIT = 1000,
  NI_LIMIT = N_WALK_LIMIT*1000,
  NJ_LIMIT = N_WALK_LIMIT*10000,
  };
*/
class GravityEPI{
public:
  PS::F64vec pos;
  PS::F64 eps2;
};

class GravityEPJ{
public:
  PS::F64 mass;
  PS::F64vec pos;
  PS::F64 eps2;
};

class PS::GravitySPJ{
public:
  PS::F64 mass;
  PS::F64vec pos;
  PS::F64 eps2;
};

class Gravity{
public:
  PS::F64vec acc;
  PS::F64 pot;
  PS::F64 eta;
};

struct Epi_GravityGPU{
  DF3 pos;
  DF eps2;
  int id_walk;
};

struct Epj_GravityGPU{
  DF mass;
  DF3 pos;
  DF eps2;
};

struct Spj_GravityGPU{
  DF mass;
  DF3 pos;
  DF eps2;
};

struct GravityGPU{
  DF3 acc;
  DF pot;
  DF eta;
};

inline __device__ void GravityForce_calc(const struct Epi_GravityGPU ip,
					 const struct Epj_GravityGPU jp,
					 struct GravityGPU *gravity)
{
  DF dpx_ij = ip.pos.x - jp.pos.x;
  DF dpy_ij = ip.pos.y - jp.pos.y;
  DF dpz_ij = ip.pos.z - jp.pos.z;

  DF r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;

  DF re2_i = r2_ij + ip.eps2;
  //DF rei_i = 1. / sqrt(re2_i);
  DF rei_i = rsqrtf(re2_i);
  DF rei3_i = rei_i * rei_i * rei_i;

  DF re2_j = r2_ij + jp.eps2;
  //DF rei_j = 1. / sqrt(re2_j);
  DF rei_j = rsqrtf(re2_j);
  DF rei3_j = rei_j * rei_j * rei_j;
  
  DF dg2_ij = jp.mass * (rei3_i + rei3_j);

  gravity->pot -= jp.mass * (rei_i  + rei_j);
  gravity->acc.x -= dpx_ij * dg2_ij;
  gravity->acc.y -= dpy_ij * dg2_ij;
  gravity->acc.z -= dpz_ij * dg2_ij;
  gravity->eta += jp.mass * rei3_i;

}

__device__ void GravityForceKernel_1walk(struct Epj_GravityGPU *jpsh,
					 const struct Epi_GravityGPU ip,
					 const int id_walk,
					 const int2 *ij_disp,
					 const Epj_GravityGPU *epj,
					 struct GravityGPU *gravity)
{
  const int tid = threadIdx.x;
  const int j_head = ij_disp[id_walk  ].y;
  const int j_tail = ij_disp[id_walk+1].y;

  for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
    jpsh[tid] = ((Epj_GravityGPU *)(epj + j)) [tid];
   
    if(j_tail-j < N_THREAD_GPU){
      for(int jj=0; jj<j_tail-j; jj++){
	GravityForce_calc(ip, jpsh[jj], gravity);
      }
    }else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	GravityForce_calc(ip, jpsh[jj], gravity);
      }
    }
  }
	
  gravity->acc.x *= 0.5;
  gravity->acc.y *= 0.5;
  gravity->acc.z *= 0.5;
  gravity->pot *= 0.5;

}

__device__ void GravityForceKernel_2walk(struct Epj_GravityGPU (*jpsh)[N_THREAD_GPU],
					 const struct Epi_GravityGPU ip,
					 const int id_walk,
					 const int2 *ij_disp,
					 const Epj_GravityGPU *epj,
					 struct GravityGPU *gravity,
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

  for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_GravityGPU *)(epj + jbeg0 + j))[tid];
    jpsh[1][tid] = ((Epj_GravityGPU *)(epj + jbeg1 + j))[tid];
    if(nj_shorter-j < N_THREAD_GPU){
      for(int jj=0; jj<nj_shorter-j; jj++){	
	GravityForce_calc(ip, jpsh[mywalk][jj], gravity);
      }
    } else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	GravityForce_calc(ip, jpsh[mywalk][jj], gravity);
      }
    }
  }

  for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_GravityGPU *)(epj + jbeg_longer + j))[tid];
    int jrem = nj_longer - j;
    if(jrem < N_THREAD_GPU){
      for(int jj=0; jj<jrem; jj++){
	if(mywalk == walk_longer)
	  GravityForce_calc(ip, jpsh[0][jj], gravity);
      }
    } else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	if(mywalk == walk_longer)
	  GravityForce_calc(ip, jpsh[0][jj], gravity);
      }
    }
  }

  gravity->acc.x *= 0.5;
  gravity->acc.y *= 0.5;
  gravity->acc.z *= 0.5;
  gravity->pot *= 0.5;
}

struct Walk {
  int nj;
  int jbeg;
  int nwalk;
  int id_walk;
};

inline __device__ void swap(struct Walk *walk_x, struct Walk *walk_y)
{
  struct Walk temp = *walk_x;
  *walk_x = *walk_y;
  *walk_y = temp;
}

__device__ void GravityForceKernel_3walk(struct Epj_GravityGPU (*jpsh)[N_THREAD_GPU],
					 const struct Epi_GravityGPU ip,
					 const int id_walk,
					 const int2 *ij_disp,
					 const Epj_GravityGPU *epj,
					 struct GravityGPU *gravity,
					 const int iwalk0,
					 const int iwalk1,
					 const int iwalk2)
{
  struct Walk walk[3];
  walk[0].jbeg  = ij_disp[iwalk0].y;
  walk[1].jbeg  = ij_disp[iwalk1].y;
  walk[2].jbeg  = ij_disp[iwalk2].y;
  walk[0].nj = ij_disp[iwalk0+1].y - walk[0].jbeg;
  walk[1].nj = ij_disp[iwalk1+1].y - walk[1].jbeg;
  walk[2].nj = ij_disp[iwalk2+1].y - walk[2].jbeg;
  walk[0].nwalk = 0;
  walk[1].nwalk = 1;
  walk[2].nwalk = 2;
  walk[0].id_walk = iwalk0;
  walk[1].id_walk = iwalk1;
  walk[2].id_walk = iwalk2;

  if(walk[0].nj > walk[1].nj) swap(&walk[0], &walk[1]);
  if(walk[1].nj > walk[2].nj) swap(&walk[1], &walk[2]);
  if(walk[0].nj > walk[1].nj) swap(&walk[0], &walk[1]);

  const int nj_longer = walk[2].nj;
  const int nj_middle = walk[1].nj;
  const int nj_shorter = walk[0].nj;

  const int walk_longer = walk[2].nwalk;
  // const int walk_middle = walk[1].nwalk;
  const int walk_shorter = walk[0].nwalk;

  const int jbeg_longer = walk[2].jbeg;
  const int jbeg_middle = walk[1].jbeg;

  const int mywalk = id_walk - iwalk0; 

  // printf("%d %d %d %d %d %d\n", walk_longer, walk_middle, walk_shorter, jend0 - jbeg0, jend1 - jbeg1, jend2 - jbeg2);
  const int tid = threadIdx.x;

  for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_GravityGPU *)(epj + walk[0].jbeg + j))[tid];
    jpsh[1][tid] = ((Epj_GravityGPU *)(epj + walk[1].jbeg + j))[tid];
    jpsh[2][tid] = ((Epj_GravityGPU *)(epj + walk[2].jbeg + j))[tid];
    if(nj_shorter-j < N_THREAD_GPU){
      for(int jj=0; jj<nj_shorter-j; jj++){	
	GravityForce_calc(ip, jpsh[mywalk][jj], gravity);
      }
    } else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	GravityForce_calc(ip, jpsh[mywalk][jj], gravity);
      }
    }
  }

  for(int j = nj_shorter; j < nj_middle; j += N_THREAD_GPU){
    jpsh[walk[1].nwalk][tid] = ((Epj_GravityGPU *)(epj + jbeg_middle + j))[tid];
    jpsh[walk[2].nwalk][tid] = ((Epj_GravityGPU *)(epj + jbeg_longer + j))[tid];
    int jrem = nj_middle - j;
    if(jrem < N_THREAD_GPU){
      for(int jj=0; jj<jrem; jj++){
	if(mywalk != walk_shorter) {
	  GravityForce_calc(ip, jpsh[mywalk][jj], gravity);
	}
      }
    } else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	if(mywalk != walk_shorter)
	  GravityForce_calc(ip, jpsh[mywalk][jj], gravity);
      }
    }
  }

  for(int j=nj_middle; j<nj_longer; j+=N_THREAD_GPU){
    jpsh[walk[2].nwalk][tid] = ((Epj_GravityGPU *)(epj + jbeg_longer + j))[tid];
    int jrem = nj_longer - j;
    if(jrem < N_THREAD_GPU){
      for(int jj=0; jj<jrem; jj++){
	if(mywalk == walk_longer)
	  GravityForce_calc(ip, jpsh[walk[2].nwalk][jj], gravity);
      }
    } else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	if(mywalk == walk_longer)
	  GravityForce_calc(ip, jpsh[walk[2].nwalk][jj], gravity);
      }
    }
  }

  gravity->acc.x *= 0.5;
  gravity->acc.y *= 0.5;
  gravity->acc.z *= 0.5;
  gravity->pot *= 0.5;
}

__device__ void GravityForceKernel_multiwalk(const struct Epi_GravityGPU ip,
					     const int id_walk,
					     const int2 *ij_disp,
					     const Epj_GravityGPU *epj,
					     struct GravityGPU *gravity)
{
  const int j_head = ij_disp[id_walk  ].y;
  const int j_tail = ij_disp[id_walk+1].y;

  for(int j=j_head; j<j_tail; j++){
    const struct Epj_GravityGPU jp = epj[j];
    GravityForce_calc(ip, jp, gravity);
  }

  gravity->acc.x *= 0.5;
  gravity->acc.y *= 0.5;
  gravity->acc.z *= 0.5;
  gravity->pot *= 0.5;
  
}

__global__ void GravityForceKernel(const int2 *ij_disp,
				   const Epi_GravityGPU *epi,
				   const Epj_GravityGPU *epj,
				   GravityGPU *dev_gravity)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  int t_head = blockDim.x * blockIdx.x;
  int t_tail = t_head + N_THREAD_GPU - 1;
  int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);
  const int id_walk = epi[tid].id_walk;

  const struct Epi_GravityGPU ip = epi[tid];
  struct GravityGPU gravity;
#ifdef USE_FLOAT
  gravity.acc = make_float3(0.f, 0.f, 0.f);
#endif
#ifdef USE_DOUBLE
  gravity.acc = make_double3(0.f, 0.f, 0.f);
#endif
  gravity.pot = 0.f;
  gravity.eta = 0.f;

  __shared__ struct Epj_GravityGPU jpsh[3][N_THREAD_GPU];
      
  if(1 == nwalk_in_block) {
    GravityForceKernel_1walk(jpsh[0], ip, id_walk, ij_disp, epj, &gravity);
  } else if(2 == nwalk_in_block) {
    const int iwalk0 = epi[t_head].id_walk;
    const int iwalk1 = epi[t_tail].id_walk;
    GravityForceKernel_2walk(jpsh, ip, id_walk, ij_disp, epj, &gravity, iwalk0, iwalk1);
  } else if(3 == nwalk_in_block) {
    // const int iwalk0 = epi[t_head].id_walk;
    // const int iwalk1 = iwalk0 + 1;
    // const int iwalk2 = iwalk0 + 2;
    // GravityForceKernel_3walk(jpsh, ip, id_walk, ij_disp, epj, &gravity, iwalk0, iwalk1, iwalk2);
    GravityForceKernel_multiwalk(ip, id_walk, ij_disp, epj, &gravity);
  } else {
    GravityForceKernel_multiwalk(ip, id_walk, ij_disp, epj, &gravity);
  }
  
  //GravityForceKernel_multiwalk(ip, id_walk, ij_disp, epj, &gravity);
  
  dev_gravity[tid] = gravity;

}

static cudaPointer<Epi_GravityGPU>   dev_epi;
static cudaPointer<Epj_GravityGPU>   dev_epj;
static cudaPointer<Spj_GravityGPU>   dev_spj;
static cudaPointer<GravityGPU> dev_gravity;
static cudaPointer<int2>     ij_disp;
static bool init_call = true;

PS::S32 DispatchKernel_Gravity(const PS::S32 tag,
			       const PS::S32 n_walk,
			       const GravityEPI **epi,
			       const PS::S32 *nip,
			       const GravityEPJ **epj,
			       const PS::S32 *njp,
			       const PS::GravitySPJ **spj,
			       const PS::S32 *nsp)
{
  assert(n_walk <= N_WALK_LIMIT);

  if(init_call){
    dev_epi    .allocate(NI_LIMIT);
    dev_epj    .allocate(NJ_LIMIT);
    dev_gravity.allocate(NI_LIMIT);
    ij_disp    .allocate(N_WALK_LIMIT+2);
    init_call = false;
  }

  ij_disp[0].x = 0;
  ij_disp[0].y = 0;
  for(int k=0; k<n_walk; k++){
    ij_disp[k+1].x = ij_disp[k].x + nip[k];
    ij_disp[k+1].y = ij_disp[k].y + (njp[k] + nsp[k]);
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
  
  int iw, i, j, num_i, num_j, num_sj;
#pragma omp parallel for private(i,j, num_i, num_j, num_sj) reduction(+:ni_tot, nj_tot)
  for(iw=0; iw<n_walk; iw++){
    for(i=0; i<nip[iw]; i++){
      num_i = ij_disp[iw].x + i;
      dev_epi[num_i].pos.x = epi[iw][i].pos.x;
      dev_epi[num_i].pos.y = epi[iw][i].pos.y;
      dev_epi[num_i].pos.z = epi[iw][i].pos.z;
      dev_epi[num_i].eps2 = epi[iw][i].eps2;
      dev_epi[num_i].id_walk = iw;
      ni_tot++;
    }
    for(j=0; j<njp[iw]; j++){
      num_j = ij_disp[iw].y + j;
      dev_epj[num_j].pos.x = epj[iw][j].pos.x;
      dev_epj[num_j].pos.y = epj[iw][j].pos.y;
      dev_epj[num_j].pos.z = epj[iw][j].pos.z;
      dev_epj[num_j].mass = epj[iw][j].mass;
      dev_epj[num_j].eps2 = epj[iw][j].eps2;
      nj_tot++;
    }
    for(j=0; j<nsp[iw]; j++){
      num_sj = ij_disp[iw].y + njp[iw] + j;
      dev_epj[num_sj].pos.x = spj[iw][j].pos.x;
      dev_epj[num_sj].pos.y = spj[iw][j].pos.y;
      dev_epj[num_sj].pos.z = spj[iw][j].pos.z;
      dev_epj[num_sj].mass = spj[iw][j].mass;
      dev_epj[num_sj].eps2 = spj[iw][j].eps2;
      nj_tot++;
    }
  }

  for(int i=ni_tot; i<ni_tot_reg; i++){
    dev_epi[i].id_walk = n_walk;
  }

  dev_epi.htod(ni_tot_reg);
  dev_epj.htod(nj_tot);

  int nblocks  = ni_tot_reg / N_THREAD_GPU;
  int nthreads = N_THREAD_GPU;

  GravityForceKernel <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_gravity);

  return 0;
}

PS::S32 RetrieveKernel_Gravity(const PS::S32 tag,
			       const PS::S32 n_walk,
			       const PS::S32 *ni,
			       Gravity **gravity)
{
  int ni_tot = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }
  dev_gravity.dtoh(ni_tot);

  int n_cnt = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<ni[iw]; i++){ 
      gravity[iw][i].acc.x = dev_gravity[n_cnt].acc.x;
      gravity[iw][i].acc.y = dev_gravity[n_cnt].acc.y;
      gravity[iw][i].acc.z = dev_gravity[n_cnt].acc.z;
      gravity[iw][i].pot = dev_gravity[n_cnt].pot;
      gravity[iw][i].eta = dev_gravity[n_cnt].eta;
      n_cnt++;
    }
  }

  return 0;
}
