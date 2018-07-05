#include <particle_simulator.hpp>
#include <omp.h>
#include "cuda_pointer.h"
#include "helper_cuda.h"
#include "DF.hpp"
#include "use_gpu.hpp"

__constant__ DF ceff1_c;
/*
  enum{
  N_THREAD_GPU = 32,
  N_WALK_LIMIT = 1000,
  NI_LIMIT = N_WALK_LIMIT*1000,
  NJ_LIMIT = N_WALK_LIMIT*10000,
  };
*/
class HydroEPI{
public:
  PS::S64 id;
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64    ksr;
  PS::F64    hinv;
  PS::F64    pres;
  PS::F64    presc;
  PS::F64    bswt;
  PS::F64    dens;
  PS::F64    vsnd;
  PS::F64    alph;
  PS::F64    alphu;
  PS::F64    thrm;
  PS::F64    eta;
};

class HydroEPJ{
public:
  PS::S64 id;
  PS::F64    mass;
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64    ksr;
  PS::F64    hinv;
  PS::F64    pres;
  PS::F64    presc;
  PS::F64    bswt;
  PS::F64    dens;
  PS::F64    vsnd;
  PS::F64    alph;
  PS::F64    alphu;
  PS::F64    thrm;
  PS::F64    eta;
};

class Hydro{
public:
  PS::F64vec acch;
  PS::F64vec accg;
  PS::F64    udot;
  PS::F64    vsmx;
  PS::F64    diffu;
};

struct Epi_HydroGPU{
  long long int id;
  DF3 pos;
  DF3 vel;
  DF hinv;
  DF pres;
  DF presc;
  DF bswt;
  DF dens;
  DF vsnd;
  DF alph;
  DF alphu;
  DF thrm;
  DF eta;
  int id_walk;
};

struct Epj_HydroGPU{
  long long int id;
  DF3 pos;
  DF3 vel;
  DF hinv;
  DF pres;
  DF presc;
  DF bswt;
  DF dens;
  DF vsnd;
  DF alph;
  DF alphu;
  DF thrm;
  DF eta;
  DF mass;
};

struct HydroGPU{
  DF3 acc;
  DF udot;
  DF vsmx;
  DF diffu;
  DF3 accg;
};

inline __device__ DF kernel1st(const DF q)
{
  // const DF ceff1 = +1.336901521971920914e+01;
  const DF qmin  = ((1. - q > 0.) ? 1. - q : 0.);
  const DF qmin2 = qmin  * qmin;
  const DF qmin3 = qmin  * qmin2;
  const DF qmin4 = qmin2 * qmin2;
  return ceff1_c * (qmin4 - qmin3 * (1. + 4. * q));
}

inline __device__ void HydroForce_calc(const struct Epi_HydroGPU ip,
				       const struct Epj_HydroGPU jp,
				       DF3 *acc,
				       DF3 *accg,
				       DF *udot,
				       DF *vsmx,
				       DF *diffu,
				       const DF hi4_i)
{
  const DF dpx_ij = ip.pos.x - jp.pos.x;
  const DF dpy_ij = ip.pos.y - jp.pos.y;
  const DF dpz_ij = ip.pos.z - jp.pos.z;
  const DF dvx_ij = ip.vel.x - jp.vel.x;
  const DF dvy_ij = ip.vel.y - jp.vel.y;
  const DF dvz_ij = ip.vel.z - jp.vel.z;

  const DF r2_ij = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij;
  //const DF ri_ij = ((ip.id != jp.id) ? 1. / sqrt(r2_ij) : 0.);
  const DF ri_ij = ((ip.id != jp.id) ? rsqrtf(r2_ij) : 0.);
  const DF r1_ij = r2_ij * ri_ij;
  const DF q_i   = r1_ij * ip.hinv;
  const DF q_j   = r1_ij * jp.hinv;

  const DF hi4_j = jp.hinv * jp.hinv * jp.hinv * jp.hinv;
  const DF dw_i  = hi4_i * kernel1st(q_i);
  const DF dw_j  = hi4_j * kernel1st(q_j);
  const DF ka_ij = (dw_i + dw_j) * jp.mass;

  const DF rv_ij  = dpx_ij * dvx_ij + dpy_ij * dvy_ij + dpz_ij * dvz_ij;
  const DF w_ij   = rv_ij * ri_ij;
  const DF w0_ij  = ((w_ij < 0.) ? w_ij : 0.);
  const DF vs_ij  = ip.vsnd + jp.vsnd - 3. * w0_ij;
  const DF rhi_ij = 1. / (ip.dens + jp.dens);
  const DF av0_ij = (ip.bswt + jp.bswt) * (ip.alph + jp.alph) * vs_ij * w0_ij;
  
  *vsmx  = ((*vsmx > vs_ij) ? *vsmx : vs_ij);

  const DF ta_ij = (ip.presc + jp.presc - 0.5 * av0_ij * rhi_ij) * ka_ij * ri_ij;
  acc->x -= ta_ij * dpx_ij;
  acc->y -= ta_ij * dpy_ij;
  acc->z -= ta_ij * dpz_ij;

  const DF vsu2_ij = fabs(ip.pres - jp.pres) * rhi_ij * 2.;
  const DF vsui_ij = ((vsu2_ij != 0.) ? 1. / sqrt(vsu2_ij) : 0.);
  const DF vsu_ij  = vsu2_ij * vsui_ij;
  const DF du_ij   = ip.thrm - jp.thrm;

  *udot += ka_ij * (ip.presc * w_ij - rhi_ij * (0.25 * av0_ij * w0_ij - (ip.alphu + jp.alphu) * vsu_ij * du_ij));

  *diffu += rhi_ij * du_ij * ka_ij * ri_ij;

  const DF dg_ij = jp.mass * ri_ij * (ip.eta * dw_i + jp.eta * dw_j);
  accg->x += dg_ij * dpx_ij;
  accg->y += dg_ij * dpy_ij;
  accg->z += dg_ij * dpz_ij;
}

__device__ void HydroForceKernel_1walk(struct Epj_HydroGPU *jpsh,
				       const struct Epi_HydroGPU ip,
				       const Epj_HydroGPU *epj,
				       const int id_walk,
				       const int2 *ij_disp,
				       DF3 *acc,
				       DF3 *accg,
				       DF *udot,
				       DF *vsmx,
				       DF *diffu)
{
  const int tid = threadIdx.x;
  const int j_head = ij_disp[id_walk  ].y;
  const int j_tail = ij_disp[id_walk+1].y;

  const DF hi4_i = ip.hinv * ip.hinv * ip.hinv * ip.hinv;

  for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
    jpsh[tid] = ((Epj_HydroGPU *)(epj + j))[tid];
    if(j_tail-j < N_THREAD_GPU){
      for(int jj=0; jj<j_tail-j; jj++){
	HydroForce_calc(ip, jpsh[jj], acc, accg, udot, vsmx, diffu, hi4_i);
      }
    }else{
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	HydroForce_calc(ip, jpsh[jj], acc, accg, udot, vsmx, diffu, hi4_i);
      }
    }
  }

  acc->x *= 0.5;
  acc->y *= 0.5;
  acc->z *= 0.5;
  accg->x *= 0.5;
  accg->y *= 0.5;
  accg->z *= 0.5;
  *udot *= 0.5;
  *diffu *= 2.;

}

__device__ void HydroForceKernel_2walk(struct Epj_HydroGPU (*jpsh)[N_THREAD_GPU],
				       const struct Epi_HydroGPU ip,
				       const Epj_HydroGPU *epj,
				       const int id_walk,
				       const int2 *ij_disp,
				       DF3 *acc,
				       DF3 *accg,
				       DF *udot,
				       DF *vsmx,
				       DF *diffu,
				       int iwalk0,
				       int iwalk1)
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
 
  const DF hi4_i = ip.hinv * ip.hinv * ip.hinv * ip.hinv;

  for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_HydroGPU *)(epj + jbeg0 + j))[tid];
    jpsh[1][tid] = ((Epj_HydroGPU *)(epj + jbeg1 + j))[tid];

    if(nj_shorter-j < N_THREAD_GPU){
      for(int jj=0; jj<nj_shorter-j; jj++){	
	HydroForce_calc(ip, jpsh[mywalk][jj], acc, accg, udot, vsmx, diffu, hi4_i);
      }
    }else {
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	HydroForce_calc(ip, jpsh[mywalk][jj], acc, accg, udot, vsmx, diffu, hi4_i);
      }
    }
  }

  for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
    jpsh[0][tid] = ((Epj_HydroGPU *)(epj + jbeg_longer + j))[tid];
    int jrem = nj_longer - j;
    if(jrem < N_THREAD_GPU){
      for(int jj=0; jj<jrem; jj++){
	if(mywalk == walk_longer)
	  HydroForce_calc(ip, jpsh[0][jj], acc, accg, udot, vsmx, diffu, hi4_i);
      }
    }else {
#pragma unroll
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	if(mywalk == walk_longer)
	  HydroForce_calc(ip, jpsh[0][jj], acc, accg, udot, vsmx, diffu, hi4_i);
      }
    }
  }

  acc->x *= 0.5;
  acc->y *= 0.5;
  acc->z *= 0.5;
  accg->x *= 0.5;
  accg->y *= 0.5;
  accg->z *= 0.5;
  *udot *= 0.5;
  *diffu *= 2.;

}

__device__ void HydroForceKernel_multiwalk(const struct Epi_HydroGPU ip,
					   const Epj_HydroGPU *epj,
					   const int id_walk,
					   const int2 *ij_disp,
					   DF3 *acc,
					   DF3 *accg,
					   DF *udot,
					   DF *vsmx,
					   DF *diffu)
{ 
  const int j_head = ij_disp[id_walk  ].y;
  const int j_tail = ij_disp[id_walk+1].y;

  const DF hi4_i = ip.hinv * ip.hinv * ip.hinv * ip.hinv;
  for(int j=j_head; j<j_tail; j++){
    struct Epj_HydroGPU jp = epj[j];
    HydroForce_calc(ip, jp, acc, accg, udot, vsmx, diffu, hi4_i);
  }

  acc->x *= 0.5;
  acc->y *= 0.5;
  acc->z *= 0.5;
  accg->x *= 0.5;
  accg->y *= 0.5;
  accg->z *= 0.5;
  *udot *= 0.5;
  *diffu *= 2.;

}

__global__ void HydroForceKernel(const int2 *ij_disp,
				 const Epi_HydroGPU *epi,
				 const Epj_HydroGPU *epj,
				 HydroGPU *dev_hydro)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  int t_head = blockDim.x * blockIdx.x;
  int t_tail = t_head + N_THREAD_GPU - 1;
  int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);
#ifdef USE_FLOAT
  DF3 acc = make_float3(0.f, 0.f, 0.f);
  DF3 accg =  make_float3(0.f, 0.f, 0.f);
#endif
#ifdef USE_DOUBLE
  DF3 acc = make_double3(0.f, 0.f, 0.f);
  DF3 accg =  make_double3(0.f, 0.f, 0.f);
#endif

  DF udot=0.f;
  DF vsmx = 0.f;
  DF diffu = 0.f;

  int id_walk = epi[tid].id_walk;
  const struct Epi_HydroGPU ip = epi[tid];

  __shared__ struct Epj_HydroGPU jpsh[2][N_THREAD_GPU];
       
  if(1 == nwalk_in_block){
    HydroForceKernel_1walk(jpsh[0], ip, epj, id_walk, ij_disp, &acc, &accg, &udot, &vsmx, &diffu);       } else if(2 == nwalk_in_block){
    int iwalk0 = epi[t_head].id_walk;
    int iwalk1 = epi[t_tail].id_walk;
    HydroForceKernel_2walk(jpsh, ip, epj, id_walk, ij_disp, &acc, &accg, &udot, &vsmx, &diffu, iwalk0, iwalk1);
  } else{
    HydroForceKernel_multiwalk(ip, epj, id_walk, ij_disp, &acc, &accg, &udot, &vsmx, &diffu);
  }
  
  //HydroForceKernel_multiwalk(ip, epj, id_walk, ij_disp, &acc, &accg, &udot, &vsmx, &diffu);

  dev_hydro[tid].acc = acc;
  dev_hydro[tid].accg = accg;
  dev_hydro[tid].udot = udot;
  dev_hydro[tid].vsmx = vsmx;
  dev_hydro[tid].diffu = diffu;

}

static cudaPointer<Epi_HydroGPU>   dev_epi;
static cudaPointer<Epj_HydroGPU>   dev_epj;
static cudaPointer<HydroGPU> dev_hydro;
static cudaPointer<int2>     ij_disp;
static bool init_call = true;

PS::S32 DispatchKernel_HydroForce(const PS::S32 tag,
				  const PS::S32 n_walk,
				  const HydroEPI **epi,
				  const PS::S32 *nip,
				  const HydroEPJ **epj,
				  const PS::S32 *njp)
{
  assert(n_walk <= N_WALK_LIMIT);

  if(init_call){
    dev_epi   .allocate(NI_LIMIT);
    dev_epj   .allocate(NJ_LIMIT);
    dev_hydro .allocate(NI_LIMIT);
    ij_disp   .allocate(N_WALK_LIMIT+2);
    init_call = false;
  }

  ij_disp[0].x = 0;
  ij_disp[0].y = 0;
  for(int k=0; k<n_walk; k++){
    ij_disp[k+1].x = ij_disp[k].x + nip[k];
    ij_disp[k+1].y = ij_disp[k].y + njp[k];
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
    for(i=0; i<nip[iw]; i++){
       num_i = ij_disp[iw].x + i;
      dev_epi[num_i].id = epi[iw][i].id;
      dev_epi[num_i].pos.x = epi[iw][i].pos.x;
      dev_epi[num_i].pos.y = epi[iw][i].pos.y;
      dev_epi[num_i].pos.z = epi[iw][i].pos.z;
      dev_epi[num_i].vel.x = epi[iw][i].vel.x;
      dev_epi[num_i].vel.y = epi[iw][i].vel.y;
      dev_epi[num_i].vel.z = epi[iw][i].vel.z;
      dev_epi[num_i].hinv = epi[iw][i].hinv;
      dev_epi[num_i].pres = epi[iw][i].pres;
      dev_epi[num_i].presc = epi[iw][i].presc;
      dev_epi[num_i].bswt = epi[iw][i].bswt;
      dev_epi[num_i].dens = epi[iw][i].dens;
      dev_epi[num_i].vsnd = epi[iw][i].vsnd;
      dev_epi[num_i].alph = epi[iw][i].alph;
      dev_epi[num_i].alphu = epi[iw][i].alphu;
      dev_epi[num_i].thrm = epi[iw][i].thrm;
      dev_epi[num_i].eta = epi[iw][i].eta;
      dev_epi[num_i].id_walk = iw;
      ni_tot++;
    }
    for(j=0; j<njp[iw]; j++){
      num_j = ij_disp[iw].y + j;
      dev_epj[num_j].id = epj[iw][j].id;
      dev_epj[num_j].pos.x = epj[iw][j].pos.x;
      dev_epj[num_j].pos.y = epj[iw][j].pos.y;
      dev_epj[num_j].pos.z = epj[iw][j].pos.z;
      dev_epj[num_j].vel.x = epj[iw][j].vel.x;
      dev_epj[num_j].vel.y = epj[iw][j].vel.y;
      dev_epj[num_j].vel.z = epj[iw][j].vel.z;
      dev_epj[num_j].hinv = epj[iw][j].hinv;
      dev_epj[num_j].pres = epj[iw][j].pres;
      dev_epj[num_j].presc = epj[iw][j].presc;
      dev_epj[num_j].bswt = epj[iw][j].bswt;
      dev_epj[num_j].dens = epj[iw][j].dens;
      dev_epj[num_j].vsnd = epj[iw][j].vsnd;
      dev_epj[num_j].alph = epj[iw][j].alph;
      dev_epj[num_j].alphu = epj[iw][j].alphu; 
      dev_epj[num_j].thrm = epj[iw][j].thrm;
      dev_epj[num_j].eta = epj[iw][j].eta;
      dev_epj[num_j].mass = epj[iw][j].mass;
      nj_tot++;
    }
  }

  for(int i=ni_tot; i<ni_tot_reg; i++){
    dev_epi[i].id_walk = n_walk;
  }

  dev_epi.htod(ni_tot_reg);
  dev_epj.htod(nj_tot);

  const DF ceff1 = +1.336901521971920914e+01;

  cudaMemcpyToSymbol(ceff1_c, &ceff1, sizeof(DF));

  int nblocks  = ni_tot_reg / N_THREAD_GPU;
  int nthreads = N_THREAD_GPU;
  
  // cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual);  
  HydroForceKernel <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_hydro);

  return 0;
}

PS::S32 RetrieveKernel_HydroForce(const PS::S32 tag,
				  const PS::S32 n_walk,
				  const PS::S32 *ni,
				  Hydro         **hydro)
{
  int ni_tot = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }
  dev_hydro.dtoh(ni_tot);

  int n_cnt = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<ni[iw]; i++){ 
      hydro[iw][i].acch.x = dev_hydro[n_cnt].acc.x;
      hydro[iw][i].acch.y = dev_hydro[n_cnt].acc.y;
      hydro[iw][i].acch.z = dev_hydro[n_cnt].acc.z;
      hydro[iw][i].accg.x = dev_hydro[n_cnt].accg.x;
      hydro[iw][i].accg.y = dev_hydro[n_cnt].accg.y;
      hydro[iw][i].accg.z = dev_hydro[n_cnt].accg.z;
      hydro[iw][i].udot = dev_hydro[n_cnt].udot;
      hydro[iw][i].vsmx = dev_hydro[n_cnt].vsmx;
      hydro[iw][i].diffu = dev_hydro[n_cnt].diffu;
      n_cnt++;
    }
  }

  return 0;
}
