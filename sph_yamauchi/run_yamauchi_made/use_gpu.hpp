#pragma once

class DensityEPI;
class DensityEPJ;
class Density;
class HydroEPI;
class HydroEPJ;
class Hydro;
class GravityEPI;
class GravityEPJ;
namespace ParticleSimulator{
class GravitySPJ;
};
class Gravity;

void setup_gpu(int my_proc);

void reset_gpu(int my_proc);

enum{
  N_THREAD_GPU = 32,
  N_WALK_LIMIT = 1000,
  NI_LIMIT = N_WALK_LIMIT*1000,
  NJ_LIMIT = N_WALK_LIMIT*10000,
  NumberOfDimension = 3
};

PS::S32 DispatchKernel_Dens(const PS::S32 tag,
			    const PS::S32 n_walk,
			    const DensityEPI **epi,
			    const PS::S32 *nip,
			    const DensityEPJ **epj,
			    const PS::S32 *njp);

PS::S32 RetrieveKernel_Dens(const PS::S32 tag,
			    const PS::S32 n_walk,
			    const PS::S32 *ni,
			    Density **density);

PS::S32 DispatchKernel_HydroForce(const PS::S32 tag,
				  const PS::S32 n_walk,
				  const HydroEPI **epi,
				  const PS::S32 *nip,
				  const HydroEPJ **epj,
				  const PS::S32 *njp);

PS::S32 RetrieveKernel_HydroForce(const PS::S32 tag,
				  const PS::S32 n_walk,
				  const PS::S32 *ni,
				  Hydro **hydro);

PS::S32 DispatchKernel_Gravity(const PS::S32 tag,
			       const PS::S32 n_walk,
			       const GravityEPI **epi,
			       const PS::S32 *nip,
			       const GravityEPJ **epj,
			       const PS::S32 *njp,
			       const PS::GravitySPJ **spj,
			       const PS::S32 *nsp);

PS::S32 RetrieveKernel_Gravity(const PS::S32 tag,
			       const PS::S32 n_walk,
			       const PS::S32 *ni,
			       Gravity **gravity);

