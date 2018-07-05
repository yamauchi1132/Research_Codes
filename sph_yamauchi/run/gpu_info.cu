#include <stdio.h>
#include <cuda.h>
#include <particle_simulator.hpp>
#include "use_gpu.hpp"

void setup_gpu(int my_proc)
{
  if(my_proc == 0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "###################\n");
    fprintf(stderr, "### Using GPU ! ###\n");
    fprintf(stderr, "###################\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "Setting up GPU...\n");
  }

  //getting the number of GPU per node
  int ngpus;
  cudaGetDeviceCount(&ngpus);

  //set up device for each process
  int dev = my_proc % ngpus;
  //fprintf(stderr, "%d %d %d test\n", dev, ngpus, my_proc);
  cudaSetDevice(dev);

  //getting information of gpu
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  if(my_proc == 0) {
    fprintf(stderr, "%d GPU per node available\n", ngpus);
    fprintf(stderr, "GPU name : %s\n\n", deviceProp.name);
  }
}

void reset_gpu(int my_proc)
{
  if(my_proc == 0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Resetting GPU...\n");
  }
  cudaDeviceReset();
}
