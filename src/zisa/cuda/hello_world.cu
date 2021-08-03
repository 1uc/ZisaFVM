#include <zisa/cuda/hello_world.hpp>

#include <stdio.h>

namespace zisa {

// TODO Remove once there is real CUDA code to compile.
__global__ void hello_world_kernel() { printf("hello world.\n"); }

void hello_world() {
  hello_world_kernel<<<1, 1>>>();
  cudaDeviceSynchronize();
}
}
