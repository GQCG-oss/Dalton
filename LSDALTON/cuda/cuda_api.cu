#include <cuda.h>
#include <cuda_runtime.h>

extern "C" {
void get_dev_mem(size_t& total, size_t& free) 
{
    cuMemGetInfo(&free, &total);
}
}
