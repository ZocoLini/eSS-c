#include <method_module/structure_paralleltestbed.h>
#include <bbob/bbobStructures.h>

int parallel_asynchronous_DE(void *, void *, void  *,  void*(*fitnessfunction)(double*, void *), void *, long , double , int );
int parallel_synchronous_DE(void *, void *, void  *,  void*(*fitnessfunction)(double*, void *), void *, long , double , int );