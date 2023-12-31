#include <stdlib.h>
#include <stdio.h>
#include <threads.h>


int szloc;

int
main_thrd(
    void *args
    )
{
  //now the pointer passing is the 2D array's addr
  double *v = ((double**)args)[0];
  //so as you can see, this is for calculation of one (the 1st)  row in args[][]
  double sum = 0.;
  for ( int ix = 0; ix < szloc; ++ix )
    sum += v[ix];
  *((double**)args)[1] = sum;

  
  // Memory can be freed in a different thread than where it was allocated, but
  // it has to be freed exactly once.
  free(args);

  return 0;
}

int
main()
{
  const int nthrds = 8;

  //------------------------------//
  const int sz = 1000;
  double *v = (double*) malloc(sz*nthrds*sizeof(double));

  for ( int ix = 0; ix < sz; ++ix )
    v[ix] = ix+1;

  thrd_t thrds[nthrds];
  
  double sums[nthrds];
  szloc = sz; 
  //------------------------------//
  
  for ( int tx = 0; tx < nthrds; ++tx ) {
    // To pass arguments to threads we cannot use array allocation
    // double args[2];
    // inside of the for loop, since this would be freed one the scope of the
    // for-loop's body ends. It is possible to use array allocation places
    // outside of the for-loop.
    double** args = (double**) malloc(2*sizeof(double*));
    args[0] = v + tx*sz;
    args[1] = sums + tx;




    //create threads
    int r;
    r = thrd_create(thrds+tx, main_thrd, (void*)args);
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }
  //fork threads
  for ( int tx = 0; tx < nthrds; ++tx ) {
    int r;
    thrd_join(thrds[tx], &r);
  }




  //
  double sum = sums[0];
  for ( int tx = 1; tx < nthrds; ++tx )
    sum += sums[tx];

  printf("sum: %f\n", sum);

  return 0;
}
