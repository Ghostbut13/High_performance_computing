#include <stdlib.h>
#include <stdio.h>
#include <threads.h>
//------------------------------
// To gain flexibility one can define a structure that groups further variables
// and then submit to the created thread a pointer to one or several of these
// structures. This is as powerful as the usual passing of arguments to
// functions.
//------------------------------
//so i mean: since main_thread only permits one parameter( and only void*) passing
//, we can create a structure to pass multi parameters
typedef struct {
  const double *v;
  int ib;
  int ie;
  double *sum;
  mtx_t *mtx;
} thrd_info_t;
// The name thrd_info_t of any standard, so any other name can be used here.

int
main_thrd(
    void *args
    )
{
  //in thread, we initialize a structure to carry the args's structure
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  const double* v = thrd_info->v;
  const int ib = thrd_info->ib;
  const int ie = thrd_info->ie;

  double sum = 0.;
  for ( int ix = ib; ix < ie; ++ix )
    sum += v[ix];

  // When writing to memory that other threads might be access as well, it has
  // to be "protected". A typical technique is to lock a mutex, which can be
  // locked by at most one thread at a time.
  mtx_lock(thrd_info->mtx);
  *thrd_info->sum += sum;
  mtx_unlock(thrd_info->mtx);

  return 0;
}

int
main()
{
  const int sz = 1783;
  double *v = (double*) malloc(sz*sizeof(double));

  for ( int ix = 0; ix < sz; ++ix )
    v[ix] = ix+1;

  const int nthrds = 8;
  thrd_t thrds[nthrds];
  // Here we array allocate what later will be used as arguments to the thread
  // main functions.
  //------------------------------
  //initialize each thread_info for each thread
  thrd_info_t thrds_info[nthrds];

  //********************create a group thrds array******************************
  //**************************************************
  const int szloc = sz / nthrds;//this is interval (step)
  for ( int tx = 0, ib = 0; tx < nthrds; ++tx, ib += szloc ) {
    thrds_info[tx].v = v;    // v is vector input( we want to use it)
    thrds_info[tx].ib = ib;  // ib is begin of this v (we cut 'v' into pieces, different pieces have different begin location)
    //acutully, if sz cannot be perfectly divided by nthreads like, 97/4, we need this
    thrds_info[tx].ie = tx != nthrds - 1 ? ib + szloc : sz;// ib is end of this v (we cut 'v' into pieces, different pieces have different end location)
    //if not, i prefer: thrds_info[tx].ie = ib+szloc;

    //*********** KEY KEY KEY *******************//
    //******************************//
    // Both the sum and the mutex variable point to the same memory in all
    // threads, and are therefore shared memory.
    double sum = 0.;
    // A mutex is a construct to synchronize threads. A mutex can be "locked" by
    // at most one thread at a time. As allocated memory must be freed, mutices
    // must be destroyed after use.
    mtx_t mtx;
    mtx_init(&mtx, mtx_plain);
    thrds_info[tx].sum = &sum;
    thrds_info[tx].mtx = &mtx;
    //******************************//

    
    //CREATE threads
    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }
  //**************************************************


  
  for ( int tx = 0; tx < nthrds; ++tx ) {
    int r;
    thrd_join(thrds[tx], &r);
  }


  // destroy mutex
  mtx_destroy(&mtx);


  
  printf("sum diff: %f\n", sum - sz * (double)(sz+1) / 2);

  return 0;
}
