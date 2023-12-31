#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <threads.h>

int
main_thrd(
    void *args
    )
{
  if ( args != NULL )
    return 1;

  return 0;
}

int
main()
{
  const int nthrds = 8;

  // On creation, a thread initialized a corresponing structure of type thrd_t.
  // This structure can later be used to wait for the threads termination. Very
  // little control beyond this is possible throught the thrd_t structure.
  thrd_t thrds[nthrds];
  
  for ( int tx = 0; tx < nthrds; ++tx ) {
    int r;
    if ( tx & 1 ) // when tx ==0 , i,e., it is number of thread; and it is 0
      r = thrd_create(thrds+tx, main_thrd, (void*)main);// we will regard 0 thread as main thread
    else
      r = thrd_create(thrds+tx, main_thrd, NULL);
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
    // anyway, i prefer just : r = thrd_create(thrds+tx, main_thrd, NULL);
  }

  for ( int tx = 0; tx < nthrds; ++tx ) {
    int r;
    thrd_join(thrds[tx], &r);
    if ( tx & 1 )
      assert( r == 1 );
    else 
      assert( r == 0 );
  }

  return 0;
}
