#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <threads.h>
#include <complex.h>

int threads = 0;
int size = 0;
int degree_global = 0;

float epsilon=0.0001;
double epsilon2=0.000001;
long int upper_limit=10000000000;


complex double roots[10][10];
complex double z_global[size][size];



typedef struct {
  //data
  float complex **z;
  int **iter;
  int **attr;
  //thread controll
  int ib;
  int istep;
  int sz;
  int tx;
  //mtx_ cnd_
  mtx_t *mtx;
  cnd_t *cnd;
} thrd_info_t;



int main_thrd( void *args ){
  //structure
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  //data
  float complex **z = thrd_info->z;
  int **attr = thrd_info->attr;
  int **iter = thrd_info->iter;
  //thread control
  const int ib = thrd_info->ib;
  const int istep = thrd_info->istep;
  const int sz = thrd_info->sz;
  const int tx = thrd_info->tx;
  //mutex + cnd
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;

  //
  int degree = degree_global;


  // one row
  for ( int ix = ib; ix < sz; ix += istep ) {
    // We allocate the rows of the result before computing, and free them in another thread.
    float complex *input = z_global[ix];
    float complex *zix = z[ix];
    int *attr_ix = attr[ix];
    int *iter_ix = iter[ix];
    
    for ( int col = 0; col < sz; ++col ) {
      int conv;
      // newton's iteration
      for ( conv = 0, attr_ix[col] = -1 ; conv<128; ++conv ) { 
	if ( creal(zix[col])*creal(zix[col])+cimag(zix[col])*cimag(zix[col]) <= 1e-6 ) {
	  attr_ix[col] = -1;
	  break;
	}
	if ( fabs(creal(zix[col])) > upper_limit || fabs(cimag(zix[col])) >upper_limit ) {
	  attr_ix[col] = -1;
	  break;
	}
	for ( int ix_root=0; ix_root < degree; ix_root++ ){
	  if ( fabs(zix[col]-roots[degree-1][ix_root]) < 1e-3 ) {
	    attr_ix[col] = ix_root;
	    break;
	  }
	}
	if ( attr_ix[col] != -1 )
	  break;

	
	// computation
	switch ( degree ) {
	case 1:
	  //STATEMENTS FOR DEGREE 1;
	  zix[col]=zix[col]-(zix[col]-1);
	  break;
	case 2:
	  //STATEMENTS FOR DEGREE 2;	
	  zix[col]=zix[col]-((zix[col]*zix[col]-1)/(2*zix[col]));
	  break;
	case 3:
	  //STATEMENTS FOR DEGREE 3;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]-1)/(3*zix[col]*zix[col]));
	  break;
	case 4:
	  //STATEMENTS FOR DEGREE 4;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]*zix[col]-1)/(4*zix[col]*zix[col]*zix[col]));
	  break;
	case 5:
	  //STATEMENTS FOR DEGREE 5;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]*zix[col]*zix[col]-1)/(5*zix[col]*zix[col]*zix[col]*zix[col]));
	  break;
	case 6:
	  //STATEMENTS FOR DEGREE 6;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]-1)/(6*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]));
	  break;
	case 7:
	  //STATEMENTS FOR DEGREE 7;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]-1)/(7*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]));
	  break;
	case 8:
	  //STATEMENTS FOR DEGREE 8;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]-1)/(8*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]));
	  break;
	case 9:
	  //STATEMENTS FOR DEGREE 9;
	  zix[col]=zix[col]-((zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]-1)/(9*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]));
	  break;
	case 10:
	  //STATEMENTS FOR DEGREE 10;
	  zix[col]=zix[col]-(( zix[col] * zix[col] * zix[col] * zix[col] * zix[col] * zix[col] * zix[col] * zix[col] * zix[col] *zix[col]-1)/    (10  *zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]*zix[col]));
	  break;
	  // insert further cases

	default:
	  fprintf(stderr, "unexpected degree\n");
	  exit(1);
	}

	
      }// newton iter loop end
      
      iter_ix[col]=conv;

    }// column loop end

    
    mtx_lock(mtx);
    iter[ix] = iter_ix;
    attr[ix] = attr_ix;
    printf("Thread %d : ", tx);
    for (int i = 0; i < sz; i++) {
      printf("%lf+j%lf  ",creal(input[i]),cimag(input[i]));
      printf("%lf+j%lf %d %d ", creal(zix[i]), cimag(zix[i]), attr_ix[i] , iter_ix[i]);

      /* printf("attr : %d \n", attr_ix[i]); */
      /* printf("iter : %d \n", iter_ix[i]); */
    }
    printf("\n");
    //status[tx].val = ix + istep;
    mtx_unlock(mtx);
    cnd_signal(cnd);

    // In order to illustrate thrd_sleep and to force more synchronization
    // points, we sleep after each line for one micro seconds.
    thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=1000}, NULL);
    
  }// row in thread

  return 0;
}




int main(int argc, char *argv[]){
  // roots for x - 1
  roots[0][0] = 1 + 0 * I;
  // roots for x^2 - 1
  roots[1][0] = 1 + 0 * I;
  roots[1][1] = -1 + 0 * I;
  // roots for x^3 - 1
  roots[2][0] = 1 + 0 * I;
  roots[2][1] = -0.5 + 0.86603 * I;
  roots[2][2] = -0.5 - 0.86606 * I;
  // roots for x^4 - 1
  roots[3][0] = 1 + 0 * I;
  roots[3][1] = 0 + 1 * I;
  roots[3][2] = -1 + 0 * I;
  roots[3][3] = 0 - 1 * I;
  // roots for x^5 - 1
  roots[4][0] = 1 + 0 * I;
  roots[4][1] = 0.309017 + 0.951057 * I;
  roots[4][2] = -0.809017 + 0.587785 * I;
  roots[4][3] = -0.809017 - 0.587785 * I;
  roots[4][4] = 0.309017 - 0.951057 * I;
  // roots for x^6 - 1
  roots[5][0] = 1 + 0 * I;
  roots[5][1] = 0.5 + 0.866025 * I;
  roots[5][2] = -0.5 + 0.866025 * I;
  roots[5][3] = -1 - 0 * I;
  roots[5][4] = -0.5 - 0.866025 * I;
  roots[5][5] = 0.5 - 0.866025 * I;
  // roots for x^7 - 1
  roots[6][0] = 1 + 0 * I;
  roots[6][1] = 0.62349 + 0.781831 * I;
  roots[6][2] = -0.222521 + 0.974928 * I;
  roots[6][3] = -0.900969 + 0.433884 * I;
  roots[6][4] = -0.900969 - 0.433884 * I;
  roots[6][5] = -0.222521 - 0.974928 * I;
  roots[6][6] = 0.62349 - 0.781831 * I;
  // roots for x^8 - 1
  roots[7][0] = 1 + 0 * I;
  roots[7][1] = 0.707107 + 0.707107 * I;
  roots[7][2] = 0 + 1 * I;
  roots[7][3] = -0.707107 + 0.707107 * I;
  roots[7][4] = -1 + 0 * I;
  roots[7][5] = -0.707107 - 0.707107 * I;
  roots[7][6] = 0 - 1 * I;
  roots[7][7] = 0.707107 - 0.707107 * I;
  // roots for x^9 - 1
  roots[8][0] = 1 + 0 * I;
  roots[8][1] = 0.766044 + 0.642788 * I;
  roots[8][2] = 0.173648 + 0.984808 * I;
  roots[8][3] = -0.5 + 0.866025 * I;
  roots[8][4] = -0.939693 + 0.34202 * I;
  roots[8][5] = -0.939693 - 0.34202 * I;
  roots[8][6] = -0.5 - 0.866025 * I;
  roots[8][7] = 0.173648 - 0.984808 * I;
  roots[8][8] = 0.766044 - 0.642788 * I;

  FILE *fp = fopen("numbers.txt", "w");
  if (fp == NULL) {
    perror("Error opening the file");
    return 1;
  }


  /*for (int i = 0; i < size; i++) {
    real = -2.0 + i * step;

    for (int j = 0; j < size; j++) {
    imag = -2.0 + j * step;

    fprintf(fp, "Complex Number: %lf + j%lf\n", real, imag);
    }
    }

    fclose(fp);*/ 

  for (int i = 1; i < 4; i++) {
    if (sscanf(argv[i], "-t%d", &threads) == 1) {
      //printf("Extracted value: %d\n", threads);
    } else if (sscanf(argv[i], "-l%d", &size) == 1) {
      //printf("Extracted value: %d\n", size);
    }
    else if (sscanf(argv[i], "-d%d", &degree_global) == 1) {
      //printf("Extracted value: %d\n", degree_global);
    }
    else {
      printf("Invalid input.\n");
      return 1;
    }
  }

  printf("\nthreads : %d\nsize : %d\ndegree_global : %d\n\n", threads, size, degree_global);

  const int sz = size;
   
  /* float **re = (float**) malloc(sz*sizeof(float*)); */
  /* float **im = (float**) malloc(sz*sizeof(float*)); */
  float complex **z = (float complex**) malloc(sz*sizeof(float complex * ));
  /* float **v = (float**) malloc(sz*sizeof(float*)); */
  /* float **w = (float**) malloc(sz*sizeof(float*)); */
  
  int **attr = (int**) malloc(sz*sizeof(int*));
  int **iter = (int**) malloc(sz*sizeof(int*));
  
  /* float *reentries = (float*) malloc(sz*sz*sizeof(float)); */
  /* float *imentries = (float*) malloc(sz*sz*sizeof(float)); */
  float complex *zentries = (float complex*) malloc(sz*sz*sizeof(float complex));
  int *attr_entries = (int*) malloc(sz*sz*sizeof(int));
  int *iter_entries = (int*) malloc(sz*sz*sizeof(int));
  // The entries of attr and iter will be allocated in the computation threads are freed in the check thread.

  for ( int ix = 0, jx = 0; ix < sz; ++ix, jx += sz ){
    attr[ix] = attr_entries + jx;
    iter[ix] = iter_entries + jx;
    z[ix] = zentries + jx;
  }
  for ( int ix = 0; ix < sz*sz; ++ix ){
    attr_entries[ix] = 0;
    iter_entries[ix] = 0;
    zentries[ix] = 0;
  }


  
  double step = 4.0 / (sz - 1); 
  double real, imag;
  /* for (int i = 0; i < sz; i++) { */
  /*   //real = -2.0 + i * step; */
  /*   for (int j = 0; j < sz; j++) { */
  /*     //imag = -2.0 + j * step; */
  /*     reentries[i * sz + j] = -2.0 + i * step; // not sure  */
  /* 						    imentries[i * sz + j] = 2.0 - j * step; // about this */
  /*     //fprintf(fp, "Complex Number: %lf + j%lf\n", reentries[i], imentries[j]); */
  /*   } */
  /* } */

  for (int i = 0; i < sz; i++) {
    for (int j = 0; j < sz; j++) {
      double real = -2.0 + i * step;
      double imag = 2.0 - j * step;
      z[i][j] = real+_Complex_I*imag;//CMPLX(real, imag);
      z_global[i][j] = real+_Complex_I*imag;//CMPLX(real, imag);
      // Print the complex number using creal and cimag
      fprintf(fp, "Complex Number: %lf + j%lf\n", creal(z[i][j]), cimag(z[i][j]));
    }
  }


  //thread create
  const int nthrds = threads;
  thrd_t thrds[nthrds];
  //structure create
  thrd_info_t thrds_info[nthrds];

  
  /* thrd_t thrd_check; */
  /* thrd_info_check_t thrd_info_check; */

  mtx_t mtx;
  mtx_init(&mtx, mtx_plain);

  cnd_t cnd;
  cnd_init(&cnd);

  //  int_padded status[nthrds];

  for ( int tx = 0; tx < nthrds; ++tx ) {
    /* thrds_info[tx].re = (const float**) re; */
    /* thrds_info[tx].im = (const float**) im; */
    thrds_info[tx].z = (float complex**) z;
    thrds_info[tx].attr = attr;
    thrds_info[tx].iter = iter;

    //
    thrds_info[tx].ib = tx;
    thrds_info[tx].istep = nthrds;
    thrds_info[tx].sz = sz;
    thrds_info[tx].tx = tx;
    thrds_info[tx].mtx = &mtx;
    thrds_info[tx].cnd = &cnd;
    /* thrds_info[tx].status = status; */
    /* status[tx].val = -1; */

    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
    thrd_detach(thrds[tx]);
  }

  /* { */
  /*   thrd_info_check.re = (const float**) re; */
  /*   thrd_info_check.im = (const float**) im; */
  /*   thrd_info_check.z = (const float complex**) z; */
  /*   thrd_info_check.v = (const float**) v; */
  /*   thrd_info_check.w = w; */
  /*   thrd_info_check.sz = sz; */
  /*   thrd_info_check.nthrds = nthrds; */
  /*   thrd_info_check.mtx = &mtx; */
  /*   thrd_info_check.cnd = &cnd; */
  /*   // It is important that we have initialize status in the previous for-loop, */
  /*   // since it will be consumed by the check threads. */
  /*   thrd_info_check.status = status; */

  /*   int r = thrd_create(&thrd_check, main_thrd_check, (void*) (&thrd_info_check)); */
  /*   if ( r != thrd_success ) { */
  /*     fprintf(stderr, "failed to create thread\n"); */
  /*     exit(1); */
  /*   } */
  /* } */

  /* { */
  /*   int r; */
  /*   thrd_join(thrd_check, &r); */
  /* } */

  free(zentries);
  free(attr_entries);
  free(iter_entries);
  free(z);
  free(iter);
  free(attr);

  mtx_destroy(&mtx);
  cnd_destroy(&cnd);

  fclose(fp);

  return 0;
}
 
