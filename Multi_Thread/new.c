#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <threads.h>
#include <complex.h>
#include <time.h>

int threads = 0;
int size = 0;
int degree_global = 0;

float epsilon=1e-3;
double epsilon2=1e-6;
long int upper_limit=10000000000;

struct timespec start, finish;
double elapsed;

complex double roots[10][10];

typedef struct {
  float complex **z;
  int **iter;
  int **attr;
  int ib;
  int istep;
  int sz;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int *finish_flag;
} thrd_info_t;

typedef struct{
  int **iter;
  int **attr;
  int sz;
  mtx_t *mtx;
  cnd_t *cnd;
  FILE *colorful;
  FILE *black;
  int *finish_flag;
} thrd_write_info_t;





int write_thrd(void* args){
  const thrd_write_info_t *thrd_write_info = (thrd_write_info_t *) args;
  int **iter =  thrd_write_info->iter;
  int **attr =  thrd_write_info->attr;
  int sz = thrd_write_info->sz;
  mtx_t *mtx =  thrd_write_info->mtx;
  cnd_t *cnd =  thrd_write_info->cnd;
  FILE *colorful = thrd_write_info->colorful;
  FILE *black = thrd_write_info->black;
  int *finish_flag = thrd_write_info->finish_flag;

	clock_gettime(CLOCK_MONOTONIC, &start);
  
  for(int ix=0; ix<sz; ix++){
    int *attr_ix = attr[ix];
    int *iter_ix = iter[ix];
  
    mtx_lock(mtx);
    while(finish_flag[ix]==0){
      cnd_wait(cnd,mtx);
    }
        
    for(int col=0; col<sz; col++){
      fprintf(black, "%d %d %d ", iter_ix[col]*2, iter_ix[col]*2, iter_ix[col]*2);
      switch(attr_ix[col]){
      case 0:
	fprintf(colorful, "%d %d %d ",180 ,0 ,30 );
	break;
      case 1:
	fprintf(colorful, "%d %d %d ",0 ,180 ,30 );
	break;
      case 2:
	fprintf(colorful, "%d %d %d ",0 ,30 ,80 );
	break;
      case 3:
	fprintf(colorful, "%d %d %d ",0 ,190 ,180 );
	break;
      case 4:
	fprintf(colorful, "%d %d %d ",180 ,0 ,175 );
	break;
      case 5:
	fprintf(colorful, "%d %d %d ",180 ,255 ,0 );
	break;
      case 6:
	fprintf(colorful, "%d %d %d ",155 ,170 ,180 );
	break;
      case 7:
	fprintf(colorful, "%d %d %d ",70 ,50 ,0 );
	break;
      case 8:
	fprintf(colorful, "%d %d %d ",150 ,60 ,0 );
	break;
      case 9:
	fprintf(colorful, "%d %d %d ",0 ,150 ,60 );
	break;
      default://about attr_ix[col]=-1
	fprintf(colorful, "%d %d %d ",255 ,255 ,255 );
	break;
      }
    }
    fprintf(black, "\n");
    fprintf(colorful, "\n");
    mtx_unlock(mtx);
  }

	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	printf("It took %f seconds to write to files\n", elapsed);

  return 0;
}





int main_thrd( void *args ){
  //------------------------------
  //structure
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  float complex **z = thrd_info->z;
  int **attr = thrd_info->attr;
  int **iter = thrd_info->iter;
  const int ib = thrd_info->ib;
  const int istep = thrd_info->istep;
  const int sz = thrd_info->sz;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int *finish_flag = thrd_info->finish_flag;
  int degree = degree_global;

  //------------------------------
  // one row
  for ( int ix = ib; ix < sz; ix += istep ) {
    float complex *zix = z[ix];
    int *attr_ix = attr[ix];
    int *iter_ix = iter[ix];
    
    for ( int col = 0; col < sz; ++col ) {
      int conv;
      //------------------------------
      // newton's iteration
      for ( conv = 0, attr_ix[col] = -1 ; conv<127; ++conv ) { 
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

	//------------------------------
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
	  printf("unexpected degree\n");
	  exit(1);
	}
      }// newton iter loop end
      iter_ix[col]=conv;
    }// column loop end

    
    //------------------------------
    mtx_lock(mtx);
    iter[ix] = iter_ix;
    attr[ix] = attr_ix;
    finish_flag[ix] = 1;
    /* printf("Thread %d : ", tx); */
    /* for (int i = 0; i < sz; i++) { */
    /*   //      printf("%lf+j%lf  ",creal(input[i]),cimag(input[i])); */
    /*   printf("%lf+j%lf %d %d ", creal(zix[i]), cimag(zix[i]), attr_ix[i] , iter_ix[i]); */

    /*   /\* printf("attr : %d \n", attr_ix[i]); *\/ */
    /*   /\* printf("iter : %d \n", iter_ix[i]); *\/ */
    /* } */
    /* printf("\n"); */
    /* //status[tx].val = ix + istep; */
    mtx_unlock(mtx);
    cnd_signal(cnd);

    
    // In order to illustrate thrd_sleep and to force more synchronization
    // points, we sleep after each line for one micro seconds.
    //thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=1000}, NULL);
    
  }// row in thread

  return 0;
}





int main(int argc, char *argv[]){
  //------------------------------
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

  //------------------------------
  //read date from shell
  for (int i = 1; i < 4; i++) {
    if (sscanf(argv[i], "-t%d", &threads) == 1) {
      //printf("Extracted value: %d\n", threads);
    } else if (sscanf(argv[i], "-l%d", &size) == 1) {
      //printf("Extracted value: %d\n", size);
    }
    else if (sscanf(argv[i], " %d", &degree_global) == 1) {
      //printf("Extracted value: %d\n", degree_global);
    }
    else {
      printf("Invalid input.\n");
      return 1;
    }
  }
  printf("\nthreads : %d\nsize : %d\ndegree_global : %d\n\n", threads, size, degree_global);

  //------------------------------
  //file
  FILE *black = fopen("black.ppm", "w");
  FILE *colorful = fopen("colorful.ppm", "w");
  if (black == NULL || black == NULL) {
    perror("Error opening the file");
    return 1;
  }
  fprintf(black, "P3\n");
  fprintf(black, "%d %d \n", size, size);
  fprintf(black, "255\n");
  fprintf(colorful, "P3\n");
  fprintf(colorful, "%d %d \n", size, size);
  fprintf(colorful, "255\n");

  
  //------------------------------
  const int sz = size;
  int *finish_flag = (int *)malloc(sz*sizeof(int));
  float complex **z = (float complex**) malloc(sz*sizeof(float complex *));
  int **attr = (int**) malloc(sz*sizeof(int*));
  int **iter = (int**) malloc(sz*sizeof(int*));
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
  for ( int ix = 0; ix < sz; ++ix ){
    finish_flag[ix]=0;
  }

  //------------------------------
  //initialize
  double step = 4.0 / (sz - 1); 
  double real, imag;
  for (int i = 0; i < sz; i++) {
    for (int j = 0; j < sz; j++) {
      double real = -2.0 + i * step;
      double imag = 2.0 - j * step;
      z[i][j] = real+_Complex_I*imag;//CMPLX(real, imag);
      }
  }

  //------------------------------
  //thread initial
  const int nthrds = threads;
  thrd_t thrds[nthrds];
  thrd_info_t thrds_info[nthrds];
  
  thrd_t thrd_write;
  thrd_write_info_t thrd_write_info;

  //mutex initial and condition initial
  mtx_t mtx;
  mtx_init(&mtx, mtx_plain);
  cnd_t cnd;
  cnd_init(&cnd);

  //------------------------------
  //computation thread create
  for ( int tx = 0; tx < nthrds; ++tx ) {
    thrds_info[tx].z = (float complex**) z;
    thrds_info[tx].attr = attr;
    thrds_info[tx].iter = iter;

    //row switching para
    thrds_info[tx].ib = tx;
    thrds_info[tx].istep = nthrds;
    thrds_info[tx].sz = sz;
    thrds_info[tx].tx = tx;

    //mutex and condition
    thrds_info[tx].mtx = &mtx;
    thrds_info[tx].cnd = &cnd;

    //
    thrds_info[tx].finish_flag = finish_flag;
    
    //create
    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));
    if ( r != thrd_success ) {
      printf("failed to create thread\n");
      exit(1);
    }
    //thrd_detach(thrds[tx]);
  }

  //------------------------------
  //write thread create
  thrd_write_info.iter = iter;
  thrd_write_info.attr = attr;
  thrd_write_info.sz   = sz;
  thrd_write_info.mtx  = &mtx;
  thrd_write_info.cnd  = &cnd;
  thrd_write_info.colorful = colorful;
  thrd_write_info.black    = black;
  thrd_write_info.finish_flag = finish_flag;

  int x = thrd_create(&thrd_write, write_thrd, (void*) (&thrd_write_info));
  if ( x != thrd_success ) {
    printf("failed to create thread\n");
    exit(1);
  }


  //------------------------------
  //thread join 
  for(int t;t<nthrds;t++){
    thrd_join(thrds[t], NULL);
  }
  thrd_join(thrd_write, NULL);


  //------------------------------
  //free (but i dont know the atter entry and iter entry, maybe we can use them in thrad)
  free(zentries);
  free(attr_entries);
  free(iter_entries);
  free(z);
  free(iter);
  free(attr);

  mtx_destroy(&mtx);
  cnd_destroy(&cnd);

  fclose(black);
  fclose(colorful);

  return 0;
}
 

