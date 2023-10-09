#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define N_BYTE_OF_DOUBLE 8
#define N_BYTE_OF_HEAP 5
//#define VS2022_DEBUG 1
#define HOW_MANY_DISTANCE 3464
//#define N_DATA_INPUT 10000

#ifdef VS2022_DEBUG
/* #define PATH_INPUT "E:/cpp/github_repository/High_performance_computing/OpenMP/sim/input_1e5.txt" */
/* #define PATH_OUTPUT "E:/cpp/github_repository/High_performance_computing/OpenMP/sim/output.txt" */
#define PATH_INPUT "/home/hpcuser053/work/High_performance_computing/OpenMP/sim/input_1e4.txt"
#define PATH_OUTPUT "/home/hpcuser053/work/High_performance_computing/OpenMP/sim/output.txt"

#else
//#define PATH_INPUT "/home/hpcuser053/test_data/cells_1e4"
//#define PATH_INPUT "/home/hpcuser053/work/High_performance_computing/OpenMP/sim/input.txt"
//#define PATH_INPUT "cells"
#define PATH_INPUT "/home/hpcuser053/work/High_performance_computing/OpenMP/extracted/cell_distance/cells"

#define PATH_OUTPUT "/home/hpcuser053/work/High_performance_computing/OpenMP/sim/output.txt"

#endif // VS2022_DEBUG


int main(int argc, char * argv[]) {

  //printf("%d\n",5);
  FILE* FP = fopen("cells", "r");
  if (FP == NULL) {
    printf("No Such File !! ");
    return -1;
  }
  
  int N_DATA_INPUT=0,c;
  while((c = fgetc(FP))!=EOF){
    if(c=='\n')
      N_DATA_INPUT++;
  }
  
  
  double start[3], end[3], timee[3];
  double* ThreeD = (double*)malloc(N_DATA_INPUT * 3 * sizeof(double));
  double** as = (double**)malloc(N_DATA_INPUT * sizeof(double*));
  double* OUT = (double*)malloc(HOW_MANY_DISTANCE * sizeof(double));
  int* cnt_d = (int*)malloc(HOW_MANY_DISTANCE * sizeof(int));


  //heap manage
  //----------------------------------
  start[0] = omp_get_wtime();
  for (int ix = 0, jx = 0; ix < N_DATA_INPUT; ix = ix + 1, jx = jx + 3) {
    as[ix] = ThreeD + jx;
  }
  for (size_t ix = 0; ix < N_DATA_INPUT; ++ix)
    for (size_t jx = 0; jx < 3; ++jx)
      as[ix][jx] = 0;
  for (int ix = 0; ix < HOW_MANY_DISTANCE; ix++) {
    OUT[ix] = ((double)ix / 100);
    cnt_d[ix] = 0;
  }

  end[0] = omp_get_wtime();
  timee[0] = end[0] - start[0];
  //----------------------------------


  
  //parse file
  //----------------------------------
  start[1] = omp_get_wtime();
  fseek(FP,0,SEEK_SET);
  //#pragma omp parallel
  for (int ix = 0; ix < N_DATA_INPUT * 3; ix += 1) {
    fscanf(FP, "%lf", &ThreeD[ix]);
    /*if (ix<100) {
      printf("%2.3f\n", ThreeD[ix]);
      }*/
  }
  end[1] = omp_get_wtime();
  timee[1] = end[1] - start[1];
  //----------------------------------


  
  //calculate
  //----------------------------------
  int threads_num;
  if(sscanf(argv[1],"%*[^0-9]%d",&threads_num)==1){
    //printf("aaaaaaa\n");
  }
  start[2] = omp_get_wtime();
  omp_set_num_threads(threads_num);
  int distance;
  //#ifndef VS2022_DEBUG
#pragma omp parallel for reduction(+:cnt_d[:3464]) //collapse(2) //
  //#endif //VS2022_DEBUG
  for (int ix=0; ix < N_DATA_INPUT - 1; ix ++) {
    for (int jx = ix + 1; jx < N_DATA_INPUT; jx++) {
      distance = (int)(sqrt(
			     (as[ix][1] - as[jx][1]) * (as[ix][1] - as[jx][1]) +
			     (as[ix][2] - as[jx][2]) * (as[ix][2] - as[jx][2]) +
			     (as[ix][0] - as[jx][0]) * (as[ix][0] - as[jx][0])
			      ) * 100);
      //#pragma omp critical
      ++cnt_d[distance];
    }
  }
  end[2] = omp_get_wtime(); 
  timee[2] = end[2] - start[2];

  
/*   //----------------------------------  //---------------------------------- */
/* #pragma omp parallel  */
/*   { */
/*     int id,nthread; */
/*     int distance; */
/*     id=omp_get_thread_num(); */
/*     nthread=omp_get_num_threads(); */
/*     #pragma omp for */
/*     for (int ix=id; ix < N_DATA_INPUT - 1; ix += nthread) { */
/*       for (int jx = ix + 1; jx < N_DATA_INPUT; jx++) { */
/* 	distance = (int)((sqrt( */
/* 			       (as[ix][1] - as[jx][1]) * (as[ix][1] - as[jx][1]) + */
/* 			       (as[ix][2] - as[jx][2]) * (as[ix][2] - as[jx][2]) + */
/* 			       (as[ix][0] - as[jx][0]) * (as[ix][0] - as[jx][0]) */
/* 			       )) * 100); */
/* 	//#pragma omp atomic */
/* 	cnt_d[distance]++; */
/*       } */
/*     } */
/*     printf("%d       %d\n",id,nthread); */
/*   } */
/*   end[2] = omp_get_wtime();  */
/*   timee[2] = end[2] - start[2]; */
/*   //----------------------------------  //---------------------------------- */
  
  FILE* out_FP = fopen(PATH_OUTPUT, "w");
  if (out_FP != NULL) {
    for (int ixx = 0; ixx < HOW_MANY_DISTANCE; ixx++){
      //if(cnt_d[ixx]!=0)
	fprintf(out_FP, "%05.2lf %d\n", OUT[ixx], cnt_d[ixx]);
    }
  }

 
  for (int ixx = 0; ixx < HOW_MANY_DISTANCE; ixx++) {
    //if( OUT[ixx]<10 )
      printf("%05.2f %d\n", OUT[ixx], cnt_d[ixx]);
    //else
    //printf("%.2lf %d\n", OUT[ixx], cnt_d[ixx]);
  }

  //----------------------------------  //----------------------------------

  
  //checking and dirtywork
  //printf("%s\n", PATH_INPUT);
  //printf("heap manage: %lf     ,parse file: %lf,    calculate : %lf      all: %lf\n", timee[0], timee[1], timee[2], timee[0] + timee[1] + timee[2]);
  free(as);
  free(ThreeD);
  fclose(FP);
  //fclose(out_FP);
  return 0;
}
