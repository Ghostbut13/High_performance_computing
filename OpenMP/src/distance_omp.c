#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define N_BYTE_OF_DOUBLE 8
#define N_BYTE_OF_HEAP 5
//#define VS2022_DEBUG 1
#define HOW_MANY_DISTANCE 3464
#define N_DATA_INPUT 100000

#ifdef VS2022_DEBUG
#define PATH_INPUT "E:/cpp/github_repository/High_performance_computing/OpenMP/sim/input_1e5.txt"
#define PATH_OUTPUT "E:/cpp/github_repository/High_performance_computing/OpenMP/sim/output.txt"

#else
#define PATH_INPUT "/home/hpcuser053/test_data/cells_1e5"
#define PATH_OUTPUT "/home/hpcuser053/work/High_performance_computing/OpenMP/sim/output.txt"

#endif // VS2022_DEBUG


int main() {

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


  
  FILE* FP = fopen(PATH_INPUT, "r+");
  if (FP == NULL) {
    printf("No Such File !! ");
    system("pause");
    free(as);
    free(ThreeD);
    return 0;
  }
  
  //parse file
  //----------------------------------
  start[1] = omp_get_wtime();
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
  start[2] = omp_get_wtime();
  omp_set_num_threads(100);
  int distance;
#ifndef VS2022_DEBUG
#pragma omp parallel for collapse(2) reduction(+:cnt_d[:3464])
#endif //VS2022_DEBUG
  for (int ix=0; ix < N_DATA_INPUT - 1; ix ++) {
    for (int jx = ix + 1; jx < N_DATA_INPUT; jx++) {
      distance = (int)(sqrt(
			     (as[ix][1] - as[jx][1]) * (as[ix][1] - as[jx][1]) +
			     (as[ix][2] - as[jx][2]) * (as[ix][2] - as[jx][2]) +
			     (as[ix][0] - as[jx][0]) * (as[ix][0] - as[jx][0])
			      ) * 100);
      //#pragma omp atomic
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
      for (int ixx = 0; ixx < HOW_MANY_DISTANCE; ixx++)
          fprintf(out_FP, "disdance : %lf,     cnt: %d\n", OUT[ixx], cnt_d[ixx]);
  }


 for (int ixx = 0; ixx < HOW_MANY_DISTANCE; ixx++) {
    printf("disdance : %lf,     cnt: %d\n", OUT[ixx], cnt_d[ixx]);
  }

  //----------------------------------  //----------------------------------

  
  //checking and dirtywork
  printf("%s\n", PATH_INPUT);
  printf("heap manage: %lf     ,parse file: %lf,    calculate : %lf      all: %lf\n", timee[0], timee[1], timee[2], timee[0] + timee[1] + timee[2]);
  free(as);
  free(ThreeD);
  fclose(FP);
  fclose(out_FP);
#ifdef VS2022_DEBUG
  system("pause");
#endif
  return 0;
}
