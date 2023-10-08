#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define N_BYTE_OF_DOUBLE 8
#define N_BYTE_OF_HEAP 5
#define DEBUG 1
#define HOW_MANY_DISTANCE 3464

#ifdef DEBUG
#define PATH_INPUT "E:/cpp/openmp/sim_cell_distance/input_1e5.txt"
#define N_DATA_INPUT 100000
#endif // DEBUG


int main() {
#pragma warning(disable : 6031)

	double start[3], end[3], timee[3];
	int distance;
	int ix, jx;
	double* OUT = (double*)malloc(HOW_MANY_DISTANCE * sizeof(double));
	int* cnt_d = (int*)malloc(HOW_MANY_DISTANCE * sizeof(int));
	for (int ix = 0; ix < HOW_MANY_DISTANCE; ix++) {
		OUT[ix] = ((double)ix / 100);
		cnt_d[ix] = 0;
	}
	double* ThreeD = (double*)malloc(N_DATA_INPUT * 3 * sizeof(double));
	double** as = (double**)malloc(N_DATA_INPUT * sizeof(double*));

	//heap manage
	start[0] = omp_get_wtime();
	for (int ix = 0, jx = 0; ix < N_DATA_INPUT; ix = ix + 1, jx = jx + 3) {
		as[ix] = ThreeD + jx;
	}
	for (size_t ix = 0; ix < N_DATA_INPUT; ++ix)
		for (size_t jx = 0; jx < 3; ++jx)
			as[ix][jx] = 0;
	end[0] = omp_get_wtime();
	timee[0] = end[0] - start[0];


	FILE* FP = fopen(PATH_INPUT, "r+");
	if (FP == NULL) {
		printf("No Such File !! ");
		system("pause");
		free(as);
		free(ThreeD);
		return 0;
	}
	//parse file
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


	//calculate
	start[2] = omp_get_wtime();
	//omp_set_num_threads(5);
	const int num_threads = omp_get_max_threads();
	omp_set_num_threads(num_threads);
//#pragma omp parallel for schedule(dynamic)
//	for (ix=0; ix < N_DATA_INPUT - 1; ix ++) {
//		for (jx = ix + 1; jx < N_DATA_INPUT; jx++) {
			//distance = (int)((sqrt( 
			//	(as[ix][1] - as[jx][1]) * (as[ix][1] - as[jx][1]) +
			//	(as[ix][2] - as[jx][2]) * (as[ix][2] - as[jx][2]) +
			//	(as[ix][0] - as[jx][0]) * (as[ix][0] - as[jx][0])
			//	)) * 100);
////#pragma omp atomic
//			cnt_d[distance]++;
//		}
//	}


#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (ix = 0; ix < N_DATA_INPUT - 1; ix++) {
			for (jx = ix + 1; jx < N_DATA_INPUT; jx++) {
				distance = (int)(sqrt((as[ix][1] - as[jx][1]) * (as[ix][1] - as[jx][1]) +
					(as[ix][2] - as[jx][2]) * (as[ix][2] - as[jx][2]) +
					(as[ix][0] - as[jx][0]) * (as[ix][0] - as[jx][0])) * 100);

//#pragma omp atomic
				cnt_d[distance]++;
			}
		}
	}
	end[2] = omp_get_wtime();
	timee[2] = end[2] - start[2];
	

	for (int ixx = 0; ixx < HOW_MANY_DISTANCE; ixx++) {
		printf("disdance : %lf,     cnt: %d\n", OUT[ixx], cnt_d[ixx]);
	}


	//checking and dirtywork
	printf("%s\n", PATH_INPUT);
	printf("heap manage: %lf     ,parse file: %lf,    calculate : %lf      all: %lf\n", timee[0], timee[1], timee[2], timee[0] + timee[1] + timee[2]);
	system("pause");
	free(as);
	free(ThreeD);
	fclose(FP);
	return 0;
}