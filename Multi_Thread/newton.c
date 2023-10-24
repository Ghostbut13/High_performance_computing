#include <stdlib.h>
#include <stdio.h>
#include <thread.h>
#include <complex.h>

double newton_iter(double a);//functions for calculation 

typedef struct {
	const double complex** z0;
	const double complex** root;
	mtx_t* mtx_root;
	cnd_t* cnd_root;
	int row_i_begin;
	int row_i_end;
	int row_i_step;
	int** iteration;
	int col_z;//this is size of array
} thrd_info_t;


typedef struct {
	double complex** z0;
	double complex** root;
	int** iteration;
	int col_z;
	mtx_t* mtx_root;
	cnd_t* cnd_root;
} thrd_dumpppm_info_t;

int main_computing_thrd(void* argc) {
	const thrd_info_t* thrd_info = (thrd_info_t*)argc;
	double complex** z0 = thrd_info->z0;
	double complex** root = thrd_info->root;
	mtx_t* mtx_root	= thrd_info->mtx_root;
	cnd_t* cnd_root	= thrd_info->cnd_root;
	int row_i_begin	= thrd_info->row_i_begin;
	int row_i_end	= thrd_info->row_i_end;
	int row_i_step	= thrd_info->row_i_step;
	int col_z		= thrd_info->col_z;
	int** iteration  = thrd_info->iteration;

	for (int ix = row_i_begin; ix < row_i_end; ix += row_i_step) {

		//calculate each element
		double* root_one_row = (double*)malloc(col_z * sizeof(double));
		double* iteration_one_row = (double*)malloc(col_z * sizeof(double));
		for (int jx = 0; jx < col_z; jx++) {
			root_one_row[jx] = newton_iter(z0[ix][jx]);
			iteration_one_row[jx] = jx++;
		}

		//---------------------------------------------




		//---------------------------------------------

		//mutex
		mtx_lock(mtx_root);
		root[ix] = root_one_row;
		
		mtx_lock(mtx_root);

		//condition variables
		cnd_signal(cnd_root);
	}

	return 0;
}


int main_dumpppm_thrd(void* argc) {
	const thrd_dumpppm_info_t* thrd_dumpppm_info = (thrd_dumpppm_info_t*)argc;
	double complex** z0 = thrd_dumpppm_info->z0;
	double complex** root = thrd_dumpppm_info->root;
	int** iteration = thrd_dumpppm_info->iteration;
	int col_z = thrd_dumpppm_info->col_z;
	mtx_t* mtx_root = thrd_info->mtx_root;
	cnd_t* cnd_root = thrd_info->cnd_root;


	return 0;
}
//------------------------------------------------------------//

int main(int argc, char* agrv[]) {
	// variables
	const int sz;
	// alloc array
	const double complex** z0 = (double complex**)malloc(sizeof(double complex*) * sz);// maybe this format not for complex array, we need change it
	double complex* z0_entries = (double complex*)malloc(sizeof(double complex) * sz * sz);// maybe this format not for complex array
	for (int ix = 0, jx = 0; ix < sz; jx += size, ix++) {
		z0[ix] = z0_entries + jx;
	}
	const double complex** root = (double complex**)malloc(sizeof(double complex*) * sz);
	
	// initialize z0
	double re_interval = 4.0 / sz;
	double im_interval = 4.0 / sz;
	double complex z = -2.0 + 2.0 * _Complex_I;
	for (int ix = 0; ix < sz; ix++) {
		for (int jx = 0; jx < sz; jx++) {
			z += re_interval;
		}
		z=-2.0+(2.0-im_interval)*_Complex_I
	}

	// thrd number and create
	const int nthrd = 8;
	thrd_t thrds[nthrd];
	thrd_t thrd_dumpppm;

	// information struct create
	thrd_info_t thrds_info[nthrd];
	thrd_dumpppm_info_t thrd_dumpppm_info;

	// mutex and condition create and init
	mtx_t mtx_root;
	cnd_t cnd_root;
	mtx_init(&mtx_root, mtx_plain);
	cnd_init(&cnd_root);

	// create main_computing_thrd
	for (int thrd_number = 0; thrd_number < nthrds; ++thrd_number) {
		//structure member 
		thrds_info[thrd_number].z0 = z0;
		thrds_info[thrd_number].root = root;
		thrds_info[thrd_number].row_i_begin = thrd_number;
		thrds_info[thrd_number].row_i_end = sz;
		thrds_info[thrd_number].row_i_step = nthrds;
		thrds_info[thrd_number].col_z = sz;
		thrds_info[thrd_number].mtx_root = &mtx_root;
		thrds_info[thrd_number].cnd_root = &cnd_root;

		
		//create thrd
		int r = thrd_create(thrds + thrd_number, main_computing_thrd, (void*)(thrds_info + thrd_number));
		if (r != thrd_success) {
			fprintf(stderr, "failed to create thread\n");
			exit(1);
		}
		thrd_detach(thrds[thrd_number]);
	}


	// create dump_ppm_thrd


	// threads joining
	// 
	

	//free
	free(z0_entries);
	free(z0);
	free(root);



	// distory mutex and condition
	mtx_destroy(&mtx_root);
	cnd_destroy(&cnd_root);

	return 0;
}