#include <stdlib.h>
#include <stdio.h>
#include <thread.h>
#include <complex.h>

double func_re(double a);//functions for calculation 
double func_im(double a);//maybe we only need one for complex number, not split them

typedef struct {
	int a;
} status_t; //maybe we need this to document key information like iteration number

typedef struct {
	const double** x0_im;
	const double** x0_re;//maybe we dont need split re and im
	const double** root_im;
	const double** root_re;
	mtx_t* mtx_root;
	cnd_t* cnd_root;
	int row_i_begin;
	int row_i_end;
	int row_i_step;
	int col_z;//this is size of array
	status_t* status_root;
} thrd_info_t;

typedef struct {
	int a;
} thrd_info_check_t;	//add more

typedef struct {
	int a;
} thrd_info_dumpppm_t;	//add more

int main_computing_thrd(void* argc) {
	const thrd_info_t* thrd_info = (thrd_info_t*)argc;
	double** x0_im	= thrd_info->x0_im;
	double** x0_re	= thrd_info->x0_re;   //maybe we dont need split re and im
	double** root_im = thrd_info->root_im;
	double** root_re = thrd_info->root_re;
	mtx_t* mtx_root	= thrd_info->mtx_root;
	cnd_t* cnd_root	= thrd_info->cnd_root;
	int row_i_begin	= thrd_info->row_i_begin;
	int row_i_end	= thrd_info->row_i_end;
	int row_i_step	= thrd_info->row_i_step;
	int col_z		= thrd_info->col_z;
	status_t* status_root = thrd_info->status_root;


	for (int ix = row_i_begin; ix < row_i_end; ix += row_i_step) {

		//calculate each element
		double* root_im_one_row = (double*)malloc(col_z * sizeof(double));
		double* root_re_one_row = (double*)malloc(col_z * sizeof(double));
		for (int jx = 0; jx < col_z; jx++) {
			root_im_one_row[jx] = func_im(x0_im[ix][jx]);
			root_re_one_row[jx] = func_re(x0_re[ix][jx]);
		}

		//mutex
		mtx_lock(mtx_root);
		root_im[ix] = root_im_one_row;
		root_re[ix] = root_re_one_row;
		//maybe some information in mutex
		mtx_lock(mtx_root);

		//condition variables
		cnd_signal(cnd_root);
	}

	return 0;
}

int main_check_thrd(void* argc) {
	return 0;// we need this, we need cnd_wait here and also free()
}

int dump_ppm_thrd(void* argc) {
	return 0;//also we need this, maybe we also need cnd_wait here
}
//------------------------------------------------------------//

int main(int argc, char* agrv[]) {
	// variables
	const int sz;
	// alloc array
	const double complex** x0 = (double complex**)malloc(sizeof(double complex*) * sz);// maybe this format not for complex array, we need change it
	double complex* x0_entries = (double complex*)malloc(sizeof(double complex) * sz * sz);// maybe this format not for complex array
	for (int ix = 0, jx = 0; ix < sz; jx += size, ix++) {
		x0[ix] = x0_entries + jx;
	}
	const double complex** root = (double complex**)malloc(sizeof(double complex*) * sz);
	
	// initialize x0
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
	thrd_t thrd_check;

	// information struct create
	thrd_info_t thrds_info[nthrd];
	thrd_info_check_t thrd_info_check;//now we dont have them
	//thrd_info_dumpppm_t thrd_info_dumpppm;
	//status_t status;   //now we dont have them

	// mutex and condition create and init
	mtx_t mtx_root;
	cnd_t cnd_root;
	mtx_init(&mtx_root, mtx_plain);
	cnd_init(&cnd_root);

	// create main_computing_thrd
	for (int thrd_number = 0; thrd_number < nthrds; ++thrd_number) {
		//structure member 
		thrds_info[thrd_number].x0_im = cimag(x0);
		thrds_info[thrd_number].x0_re = creal(x0);//maybe the creal is not ok here, also maybe pass whole complex number is better than individually im and re
		thrds_info[thrd_number].root_im = cimag(root);
		thrds_info[thrd_number].root_re = creal(root);//maybe the creal is not ok here, also maybe pass whole complex number is better than individually im and re
		thrds_info[thrd_number].row_i_begin = thrd_number;
		thrds_info[thrd_number].row_i_end = sz;
		thrds_info[thrd_number].row_i_step = nthrds;
		thrds_info[thrd_number].col_z = sz;
		thrds_info[thrd_number].mtx_root = &mtx_root;
		thrds_info[thrd_number].cnd_root = &cnd_root;
		//thrds_info[thrd_number].status_root = status; // now we dont have status
		
		//create thrd
		int r = thrd_create(thrds + thrd_number, main_computing_thrd, (void*)(thrds_info + thrd_number));
		if (r != thrd_success) {
			fprintf(stderr, "failed to create thread\n");
			exit(1);
		}
		thrd_detach(thrds[thrd_number]);
	}

	// create main_computing_thrd

	// create dump_ppm_thrd

	//free
	free(x0_entries);
	free(x0);
	free(root);


	// distory mutex and condition
	mtx_destroy(&mtx_root);
	cnd_destroy(&cnd_root);

	return 0;
}