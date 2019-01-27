#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
extern "C" {
#include "direct.h"
#include "../directtest/gkls.h"
#include "../directtest/rnd_gen.h"
}

double tst_obj(int n, const double *xy, int *undefined_flag, void *unused)
{
	double f, *x;
	x = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		x[i] = xy[i];
	}
	f = GKLS_D_func(x);
	// f = GKLS_ND_func(x); // -- for ND-type test function 
	// f = GKLS_D2_func(x); // -- for D2-type test function 
	free(x);
	return f;
}

double *a, *b, *x; //sets the search area for task
std::chrono::steady_clock sc; //for runtime measurement

void print_error_msg(int);

int main()
{
	int error_code;    /* error codes variable */
	int func_num;      /* test function number within a class     */
	std::ofstream fp;

	std::cout << "GKLS-Generator of Classes of ND, D, and D2 Test Functions";
	std::cout << std::endl << "for Global Optimization,";
	std::cout << std::endl << "(C) 2002-2003, M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev" << std::endl << std::endl;

	GKLS_dim = 2;
	GKLS_num_minima = 10;
	if ((error_code = GKLS_domain_alloc()) != GKLS_OK)
		return error_code;
	GKLS_global_dist = 2.0 / 3.0;
	GKLS_global_radius = 0.5*GKLS_global_dist;
	GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
	if ((error_code = GKLS_parameters_check()) != GKLS_OK)
		return error_code;

	int n = GKLS_dim;

	fp.open("results.txt", std::ios::out);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}

	double eps;

	std::cout << "Set accuracy" << std::endl << "It can significantly affect on performance!" << std::endl;
	std::cin >> eps;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
	}

	a = new double[n];
	b = new double[n];
	x = new double[n];

	double minf; //result of direct method
				 /* Generate the class of 100 D-type functions */
	auto astart = sc.now(); //start time for full set of tests

	for (func_num = 1; func_num <= 100; func_num++)
	{
		if ((error_code = GKLS_arg_generate(func_num)) != GKLS_OK) {
			print_error_msg(error_code);
			return error_code;
		}

		for (int i = 0; i < n; i++) {
			a[i] = GKLS_domain_left[i];
			b[i] = GKLS_domain_right[i];
		} //set search area for test

		std::cout << std::endl << "Generating the function number " << func_num << std::endl;
		fp << "D-type function number " << func_num;
		fp << std::endl << "of the class with the following parameters:";
		fp << std::endl << "    global minimum value   = " << GKLS_global_value << ';';



		int maxfeval, maxits;
		int info;

		int force_stop = 0;

		maxfeval = 50000;
		maxits = 500;

		/* Information about global minimizers */
		if (GKLS_glob.gm_index == 0)
			fp << std::endl << "An error during the global minimum searching was occurred!";
		else {
			if (GKLS_glob.num_global_minima == 1) {
				fp << std::endl << std::endl << "There is one global minimizer.";
				fp << std::endl << "The global minimizer coordinates are: ";
			}
			else {
				fp << std::endl << std::endl << "There are " << GKLS_glob.num_global_minima << " global minimizers.";
				fp << std::endl << "The global minimizers coordinates are: ";
			}
			for (unsigned int i = 0; i < GKLS_glob.num_global_minima; i++) {
				fp << i + 1 << ":	";
				double *x = GKLS_minima.local_min[GKLS_glob.gm_index[i]];
				fp << '[';
				for (int k = 0; k < n; k++) {
					if (k != n - 1)
						fp << std::setprecision(3) << x[k] << ' ';
					else
						fp << std::setprecision(3) << x[k] << ']' << std::endl;
				}
			}
		}

		double glob = DIRECT_UNKNOWN_FGLOBAL;
		glob = GKLS_global_value;

		/* Function evaluating */
		auto start = sc.now(); //start time for one function evaluating

		info = direct_optimize(tst_obj, NULL, n, a, b, x, &minf,
			&maxfeval, maxits,
			0, 0,
			10e-4, 10e-4,
			0, -1.0,
			&force_stop,
			glob, eps,
			stdout, DIRECT_ORIGINAL);

		auto end = sc.now();
		auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

		fp << "Search complete:" << std::endl << "Obtained min: " << std::setprecision(4) << minf << std::endl << \
			"Real global minimum: " << std::setprecision(4) << glob << std::endl << \
			"Diff: " << std::setprecision(4) << 1.0*fabs(glob - minf) << std::endl;
		fp << "Function evalautions: " << maxfeval << std::endl;
		fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
		fp << std::endl << std::endl;


		/* Deallocate memory */
		GKLS_free();

	}

	std::cout << std::endl;

	auto aend = sc.now();
	auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

	fp << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << 's' << std::endl;

	/* Close files */
	fp.close();


	/* Deallocate the boundary vectors */
	delete[]a; delete[]b; delete[]x;
	GKLS_domain_free();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
	return 0;
}

/* Print an error message */
void print_error_msg(int error_code)
{
	switch (error_code)
	{
	case GKLS_OK:
		std::cerr << "\nGKLS_OK: There is no error.";
		break;
	case GKLS_DIM_ERROR:
		std::cerr << "\nGKLS_DIM_ERROR: The problem dimension is out of the valid range [1," << (unsigned int)NUM_RND << "].";
		break;
	case GKLS_NUM_MINIMA_ERROR:
		std::cerr << "\nGKLS_NUM_MINIMA_ERROR: The number of local minima must be greater than 1.";
		break;
	case GKLS_FUNC_NUMBER_ERROR:
		std::cerr << "\nGKLS_FUNC_NUMBER_ERROR: The number of the test function to be generated is out of the range [1,100].";
		break;
	case GKLS_BOUNDARY_ERROR:
		std::cerr << "\nGKLS_BOUNDARY_ERROR: The admissible region boundary vectors are not defined or ill-defined.";
		break;
	case GKLS_GLOBAL_MIN_VALUE_ERROR:
		std::cerr << "\nGKLS_GLOBAL_MIN_VALUE_ERROR: The global minimum value must be greater than " << (double)GKLS_PARABOLOID_MIN;
		break;
	case GKLS_GLOBAL_DIST_ERROR:
		std::cerr << "\nGKLS_GLOBAL_DIST_ERROR: The distance from the paraboloid vertex to the global minimizer is too great.";
		break;
	case GKLS_GLOBAL_RADIUS_ERROR:
		std::cerr << "\nGKLS_GLOBAL_RADIUS_ERROR: The radius of the attraction region of the global minimizer is too high.";
		break;
	case GKLS_MEMORY_ERROR:
		std::cerr << "\nGKLS_MEMORY_ERROR: There is not enough memory to allocate.";
		break;
	case GKLS_DERIV_EVAL_ERROR:
		std::cerr << "\nGKLS_DERIV_EVAL_ERROR: An error occurs during derivative evaluation.";
		break;
	case GKLS_FLOATING_POINT_ERROR:
	default:
		std::cerr << "\nUnknown error.";
	}
	std::cerr << std::endl;
} /* print_error_msg() */