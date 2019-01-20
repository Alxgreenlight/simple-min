#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "../solver/solver.h"

extern "C" {
#include "gkls.h"
#include "rnd_gen.h"
}


double *a, *b; //sets the search area for task

void print_error_msg(int);

std::chrono::steady_clock sc; //for runtime measurement


int main()
{
	int error_code;    /* error codes variable */
	int func_num;      /* test function number within a class     */
	std::ofstream fp;

	std::cout << "GKLS-Generator of Classes of ND, D, and D2 Test Functions";
	std::cout << std::endl << "for Global Optimization,";
	std::cout << std::endl << "(C) 2002-2003, M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev" << std::endl << std::endl;


	/* Set the input parameters */
	/*if ((error_code=GKLS_set_default()) != GKLS_OK) {
		print_error_msg(error_code);
		return error_code;
	}*/
	GKLS_dim = 2;
	GKLS_num_minima = 10;
	if ((error_code = GKLS_domain_alloc()) != GKLS_OK)
		return error_code;
	GKLS_global_dist = 2.0 / 3.0;
	GKLS_global_radius = 0.5*GKLS_global_dist;
	GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
	if ((error_code = GKLS_parameters_check()) != GKLS_OK)
		return error_code;

	solver::dim = GKLS_dim;

	fp.open("results.txt", std::ios::out);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}

	std::cout << "Set number of nodes per dimension" << std::endl << \
		"It can significantly affect on performance!" << std::endl;
	std::cin >> solver::nodes;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
		std::cin >> solver::nodes;
	}

	std::cout << "Set accuracy" << std::endl << "It can significantly affect on performance!" << std::endl;
	std::cin >> solver::eps;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
	}

	a = new double[solver::dim];
	b = new double[solver::dim];

	/* Generate the class of 100 D-type functions */
	auto astart = sc.now(); //start time for full set of tests

	for (func_num = 1; func_num <= 100; func_num++)
	{
		if ((error_code = GKLS_arg_generate(func_num)) != GKLS_OK) {
			print_error_msg(error_code);
			return error_code;
		}

		for (int i = 0; i < solver::dim; i++) {
			a[i] = GKLS_domain_left[i];
			b[i] = GKLS_domain_right[i];
		} //set search area for test

		solver::compute = GKLS_D_func;
		// solver::compute = GKLS_ND_func; // -- for ND-type test function 
		// solver::compute = GKLS_D2_func; // -- for D2-type test function 

		printf("\nGenerating the function number %d\n", func_num);

		fp << "D-type function number " << func_num;
		fp << std::endl << "of the class with the following parameters:";
		fp << std::endl << "    global minimum value   = " << GKLS_global_value << ';';

		solver::glob = GKLS_global_value;


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
				for (int k = 0; k < solver::dim; k++) {
					if (k != solver::dim - 1)
						fp << std::setprecision(3) << x[k] << ' ';
					else
						fp << std::setprecision(3) << x[k] << ']' << std::endl;
				}
			}
		}

		/* Function evaluating */
		auto start = sc.now(); //start time for one function evaluating

		int r = solver::min_search(a, b);

		if (r) {
			solver::checkErrors(r);
			fp.close();
			delete[]a; delete[]b;
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cin.get();
			GKLS_free();
			GKLS_domain_free();
			return -2;
		}

		auto end = sc.now();
		auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4) << solver::UPB << std::endl << \
			"Real global minimum: " << std::setprecision(4) << solver::glob << std::endl << "Lower bound: " << std::setprecision(4) << \
			solver::LOB << std::endl << "Diff: " << std::setprecision(4) << 1.0*fabs(solver::glob - solver::UPB) << std::endl;
		fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
		fp << std::endl << std::endl;


		/* Deallocate memory */
		GKLS_free();
	} /* for func_num*/

	std::cout << std::endl;

	auto aend = sc.now();
	auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

	fp << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << 's' << std::endl;

	/* Close files */
	fp.close();


	/* Deallocate the boundary vectors */
	delete[]a; delete[]b;
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