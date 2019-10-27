#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <limits>
#include "../solver/R_Optim_parallel.hpp"
#include "../solver/UltraEstim.hpp"

/* Including GKLS libraries, writen with C language */

extern "C" {
#include "gkls.h"
#include "rnd_gen.h"
}


double *a, *b, *x; /* sets the search area for task
x is used for saving coordinates of solve*/
unsigned long long int fevals;
unsigned long int iters;

void print_error_msg(int);	/* print error in GKLS */

std::chrono::steady_clock sc; /* for runtime measurement */
std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();

double func(const double* x) {	/* wrapper for function providing to solver */
	double *xc;
	xc = new double[GKLS_dim];
	for (unsigned int i = 0; i < GKLS_dim; i++) {
		xc[i] = x[i];
	}
	double r = GKLS_D_func(xc);	/* returns a function value in a point with coordinates provided in x[] */
					/* GKLS_ND_func; -- for ND-type test function */
					/* GKLS_D2_func; -- for D2-type test function */
	delete[]xc;
	return r;
}

/* why so terrible wrapper with memory allocation in every call?
method requires function with const double pointer
but GKLS works only with non-const pointers */
/* It can be fixed in GKLS code, but tests weren't show any acceleration, it's strange */


int main()
{
	int error_code;    /* error codes variable */
	int func_num;      /* test function number within a class     */
	double maxdiff = std::numeric_limits<double>::min();	/* maximum error among all solutions */
	std::ofstream fp;	/* output file stream */
	double eps;	/* required accuracy */
	double L;
	unsigned short errors = 0;

	/* Set parameters of GKLS */
	GKLS_dim = 4;
	GKLS_num_minima = 10;
	if ((error_code = GKLS_domain_alloc()) != GKLS_OK)
		return error_code;
	GKLS_global_dist = 2.0 / 3.0;
	GKLS_global_radius = 0.5*GKLS_global_dist;
	GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
	if ((error_code = GKLS_parameters_check()) != GKLS_OK)
		return error_code;


	/* Try to open output file */
	fp.open("results.txt", std::ios::out);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}


	/* Set parametrs of method */
	/* Interactive */

<<<<<<< HEAD
=======

>>>>>>> 044632ce294e4e5141810844821e9234eb989076
	std::cout << "Set accuracy" << std::endl << "It can significantly affect on performance!" << std::endl;
	std::cin >> eps;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
		std::cin >> eps;
	}

	/* Allocating required memory */

	a = new double[GKLS_dim];	/* For left bound of search region */
	b = new double[GKLS_dim];	/* For right bound of search region */
	x = new double[GKLS_dim];	/* For global minimum coordinates found */
	L_accurate<double> La;
	//La.stay_fixed(GKLS_dim, 10000);
	PrOptimizer_v1<double> rOpt;

	/* Generate the class of 100 D-type functions */

	for (func_num = 1; func_num <= 100; func_num++)
	{
		/* Initializing GKLS with set parameters */

		if ((error_code = GKLS_arg_generate(func_num)) != GKLS_OK) {
			print_error_msg(error_code);
			return error_code;
		}

		/* set search area for test */

		for (unsigned int i = 0; i < GKLS_dim; i++) {
			a[i] = GKLS_domain_left[i];
			b[i] = GKLS_domain_right[i];
		} 
		
		/* Output func info */

		std::cout << std::endl << "Generating the function number " << func_num << std::endl;

		fp << "D-type function number " << func_num;
		fp << std::endl << "of the class with the following parameters:";
		fp << std::endl << "    global minimum value   = " << GKLS_global_value << ';';


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
				double *xcoor = GKLS_minima.local_min[GKLS_glob.gm_index[i]];
				fp << '[';
				for (unsigned int k = 0; k < GKLS_dim; k++) {
					if (k != GKLS_dim - 1)
						fp << std::setprecision(3) << xcoor[k] << ' ';
					else
						fp << std::setprecision(3) << xcoor[k] << ']' << std::endl;
				}
			}
		}
		try {
<<<<<<< HEAD
			PrOptimizer_v1<double> opt;
			opt.init(GKLS_dim, eps);
=======
			rOpt.init(GKLS_dim, eps);
>>>>>>> 044632ce294e4e5141810844821e9234eb989076

			auto start = sc.now(); /* start time for one function evaluating */

			/* perform a search */

			double UPB = rOpt.search(GKLS_dim, x, a, b, func);

			/* Search completed */
			auto end = sc.now();
			atime_span += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			rOpt.getInfo(fevals, iters, L);
			/* Differece between obtained result and real global minimum */
			double diff = 1.0*fabs(GKLS_global_value - UPB);
			if (diff > eps) {
				errors++;
			}
		
		maxdiff = diff > maxdiff ? diff : maxdiff;

		/* Output all obtained results */

		fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4) << UPB << std::endl << \
			"At [";
		std::copy(x, x + GKLS_dim, std::ostream_iterator<double>(fp, " ")); 
		fp << "], with L = " << L << std::endl;
		
		//fp << "Real accurate estimated L: " << La.ultraoptimizer_f(a, b, x, UPB, func) << std::endl;
		fp << "Real global minimum: " << std::setprecision(4) << GKLS_global_value << std::endl << std::endl << "Diff: " << std::setprecision(4) << diff << std::endl;
		fp << "Function evaluations: " << fevals << " in " << iters << " iterations" << std::endl;
		fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
		fp << std::endl << std::endl;
		}
		catch (std::exception &e) {
			std::cout << "Error occured: " << e.what() << std::endl;
			GKLS_free();
			fp.close();
			delete[]a; delete[]b; delete[]x;
			GKLS_domain_free();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cin.get();
			return -1;
		}
		/* Deallocate memory */
		

		GKLS_free();

	} /* for func_num */
	La.unfix();

	//std::cout << std::endl;

	/* Benchmarking completed */

	/* Output total search time */
	std::cout << "Errors: " << errors << " of 100" << std::endl;
	fp << "Total evaluation time: " << atime_span.count() << "ms" << std::endl;
	fp << "Maximum difference in this set: " << maxdiff << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << "ms" << std::endl;

	/* Close files */
	fp.close();

	/* Deallocate the boundary arrays */
	delete[]a; delete[]b; delete[]x;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
	GKLS_domain_free();
	return 0;

}


/* Print an error message */
/* From GKLS tutorial */
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
