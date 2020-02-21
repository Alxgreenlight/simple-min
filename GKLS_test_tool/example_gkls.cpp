#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <limits>
#include <vector>
#include "solver/R_Optim_Pure_Parallel.hpp"
#include "util/helper.hpp"

extern "C" {
#include "GKLStest/gkls.h"
#include "GKLStest/rnd_gen.h"
}


double *a, *b, *x; //sets the search area for task

void print_error_msg(int);

std::chrono::steady_clock sc; //for runtime measurement

double func(const double* x) {	/* wrapper for function providing to solver */
	double *xc;
	xc = new double[GKLS_dim];
	for (unsigned int i = 0; i < GKLS_dim; i++) {
		xc[i] = x[i];
	}
	double r = GKLS_D_func(xc);
	/* GKLS_ND_func; -- for ND-type test function */
	/* GKLS_D2_func; -- for D2-type test function */
	delete[]xc;
	return r;
}


int main(int argc, char **argv)
{
	int error_code;    /* error codes variable */
	int func_num;      /* test function number within a class     */
	double maxdiff = std::numeric_limits<double>::min();	/* maximum error among all solutions */
	std::ofstream fp;	/* output file stream */
	double eps;	/* required accuracy */

	/* file with results of experiment */
	fp.open("getEX.csv", std::ios::app);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}

	/* Set accuracy in interactive way */
	/* Using default the best efficient nodes per dimension = 4 */
	if (argc < 2)
    {
        helper::help(benchlib);
        return 0;
    }

	eps = atof(argv[1]);
	if (eps < std::numeric_limits<double>::min()){
		std::cerr << "Accuracy defined incorrect, exit" << std::endl;
		return -1;
	}

	/* Create a solver object */
	PPrOptimizer<double> gs;
	fp << "Parallel search" << std::endl;


	/* full GKLS (100 funcs in each) tests with varying dimensionality */
	for (int dn = 2; dn < 6; dn++) {

		/* GKLS initialization */
		GKLS_dim = dn;
		GKLS_num_minima = 10;
		if ((error_code = GKLS_domain_alloc()) != GKLS_OK)
			return error_code;
		GKLS_global_dist = 2.0 / 3.0;
		GKLS_global_radius = 0.5*GKLS_global_dist;
		GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
		if ((error_code = GKLS_parameters_check()) != GKLS_OK)
			return error_code;


		int summary; //Amount of funcs with wrong answer
		std::vector<int> allnodes; //Numbers of funcs for which wrong answers obtained

		/* Allocating required memory */

		a = new double[GKLS_dim];	/* For left bound of search region */
		b = new double[GKLS_dim];	/* For right bound of search region */
		x = new double[GKLS_dim];	/* For global minimum coordinates found */

			/* Number of function evaluations and algorithm iteartions during search */
		unsigned long long int fevs;
		unsigned long int its;

		/* Generate the class of 100 D-type functions */
			/* provide settings to solver */
			gs.init(GKLS_dim, eps);
			fevs = 0;
			its = 0;
			summary = 0;
			allnodes.clear();
			unsigned long long int time = 0;

			/* run solver with all 100 test functions */
			for (func_num = 1; func_num <= 100; func_num++)
			{
				if ((error_code = GKLS_arg_generate(func_num)) != GKLS_OK) {
					print_error_msg(error_code);
					return error_code;
				}

				/* set search area for test */

				for (unsigned int i = 0; i < GKLS_dim; i++) {
					a[i] = GKLS_domain_left[i];
					b[i] = GKLS_domain_right[i];
				}

				printf("\nGenerating the function number %d\n", func_num);
				auto astart = sc.now(); //start time for full set of tests
				double UPB = gs.search(GKLS_dim, x, a, b, func);
				auto aend = sc.now();
				auto atime_span = std::chrono::duration_cast<std::chrono::milliseconds>(aend - astart);
				time += atime_span.count();
				/* Number of function evaluations and algorithm iteartions during search */
				unsigned long long int fevals;
				unsigned long int iters;
				double Lmax;

				gs.getInfo(fevals, iters, Lmax);

				fevs += fevals;
				its += iters;

				/* if answer differ from real global minimum more than eps, the answer is wrong */
				if ((1.0*fabs(GKLS_global_value - UPB)) > eps) {
					summary++;
					allnodes.push_back(func_num);
				}


				/* Deallocate memory */
				GKLS_free();
			} /* for func_num*/
			
			/* Write results to file */
			fp << "Dim = " << GKLS_dim << ';';
			for (int &d : allnodes) {
				fp << d << ' ';
			}
			fp << ";time = " << time;
			fp << ";avg.time = " << static_cast<double>(time) / 100.0;
			fp << ';' << "avg. num. of evaluations = " << static_cast<long double>(fevs) / 100.0;
			fp << ';' << "avg. num. of iterations" << static_cast<long double>(its) / 100.0 << std::endl;

		std::cout << std::endl;


		/* Close files */

		delete[]a; delete[]b; delete[]x;
		GKLS_domain_free();
	}
	fp.close();
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
