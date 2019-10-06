#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>                          
#include "../solver/R_Optim.hpp"

extern "C" {
#include "vagris.h"
}

double *a, *b, *x; //sets the search area for task
unsigned long long int fevals;
unsigned long int iters;

std::chrono::steady_clock sc; //for runtime measurement
std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();

double func(const double *xx) {
	return random_func(xx[0], xx[1]);	//random_func returns value of function in x1,x2 (from vagris.h)
}

int main()
{
	int dim = 2;
	int nf; //number of test function
	std::ofstream fp; //output file for results
	int nodes;
	double eps;
	double L;
	unsigned short errors = 0;

	fp.open("results.txt", std::ios::out);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}

	/* Interactively set parameters */
	std::cout << "Set number of nodes per dimension" << std::endl << \
		"It can significantly affect on performance!" << std::endl;
	std::cin >> nodes;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
		std::cin >> nodes;
	}

	std::cout << "Set accuracy" << std::endl << "It can significantly affect on performance!" << std::endl;
	std::cin >> eps;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
		std::cin >> eps;
	}

	/*set left and right bounds of serach regions */
	a = new double[dim];
	b = new double[dim];
	x = new double[dim];

	for (int i = 0; i < dim; i++) {
		a[i] = 0.0;
		b[i] = 1.0;
	} //set search area for all tests to [0..1] [0..1] (like in tutorial for this tests)

	try {
		rOptimizer<double> opt;
		opt.init(dim, nodes, eps);
		for (nf = 1; nf <= 100; nf++)
		{
			set_random(nf);

			fp << std::endl << "Function number " << nf;
			fp << std::endl << "Its global minimizer is (" << std::setprecision(4) << rand_minimums[(nf - 1) * 2] << ',' << std::setprecision(4) << rand_minimums[(nf - 1) * 2 + 1] << ')' << std::endl;
			double glob = random_func(rand_minimums[(nf - 1) * 2], rand_minimums[(nf - 1) * 2 + 1]);
			std::cout << std::endl << "Computing function " << nf << "..." << std::endl;

			auto start = sc.now(); //start time for one function evaluating

			double UPB = opt.search(dim, x, a, b, func);

			auto end = sc.now();
			opt.getInfo(fevals, iters, L);
			atime_span += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			double diff = 1.0*fabs(glob - UPB);
			if (diff > eps) {
				errors++;
			}
			fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4) << UPB << std::endl;
			fp << "At [";
			std::copy(x, x + dim, std::ostream_iterator<double>(fp, " ")); 
			fp << "], with L = " << L << std::endl;
			fp <<"Real global minimum: " << std::setprecision(4) << glob << std::endl \
				<< "Diff: " << std::setprecision(4) << 1.0*fabs(glob - UPB) << std::endl;
			fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
			fp << std::endl << std::endl;
		} /* for nf*/
		std::cout << std::endl;

		std::cout << "Errors: " << errors << " of 100" << std::endl;
		fp << "Total evaluation time: " << atime_span.count() << "ms" << std::endl;
		std::cout << "Total evaluation time: " << atime_span.count() << "ms" << std::endl;
	}
	catch (std::exception &e) {
		std::cout << "Error occured: " << e.what() << std::endl;
		fp.close();
		delete[]a; delete[]b; delete[]x;
	}
	fp.close();
	delete[]a; delete[]b; delete[]x;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
	return 0;
}






