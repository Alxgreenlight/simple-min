#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "../solver/gridsolver_hlp.hpp"

extern "C" {
#include "vagris.h"
}

double *a, *b, *x; //sets the search area for task

std::chrono::steady_clock sc; //for runtime measurement

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

	//GridSolver<double> gs;
	GridSolverOMP<double> gs;
	//GridSolverHLP<double> gs;

	gs.setparams(nodes, eps);

	auto astart = sc.now(); //start time for full set of tests
	for (nf = 1; nf <= 100; nf++)
	{
		set_random(nf);

		fp << std::endl << "Function number " << nf;
		fp << std::endl << "Its global minimizer is (" << std::setprecision(4) << rand_minimums[(nf - 1) * 2] << ',' << std::setprecision(4) << rand_minimums[(nf - 1) * 2 + 1] << ')' << std::endl;
		double glob = random_func(rand_minimums[(nf - 1) * 2], rand_minimums[(nf - 1) * 2 + 1]);
		std::cout << std::endl << "Computing function " << nf << "..." << std::endl;

		auto start = sc.now(); //start time for one function evaluating

		double UPB = gs.search(2,x,a,b,func);

		auto end = sc.now();
		gs.checkErrors();
		auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4)  << UPB << std::endl << \
			"Real global minimum: " << std::setprecision(4) << glob << std::endl \
		<< "Diff: " << std::setprecision(4) << 1.0*fabs(glob - UPB) << std::endl;
		fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
		fp << std::endl << std::endl;
	} /* for nf*/

	std::cout << std::endl;

	auto aend = sc.now();
	auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

	fp << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << 's' << std::endl;

	fp.close();
	delete[]a; delete[]b; delete[]x;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
	return 0;
}






