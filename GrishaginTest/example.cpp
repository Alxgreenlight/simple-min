#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "../solver/solver.h"
extern "C" {
#include "vagris.h"
}

double *a, *b; //sets the search area for task

std::chrono::steady_clock sc; //for runtime measurement

int main()
{
	
	int nf;
	std::ofstream fp;

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

	for (int i = 0; i < solver::dim; i++) {
		a[i] = 0.0;
		b[i] = 1.0;
	} //set search area for all tests to [0..1] [0..1] (like in tutorial for this tests)

	solver::compute = random_func; //random_func returns value of function in x1,x2 (from vagris.h)

	auto astart = sc.now(); //start time for full set of tests
	for (nf = 1; nf <= 100; nf++)
	{
		set_random(nf);

		fp << std::endl << "Function number " << nf;
		fp << std::endl << "Its global minimizer is (" << std::setprecision(4) << rand_minimums[(nf - 1) * 2] << ',' << std::setprecision(4) << rand_minimums[(nf - 1) * 2 + 1] << ')' << std::endl;
		solver::glob = random_func(rand_minimums[(nf - 1) * 2], rand_minimums[(nf - 1) * 2 + 1]);
		std::cout << std::endl << "Computing function " << nf << "..." << std::endl;

		auto start = sc.now(); //start time for one function evaluating

		int r = solver::min_search(a, b);

		if (r) {
			solver::checkErrors(r);
			fp.close();
			delete[]a; delete[]b;
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cin.get();
			return -2;
		}

		auto end = sc.now();
		auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4)  << solver::UPB << std::endl << \
			"Real global minimum: " << std::setprecision(4) << solver::glob << std::endl << "Lower bound: " << std::setprecision(4) <<\
			solver::LOB << std::endl << "Diff: " << std::setprecision(4) << 1.0*fabs(solver::glob - solver::UPB) << std::endl;
		fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
		fp << std::endl << std::endl;
	} /* for nf*/

	std::cout << std::endl;

	auto aend = sc.now();
	auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

	fp << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << 's' << std::endl;

	fp.close();
	delete[]a; delete[]b;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
	return 0;
}






