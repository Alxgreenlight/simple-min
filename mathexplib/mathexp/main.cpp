/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

/*
* File:   tutorial.cpp
*
* Created on November 20, 2017, 9:13 AM
*/
#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>
#include <iomanip>
#include "../testfuncs/benchmarks.hpp"
#include "../../solver_parallel/solver_omp.h"


using BM = Benchmark<double>;
BM *ptr;
std::chrono::steady_clock sc;
double maxdiff = std::numeric_limits<double>::min();
double calc(double *x)
{
	std::vector<double> xc;
	xc.assign(x, x + solver::dim);
	double v = ptr->calcFunc(xc);
	return v;
}
void findMin(const BM& bm, std::ofstream &fp) {
	solver::dim = bm.getDim();
	int r = 0;
	double* a = new double[solver::dim];
	double* b = new double[solver::dim];
	for (int i = 0; i < solver::dim; i++) {
		a[i] = bm.getBounds()[i].first;
		b[i] = bm.getBounds()[i].second;
	}
	solver::compute = calc;
	solver::glob = bm.getGlobMinY();
	auto start = sc.now();
	r = solver::min_search(a, b);
	if (r) {
		solver::checkErrors(r);
		delete[]a; delete[]b;
		fp.close();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return;
	}
	auto end = sc.now();
	auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	double diff = 1.0*fabs(solver::glob - solver::UPB);
	maxdiff = diff > maxdiff ? diff : maxdiff;
	fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4) << solver::UPB << std::endl << \
		"Real global minimum: " << std::setprecision(4) << solver::glob << std::endl << "Lower bound: " << std::setprecision(4) << \
		solver::LOB << std::endl << "Diff: " << std::setprecision(4) << diff << std::endl;
	fp << "Function evaluations: " << solver::fevals << " in " << solver::iters << " iterations" << std::endl;
	fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
	fp << std::endl << std::endl;
	delete[]a; delete[]b;
	return;
}

void testBench(const BM& bm, std::ofstream& fp) {
	fp << "*************Testing benchmark**********" << std::endl;
	fp << bm;
	findMin(bm, fp);
	fp << "****************************************" << std::endl << std::endl;
	return;
}

int main() {
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

	std::cout << "Set max dim for evaluate" << std::endl << "0 = unlimited (full benchmark)" << std::endl;
	int maxdim;
	std::cin >> maxdim;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
	}
	if (!maxdim) maxdim = std::numeric_limits<int>::max();
	auto astart = sc.now();
	/*ZakharovBenchmark<double> zb(3);
	ptr = &zb;
	testBench(zb, fp);
	std::cout << '|';*/
	Benchmarks<double> tests;
	for (auto bm : tests) {
		ptr = &(*bm);
		if (bm->getDim() <= maxdim) {
			testBench(*bm, fp);
			std::cout << '|';
		}
	}
	std::cout << std::endl;
	auto aend = sc.now();
	auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

	fp << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	fp << "Maximum difference in this set: " << maxdiff << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	fp.close();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
}
