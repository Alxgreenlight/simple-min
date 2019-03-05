#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>
#include <iomanip>
#include <iterator>
#include "../testfuncs/benchmarks.hpp"
#include "gridsolveronedim.hpp"


using BM = Benchmark<double>;
BM *ptr;
std::chrono::steady_clock sc;
double maxdiff = std::numeric_limits<double>::min();
double calc(const double *x)
{
	std::vector<double> xc;
	xc.assign(x, x + ptr->getDim());
	double v = ptr->calcFunc(xc);
	return v;
}
void findMin(const BM& bm, std::ofstream &fp, int nodes, double eps) {
	int dim = bm.getDim();
	int r = 0;
	double* a = new double[dim];
	double* b = new double[dim];
	double* x = new double[dim];
	for (int i = 0; i < dim; i++) {
		a[i] = bm.getBounds()[i].first;
		b[i] = bm.getBounds()[i].second;
	}
	GridSolverOneDim<double> gs;
	gs.setparams(nodes, eps);
	auto start = sc.now();
	double UPB;
	UPB = gs.search(dim, x, a, b, calc);
	gs.checkErrors();
	auto end = sc.now();
	auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	double diff = 1.0*fabs(bm.getGlobMinY() - UPB);
	maxdiff = diff > maxdiff ? diff : maxdiff;
	unsigned long long int fevals;
	unsigned long int iters;
	gs.getinfo(fevals, iters);
	fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4) << UPB << std::endl << \
		"At [";
	std::copy(x, x + dim, std::ostream_iterator<double>(fp, " "));
	fp << "]" << std::endl << "Real global minimum: " << std::setprecision(4) << bm.getGlobMinY() << std::endl << std::endl << "Diff: " << std::setprecision(4) << diff << std::endl;
	fp << "Function evaluations: " << fevals << " in " << iters << " iterations" << std::endl;
	fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
	fp << std::endl << std::endl;
	fp << std::endl << std::endl;
	delete[]a; delete[]b;
	return;
}

void testBench(const BM& bm, std::ofstream& fp, int nodes, double eps) {
	fp << "*************Testing benchmark**********" << std::endl;
	fp << bm;
	findMin(bm, fp, nodes, eps);
	fp << "****************************************" << std::endl << std::endl;
	return;
}

int main(int argc, char* argv[]) {
	std::ofstream fp;
	fp.open("results.txt", std::ios::out);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}

	int nodes;
	double eps;

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

	/*int maxdim, mindim = 0;
	if (argc == 1) {
		std::cout << "Set max dim for evaluate" << std::endl << "0 = unlimited (full benchmark)" << std::endl;
		std::cin >> maxdim;
		while (std::cin.fail()) {
			std::cerr << "Please, repeat input" << std::endl;
		}
		if (!maxdim) maxdim = std::numeric_limits<int>::max();
	}
	try {
		if (argc == 2) maxdim = std::stoi(argv[1]);
		if (argc == 3) {
			mindim = std::stoi(argv[1]);
			maxdim = std::stoi(argv[2]);
		}
	}
	catch (std::exception const & e)
	{
		std::cerr << "error : " << e.what() << std::endl;
	}*/

	auto astart = sc.now();
	/*ZakharovBenchmark<double> zb(3);
	ptr = &zb;
	testBench(zb, fp);
	std::cout << '|';*/
	Benchmarks<double> tests;
	for (auto bm : tests) {
		ptr = &(*bm);
		//if ((bm->getDim() <= maxdim)&&(bm->getDim() >= mindim)) {
		if (bm->getDim() == 1) {
			testBench(*bm, fp, nodes, eps);
			std::cout << '|';
		}
		//}
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
