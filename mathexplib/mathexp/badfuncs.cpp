#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>
#include <iomanip>
#include <iterator>
#include "../testfuncs/benchmarks.hpp"
#include "../../solver/gridsolver.hpp"


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
void findMin(const BM& bm, std::ofstream &fp, int nodes, double eps, unsigned long long int &time) {
	int dim = bm.getDim();
	double* a = new double[dim];
	double* b = new double[dim];
	double* x = new double[dim];
	for (int i = 0; i < dim; i++) {
		a[i] = bm.getBounds()[i].first;
		b[i] = bm.getBounds()[i].second;
	}
	GridSolverOMP<double> gs;
	gs.setparams(nodes, eps);
	double UPB;
	auto start = sc.now();
	UPB = gs.search(dim, x, a, b, calc);
	auto end = sc.now();
	auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	gs.checkErrors();
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
	time = time_span.count();
	fp << std::endl << std::endl;
	fp << std::endl << std::endl;
	delete[]a; delete[]b;
	return;
}

void testBench(const BM& bm, std::ofstream& fp, int nodes, double eps, unsigned long long int &time) {
	fp << "*************Testing benchmark**********" << std::endl;
	fp << bm;
		findMin(bm, fp, nodes, eps,time);
	fp << "****************************************" << std::endl << std::endl;
	return;
}

int main(int argc, char* argv[]) {
	std::ofstream fp;
	fp.open("parallel_results.txt", std::ios::app);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}

	int nodes;
	double eps;
	unsigned long long int time, fulltime = 0;

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

		fp << std::endl << "USING " << nodes << " NODES" << std::endl << std::endl;

		Ackley3Benchmark<double> ack;
		ptr = &ack;
		testBench(ack, fp, nodes, eps, time);
		fulltime += time;
		std::cout << "|";
		RosenbrockBenchmark<double> rs(3);
		ptr = &rs;
		testBench(rs, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		BealeBenchmark<double> bl;
		ptr = &bl;
		testBench(bl, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		GoldsteinPriceBenchmark<double> gsp;
		ptr = &gsp;
		testBench(gsp, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		BoothBenchmark<double> bb;
		ptr = &bb;
		testBench(bb, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		MatyasBenchmark<double> mb;
		ptr = &mb;
		testBench(mb, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		HimmelblauBenchmark<double> hb;
		ptr = &hb;
		testBench(hb, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		SphereBenchmark<double> sph(3);
		ptr = &sph;
		testBench(sph, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		EggHolderBenchmark<double> eg;
		ptr = &eg;
		testBench(eg, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;
		StyblinskiTangBenchmark<double> stb;
		ptr = &stb;
		testBench(stb, fp, nodes, eps, time);
		std::cout << "|";
		fulltime += time;

		fp << "Total evaluation time: " << fulltime << 'ms' << std::endl;
		fp << "Maximum difference in this set: " << maxdiff << std::endl;
		std::cout << "Total evaluation time: " << static_cast<long double>(fulltime)/1000.0 << 's' << std::endl;

	fp.close();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
}
