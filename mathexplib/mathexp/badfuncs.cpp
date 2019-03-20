#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>
#include <iomanip>
#include <iterator>
#include "../testfuncs/benchmarks.hpp"
#include "gridsolveromp.hpp"


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
	double* a = new double[dim];
	double* b = new double[dim];
	double* x = new double[dim];
	for (int i = 0; i < dim; i++) {
		a[i] = bm.getBounds()[i].first;
		b[i] = bm.getBounds()[i].second;
	}
	GridSolverOMP<double> gs;
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
	int er;
	gs.getinfo(fevals, iters, er);
	if (er == -4) {
		fp << "Set time was out, computing incomplete" << std::endl;
	}
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
	if (bm.getDim() > 5) {
		fp << "INSUFFICIENT MEMORY" << std::endl;
	}
	else
		findMin(bm, fp, nodes, eps);
	fp << "****************************************" << std::endl << std::endl;
	return;
}

int main(int argc, char* argv[]) {
	std::ofstream fp;
	fp.open("results.txt", std::ios::app);
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

	for (int nod = 35; nod <= 35; nod += 30) {

		auto astart = sc.now();
		nodes = nod;
		fp << std::endl << "USING " << nod << " NODES" << std::endl << std::endl;

		/*Alpine1Benchmark<double> ap1(3);
		ptr = &ap1;
		testBench(ap1, fp, nodes, eps);
		std::cout << "|";
		BiggsEXP4Benchmark<double> bg4;
		ptr = &bg4;
		testBench(bg4, fp, nodes, eps);
		std::cout << "|";
		BiggsEXP5Benchmark<double> bg5;
		ptr = &bg5;
		testBench(bg5, fp, nodes, eps);
		std::cout << "|";*/
		BiggsEXP6Benchmark<double> bg6;
		ptr = &bg6;
		testBench(bg6, fp, nodes, eps);
		std::cout << "|";
		BrownBenchmark<double> bb(3);
		ptr = &bb;
		testBench(bb, fp, nodes, eps);
		std::cout << "|";
		Bukin6Benchmark<double> bu6;
		ptr = &bu6;
		testBench(bu6, fp, nodes, eps);
		std::cout << "|";
		ChungReynoldsBenchmark<double> cg(3);
		ptr = &cg;
		testBench(cg, fp, nodes, eps);
		std::cout << "|";
		ColvilleBenchmark<double> cl;
		ptr = &cl;
		testBench(cl, fp, nodes, eps);
		std::cout << "|";
		DolanBenchmark<double> dl;
		ptr = &dl;
		testBench(dl, fp, nodes, eps);
		std::cout << "|";
		GriewankBenchmark<double> gw(3);
		ptr = &gw;
		testBench(gw, fp, nodes, eps);
		std::cout << "|";
		Hartman6Benchmark<double> ht;
		ptr = &ht;
		testBench(ht, fp, nodes, eps);
		std::cout << "|";
		Langerman5Benchmark<double> lm5;
		ptr = &lm5;
		testBench(lm5, fp, nodes, eps);
		std::cout << "|";
		Mishra8Benchmark<double> ms;
		ptr = &ms;
		testBench(ms, fp, nodes, eps);
		std::cout << "|";
		Mishra9Benchmark<double> ms9;
		ptr = &ms9;
		testBench(ms9, fp, nodes, eps);
		std::cout << "|";
		PinterBenchmark<double> pt(3);
		ptr = &pt;
		testBench(pt, fp, nodes, eps);
		std::cout << "|";
		PowellSingular2Benchmark<double> pw(8);
		ptr = &pw;
		testBench(pw, fp, nodes, eps);
		std::cout << "|";
		Price4Benchmark<double> pc;
		ptr = &pc;
		testBench(pc, fp, nodes, eps);
		std::cout << "|";
		Scahffer1Benchmark<double> sch;
		ptr = &sch;
		testBench(sch, fp, nodes, eps);
		std::cout << "|";
		SchafferF6Benchmark<double> sf6(3);
		ptr = &sf6;
		testBench(sf6, fp, nodes, eps);
		std::cout << "|";
		SchmidtVettersBenchmark<double> sv;
		ptr = &sv;
		testBench(sv, fp, nodes, eps);
		std::cout << "|";
		StrechedVSineWaveBenchmark<double> svs(3);
		ptr = &svs;
		testBench(svs, fp, nodes, eps);
		std::cout << "|";
		Trid10Benchmark<double> td;
		ptr = &td;
		testBench(td, fp, nodes, eps);
		std::cout << "|";
		Trid6Benchmark<double> td6;
		ptr = &td6;
		testBench(td6, fp, nodes, eps);
		std::cout << "|";
		Trigonometric2Benchmark<double> tg(3);
		ptr = &tg;
		testBench(tg, fp, nodes, eps);
		std::cout << "|";
		/*TripodBenchmark<double> tp;
		ptr = &tp;
		testBench(tp, fp, nodes, eps);
		std::cout << "|";
		WhitleyBenchmark<double> wh(3);
		ptr = &wh;
		testBench(wh, fp, nodes, eps);
		std::cout << "|";
		XinSheYang2Benchmark<double> xsz(3);
		ptr = &xsz;
		testBench(xsz, fp, nodes, eps);
		std::cout << "|";*/
		auto aend = sc.now();
		auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

		fp << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
		fp << "Maximum difference in this set: " << maxdiff << std::endl;
		std::cout << "Total evaluation time: " << atime_span.count() << 's' << std::endl;
	}

	fp.close();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
}
