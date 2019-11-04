#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include "solver/R_Optim_Pure_Parallel.hpp"
#include "OneDim_temp/miniOneDimFuncs.hpp"
#include "solver/UltraEstim.hpp"

std::chrono::steady_clock sc;
std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();
onedimopt::BasicFunc<double> *ptr;

double func(const double* x) {	/* wrapper for function providing to solver */
	return ptr->calculate(x[0]);
}

int main()
{
	std::ofstream out;
	out.open("Results.txt", std::ios::out);
	onedimopt::BMsOneDim<double> bms;
	double UpperBound;
	double eps;
	unsigned long long int fevals;
	unsigned long int iters;
	double xfound;
	double L;

	std::cout << "Set accuracy" << std::endl << "It can significantly affect on performance!" << std::endl;
	std::cin >> eps;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
		std::cin >> eps;
	}
	for (auto elem : bms) {
		ptr = &(*elem);
		std::cout << elem->detDesc() << std::endl;
		std::cout << "With bounds [" << elem->LBound() << "," << elem->RBound() << "]" << std::endl;
		std::cout << "In x = " << elem->RealminX() <<  " and real y = " << elem->RealminY() << std::endl;
		try {
			PPrOptimizer<double> rOpt;
			L_accurate<double> getL;
			rOpt.init(1, eps);
			auto start = sc.now();
			double lb = elem->LBound(), rb = elem->RBound();
			UpperBound = rOpt.search(1, &xfound, &lb, &rb, func);
			auto end = sc.now();
			auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			double diff = 1.0*fabs(elem->RealminY() - UpperBound);
			rOpt.getInfo(fevals, iters, L);
			std::cout << std::fixed << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(10) << UpperBound << std::endl << \
				"At " << xfound << ", with L = " << L << std::endl;
			std::cout << std::fixed << "Diff: " << std::setprecision(10) << diff << std::endl;
			std::cout << "Function evaluations: " << fevals << " in " << iters << " iterations" << std::endl;
			std::cout << "Evaluation time: " << time_span.count() << "ms" << std::endl;
			rOpt.clear();
			double L = getL.ultraoptimizer_1(100000000, lb, rb, xfound, UpperBound, func);
			std::cout << "Very accurate estimation: L = " << L << ", Upper bound = " << UpperBound << ", at " << xfound << std::endl;
		}
		catch (std::exception& e) {
			std::cout << e.what() << std::endl;
			out.close();
			return -1;
		}
		std::cout << std::endl << std::endl;
	}
	out.close();
	return 0;
}


