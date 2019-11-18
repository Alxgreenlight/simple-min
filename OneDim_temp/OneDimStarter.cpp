#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include "solver/R_Optim_Pure_Parallel.hpp"
#include "OneDim_temp/miniOneDimFuncs.hpp"
#include "solver/UltraEstim.hpp"
#include "util/helper.hpp"

std::chrono::steady_clock sc;
std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();
onedimopt::BasicFunc<double> *ptr;
const int dim = 1;

double func(const double* x) {	/* wrapper for function providing to solver */
	return ptr->calculate(x[0]);
}

int main(int argc, char **argv)
{
	std::ofstream out;
	onedimopt::BMsOneDim<double> bms;
	double UpperBound;
	double eps;
	unsigned long long int fevals;
	unsigned long int iters;
	double xfound;
	double L;
	int todo = 5, cur = 0;

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

	out.open("Results.txt", std::ios::out);
	helper::progress_bar(cur, todo);

	for (auto elem : bms) {
		ptr = &(*elem);
		out << elem->detDesc() << std::endl;
		out << "With bounds [" << elem->LBound() << "," << elem->RBound() << "]" << std::endl;
		out << "In x = " << elem->RealminX() <<  " and real y = " << elem->RealminY() << std::endl;
		try {
			PPrOptimizer<double> rOpt;
			L_accurate<double> getL;
			rOpt.init(dim, eps);
			double lb = elem->LBound(), rb = elem->RBound();
			auto start = sc.now();
			UpperBound = rOpt.search(1, &xfound, &lb, &rb, func);
			auto end = sc.now();
			auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			double diff = 1.0*fabs(elem->RealminY() - UpperBound);
			rOpt.getInfo(fevals, iters, L);
			out << std::fixed << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(10) << UpperBound << std::endl << \
				"At " << xfound << ", with L = " << L << std::endl;
			out << std::fixed << "Diff: " << std::setprecision(10) << diff << std::endl;
			out << "Function evaluations: " << fevals << " in " << iters << " iterations" << std::endl;
			out << "Evaluation time: " << time_span.count() << "ms" << std::endl;
			rOpt.clear();
			//double L = getL.ultraoptimizer_1(100000000, lb, rb, xfound, UpperBound, func);
			out << "Very accurate estimation: L = " << L << ", Upper bound = " << UpperBound << ", at " << xfound << std::endl;
			cur++;
			helper::progress_bar(cur, todo);
		}
		catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			out.close();
			return -1;
		}
	}
	out.close();
	return 0;
}


