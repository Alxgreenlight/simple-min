#include <cmath>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "solver/R_Optim_Pure_Parallel.hpp"
#include "util/helper.hpp"

int N = 2; //number of molecules
int m = 2; //number of dimensions
double maxr = 12.0;
double nullcase = 1.0 / (maxr * maxr * maxr);

double energy(const double *x)
{
	double fullsum = 0.0;
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			double r = 0.0;
			for (int k = 0; k < m; k++)
			{
				r += (x[m * j + k] - x[m * i + k]) * (x[m * j + k] - x[m * i + k]);
			}
			if (r < 1e-308)
			{
				fullsum = nullcase * nullcase - 2 * nullcase;
				return fullsum;
			}
			r = 1.0 / (r * r * r);
			fullsum += r * r - 2 * r;
		}
	}
	return fullsum;
}

int main(int argc, char **argv)
{
	std::chrono::steady_clock sc; //for runtime measurement
	std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();
	PPrOptimizer<double> rOpt;
	double eps, L;
	std::ofstream fp;
	unsigned long long int fevals;
	unsigned long int iters;

	if (argc < 2)
	{
		helper::help(benchlib);
		return 0;
	}

	eps = atof(argv[1]);
	if (eps < std::numeric_limits<double>::min())
	{
		std::cerr << "Accuracy defined incorrect, exit" << std::endl;
		return -1;
	}

	fp.open("results.txt", std::ios::out);

	int todo = 3, cur = 0;
	helper::progress_bar(cur, todo);
	try
	{
		for (int G = 2; G < 5; G++)
		{
			N = G;
			double *a = new double[m * N];
			double *b = new double[m * N];
			double *x = new double[m * N];
			for (int i = 0; i < m * N; i++)
			{
				a[i] = -1.0;
				b[i] = 1.0;
			}
			rOpt.init(m * N, eps);
			auto start = sc.now();
			double upb = rOpt.search(m * N, x, a, b, energy);
			auto end = sc.now();
			auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			atime_span += time_span;
			
			fp << "Found minimum energy: " << upb << std::endl;
			fp << "at coordinates: [";
			std::copy(x, x + N*m, std::ostream_iterator<double>(fp, ", "));
			fp << "]" << std::endl;
			rOpt.getInfo(fevals, iters, L);
			fp << "With L = " << L << std::endl;
			fp << "Function evaluations: " << fevals << " in " << iters << " iterations" << std::endl;
			fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
			delete[] a;
			delete[] b;
			delete[] x;
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << std::endl;
		fp.close();
		return -1;
	}
	fp << "Total evaluation time: " << atime_span.count() << std::endl;
	fp.close();
	return 0;
}
