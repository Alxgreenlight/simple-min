// EnergyCluster.cpp: определяет точку входа для консольного приложения.
//
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include "../solver/gridsolver.hpp"

int N = 2; //number of molecules
int m = 2; //number of dimensions
double maxr = 12.0;
double nullcase = 1.0 / (maxr * maxr * maxr);

double energy(const double* x) {
	double fullsum = 0.0;
	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) {
			double r = 0.0;
			for (int k = 0; k < m; k++) {
				r += (x[m*j + k] - x[m*i + k])*(x[m*j + k] - x[m*i + k]);
			}
			if (r < 1e-308) {
				fullsum = nullcase * nullcase - 2 * nullcase;
				return fullsum;
			}
			r = 1.0 / (r * r * r);
			fullsum += r * r - 2 * r;
		}
	}
	return fullsum;
}

int main()
{
	std::chrono::steady_clock sc; //for runtime measurement
	GridSolverOMP<double> gs;
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
	gs.setparams(nodes, eps);
	std::ofstream fp;
	fp.open("paralleltest.txt", std::ios::out);
	for (int G = 2; G < 5; G++) {
		N = G;
		std::cout << "Starting for " << G << " cluster" << std::endl;
		double *a = new double[m * N];
		double *b = new double[m * N];
		double *x = new double[m * N];
		for (int i = 0; i < m * N; i++) {
			a[i] = -1.0;
			b[i] = 1.0;
		}
		auto astart = sc.now();
		double upb = gs.search(m * N, x, a, b, energy);
		auto aend = sc.now();
		auto atime_span = std::chrono::duration_cast<std::chrono::milliseconds>(aend - astart);
		gs.checkErrors();
		fp << "Found minimum energy: " << upb << std::endl;
		fp << "at coordinates:" << std::endl;
		for (int i = 0; i < N; i++) {
			for (int k = 0; k < m; k++) {
				fp << x[i * m + k] << "    ";
			}
			fp << std::endl;
		}
		unsigned long long int fevals;
		unsigned long int iters;
		gs.getinfo(fevals, iters);
		fp << "Consumed evaluations: " << fevals << ", iterations: " << iters << std::endl;
		fp << "Computing completed in " << atime_span.count() << " millicseconds" << std::endl;
		delete[]a;
		delete[]b;
		delete[]x;
	}
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
    return 0;
}

