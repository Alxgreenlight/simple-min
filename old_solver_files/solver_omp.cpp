#include <cfloat>
#include <cmath>
#include <iostream>
#include "solver_omp.h"
#include "omp.h"

namespace solver {
	double eps = 100; //accuracy (100 is default value)
	double UPB, LOB, deltaL, glob = DBL_MAX; //obtained upped bound, lower bound, value of delta*L, and real global minimizer
	int nodes = 3, dim = 2, np = 0; //amount of nodes per dimension, amount of dimensions in task
	unsigned long long int fevals = 0;
	unsigned long int iters = 0;

	double(*compute)(double *x) = nullptr;

	struct part {	//description of hyperinterval
		double *a = nullptr;
		double *b = nullptr;
		double LocUP, LocLO, deltaL;
	};


	struct Partitions {		//all hyperintervals obtained on current step of algorithm
		const int chunk = 16;
		int size;			//like std::vector, but in C language with only needed methods
		int cur_alloc;		//memory already allocated for array of struct part
		struct part* base;
		Partitions() {
			size = 0;
			cur_alloc = 0;
			base = nullptr;
		}
		~Partitions() {
			for (int i = 0; i < size; i++) {
				free(base[i].a);
				free(base[i].b);
			}
			size = 0;
			if (cur_alloc != 0)
				free(base);
			cur_alloc = 0;;
		}
		void erase() {
			for (int i = 0; i < size; i++) {
				free(base[i].a);
				free(base[i].b);
			}
			size = 0;
			if (cur_alloc != 0)
				free(base);
			cur_alloc = 0;
		}
		int add(const double* toa, const double *tob) {
			if (cur_alloc == 0) {
				base = (struct part*)malloc(chunk * sizeof(struct part));
				if (base) {
					cur_alloc = chunk;
				}
				else {
					std::cerr << "Error alloc" << std::endl;
					erase();
					return -1;
				}
			}
			if (size == cur_alloc) {
				base = (struct part*)realloc(base, (cur_alloc + chunk) * sizeof(struct part));
				if (base) {
					cur_alloc += chunk;
				}
				else {
					std::cerr << "Error alloc" << std::endl;
					erase();
					return -1;
				}
			}
			base[size].a = (double*)malloc(dim * sizeof(double));
			base[size].b = (double*)malloc(dim * sizeof(double));
			if ((!base[size].a) || (!base[size].b)) {
				std::cerr << "Error alloc" << std::endl;
				erase();
				return -1;
			}
			for (int i = 0; i < dim; i++) {
				base[size].a[i] = toa[i];
				base[size].b[i] = tob[i];
			}
			size++;
			return 0;
		}
		part & operator[](int n) {
			if (n > size - 1) {
				std::cerr << "Out of size" << std::endl;
			}
			return base[n];
		}
		Partitions & operator=(Partitions& P) {
			erase();
			for (int i = 0; i < P.size; i++) {
				add(P[i].a, P[i].b);
			}
			return (*this);
		}
	};


	Partitions P, P1; //P - hyperintervals on current step, P1 - hyperintervals on next step

	double getR(const double delta) { //r value for formula for estimating Lipschitz constant
		return exp(delta);
	}

	void UpdateRecords(const double LU, const double LL) { //UPB - global upper bound, LU - local value of upper bound on current hyperinterval
		if (LU < UPB) {
			UPB = LU;
			LOB = LL;
		}
	}

	int ChooseDim(const double *a, const double *b) { //selection of coordinate for hyperinterval division
		double max = DBL_MIN, cr;
		int i, maxI = 0;
		for (i = 0; i < dim; i++) {
			cr = fabs(b[i] - a[i]);
			if (cr > max) {
				max = cr;
				maxI = i;
			}
		}
		return maxI;
	}

	void GridEvaluator(double *a, double *b, double *Frp, double *LBp, double *dL) {
		double Fr = DBL_MAX, L = DBL_MIN, R, LB, delta = DBL_MIN;
		double *step, *Fvalues;
		double *x;
		step = (double*)malloc(dim * sizeof(double)); //step of grid in every dimension
		if ((step == nullptr)) {
			*dL = -100;
			return;
		}
		for (int i = 0; i < dim; i++) {
			step[i] = fabs(b[i] - a[i]) / (nodes - 1);
			delta = step[i] > delta ? step[i] : delta;
		}
		delta = delta * 0.5 * dim;
		R = getR(delta);	//value of r
		int allnodes = (int)pow(nodes, dim);
		Fvalues = (double*)malloc(allnodes * sizeof(double));
		if (Fvalues == nullptr) {
			*dL = -100;
			return;
		}
#pragma omp parallel private(x) shared(dim, nodes,  step, a, Fvalues)
		{
			x = (double*)malloc(dim * sizeof(double));
#pragma omp for reduction(min:Fr)
				for (int j = 0; j < allnodes; j++) {
					int point = j;
					for (int k = dim - 1; k >= 0; k--) {
						int t = point % nodes;
						point = (int)(point / nodes);
						x[k] = a[k] + t * step[k];
					}
					double rs = compute(x);
					Fvalues[j] = rs;
					if (rs < Fr) {
						Fr = rs;
					}
				}
				free(x);
		}
		fevals += allnodes;

#pragma omp parallel for shared(dim,allnodes,Fvalues) reduction (max:L)
		for (int j = 0; j < allnodes; j++) {
			int neighbour;
			double loc = DBL_MIN;
			for (int k = 0; k < dim; k++) {
				int board = (int)pow(nodes, k + 1);
				neighbour = j + board / nodes;
				if ((neighbour < allnodes) && ((j / board) == (neighbour / board))) {
					loc = fabs(Fvalues[j] - Fvalues[neighbour]) / step[dim - 1 - k];
				}
				if (loc > L) {
					L = loc;
				}
			}
		}
		//final calculation
		LB = R * L * delta;
		*dL = LB;
		LB = Fr - LB;
		free(step);
		free(Fvalues);
		*Frp = Fr;
		*LBp = LB;
	}

	int min_search(const double *a, const double* b) {
		fevals = 0;
		iters = 0;
		UPB = DBL_MAX;	//upper bound
		if (glob == DBL_MAX) {	//check if uninitialized values
			return -1;
		}
		if (compute == nullptr) {
			return -2;
		}

		if (np <= 0) {
			np = omp_get_num_procs();
		}
		omp_set_dynamic(0);
		omp_set_num_threads(np);

		P.add(a, b);	//add first hyperinterval

		double *a1, *b1;
		//if hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates
		a1 = (double*)malloc(dim * sizeof(double));
		b1 = (double*)malloc(dim * sizeof(double));

		if ((a1 == nullptr) || (b1 == nullptr)) {
			std::cerr << "Problem with memory allocation" << std::endl;
			return -3;
		}

		while (P.size != 0) {	//each hyperinterval can be subdivided or pruned (if non-promisiable or fits accuracy)
			unsigned int parts = P.size, i;
			iters += parts;

			for (i = 0; i < parts; i++) {
				double lUPB, lLOB, ldeltaL;
				GridEvaluator(P[i].a, P[i].b, &lUPB, &lLOB, &ldeltaL);
				P[i].LocLO = lLOB;
				P[i].LocUP = lUPB;
				P[i].deltaL = ldeltaL;
				UpdateRecords(lUPB, lLOB);
			}

			for (i = 0; i < parts; i++) {
				if (!((P[i].LocLO >(UPB - eps)) || (P[i].deltaL < eps))) {
					int choosen = ChooseDim(P[i].a, P[i].b);

					for (int j = 0; j < dim; j++) { //make new edges for 2 new hyperintervals:
						if (j != choosen) {						//[a .. b1] [a1 .. b]
							a1[j] = P[i].a[j];			//where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ]
							b1[j] = P[i].b[j];			//and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ]
						}
						else {
							a1[j] = P[i].a[j] + fabs(P[i].b[j] - P[i].a[j]) / 2.0;
							b1[j] = a1[j];
						}
					}
					P1.add(P[i].a, b1);
					P1.add(a1, P[i].b);
				}
			}

			P = P1;
			P1.erase();
		}
		free(a1); free(b1);
		P.erase();
		return 0;
	}

	void checkErrors(int code) {
		switch (code) {
		case -1:
			std::cerr << "glob value have not been initialized!" << std::endl;
			break;
		case -2:
			std::cerr << "Pointer to computing function (*compute) have not been initialized!" << std::endl;
			break;
		case -3:
			std::cerr << "Sorry, amount of RAM on your device insufficient to solve this task, please upgrade :)" << std::endl;
			break;
		default:
			break;
		}
	}
}