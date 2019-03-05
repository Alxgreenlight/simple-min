#ifndef GRIDSOLVER_HPP
#define GRIDSOLVER_HPP

#include <math.h>
#include <algorithm>
#include <limits>
#include "common/bbsolver.hpp"


template <class T>
struct Partitions;


template <class T> class GridSolver : public BlackBoxSolver <T> {
public:

	/**
	* Constructor
	*/
	GridSolver() {
		errcode = 0;
		eps = 100;
		nodes = 2;
	}

	virtual void setparams(int n, T e) {
		eps = e;
		nodes = n;
	}

	virtual T search(int n, T* x, const T * const a, const T * const b, const std::function<T(const T * const)> &f) {
		fevals = 0;
		iters = 0;
		Partitions<T> P(n), P1(n);
		UPB = std::numeric_limits<T>::max();	//upper bound
		if (f == nullptr) {
			errcode = -1;
			return -100e10;
		}

		P.add(a, b);	//add first hyperinterval

		T *a1, *b1;
		//if hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates
		a1 = new T[n];
		b1 = new T[n];
		T *xs = new T[n];

		if ((a1 == nullptr) || (b1 == nullptr)) {
			std::cerr << "Problem with memory allocation" << std::endl;
			errcode = -2;
			return -100e10;
		}

		while (P.size != 0) {	//each hyperinterval can be subdivided or pruned (if non-promisiable or fits accuracy)
			unsigned int parts = P.size;
			iters += parts;

			for (unsigned int i = 0; i < parts; i++) {
				T lUPB, lLOB, ldeltaL;
				GridEvaluator(n,P[i].a, P[i].b, xs, &lUPB, &lLOB, &ldeltaL, f);
				P[i].LocLO = lLOB;
				P[i].LocUP = lUPB;
				P[i].deltaL = ldeltaL;
				UpdateRecords(lUPB, x, xs, n);
			}

			for (unsigned int i = 0; i < parts; i++) {
				if (!((P[i].LocLO > (UPB - eps)) || (P[i].deltaL < eps))) {
					int choosen = ChooseDim(n, P[i].a, P[i].b);

					for (int j = 0; j < n; j++) { //make new edges for 2 new hyperintervals:
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
		delete[]a1; delete[]b1; delete[]xs;
		P.erase();
		return UPB;
	}
	virtual void checkErrors() {
		switch (errcode) {
		case -1:
			std::cerr << "Pointer to computing function (*compute) have not been initialized!" << std::endl;
			break;
		case -2:
			std::cerr << "Sorry, amount of RAM on your device insufficient to solve this task, please upgrade :)" << std::endl;
			break;
		default:
			break;
		}
	}
	virtual void getinfo(unsigned long long int &evs, unsigned long int &its) {
		evs = fevals;
		its = iters;
	}

protected:
	unsigned long long fevals;
	unsigned long iters;
	T eps;
	int errcode;
	int nodes;
	T UPB, LOB;
	virtual double getR(const T delta) { //r value for formula for estimating Lipschitz constant
		return static_cast<double>(exp(delta));
	}
	virtual void UpdateRecords(const T LU, T* x, const T *xs, int n) { //UPB - global upper bound, LU - local value of upper bound on current hyperinterval
		if (LU < UPB) {
			UPB = LU;
			for (int i = 0; i < n; i++) {
				x[i] = xs[i];
			}
		}
	}

	virtual int ChooseDim(int dim, const T *a, const T *b) { //selection of coordinate for hyperinterval division
		T max = std::numeric_limits<T>::min(), cr;
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

	virtual void GridEvaluator(int dim, const T *a, const T *b, T* xfound, T *Frp, T *LBp, T *dL, const std::function<T(const T * const)> &compute) {
		T Fr = std::numeric_limits<T>::max(), L = std::numeric_limits<T>::min(), delta = std::numeric_limits<T>::min(), LB;
		double R;
		T *step, *x, *Fvalues;
		step = new T[dim]; //step of grid in every dimension
		x = new T[dim]; //algorithm now needs to evaluate function in two adjacent point at the same time
		if ((step == nullptr) || (x == nullptr)) {
			errcode =-2;
			return;
		}
		for (int i = 0; i < dim; i++) {
			step[i] = fabs(b[i] - a[i]) / (nodes - 1);
			delta = step[i] > delta ? step[i] : delta;
		}
		delta = delta * 0.5 * dim;
		R = getR(delta);	//value of r
		int allnodes = static_cast<int>(pow(nodes, dim)), node;
		Fvalues = new T[allnodes];
		if (Fvalues == nullptr) {
			errcode = -2;
			return;
		}
		for (int j = 0; j < allnodes; j++) {
			int point = j;
			for (int k = dim - 1; k >= 0; k--) {
				int t = point % nodes;
				point = (int)(point / nodes);
				x[k] = a[k] + t * step[k];
			}
			T rs = compute(x);
			Fvalues[j] = rs;
			if (rs < Fr) {
				Fr = rs;
				node = j;
			}
		}
		fevals += allnodes;
		for (int k = dim - 1; k >= 0; k--) {
			int t = node % nodes;
			node = (int)(node / nodes);
			xfound[k] = a[k] + t * step[k];
		}

		for (int j = 0; j < allnodes; j++) {
			int neighbour;
			for (int k = 0; k < dim; k++) {
				int board = (int)pow(nodes, k + 1);
				neighbour = j + board / nodes;
				if ((neighbour < allnodes) && ((j / board) == (neighbour / board))) {
					T loc = fabs(Fvalues[j] - Fvalues[neighbour]) / step[dim - 1 - k];
					L = loc > L ? loc : L;
				}
			}
		}
		//final calculation
		LB = R * L * delta;
		*dL = LB;
		LB = Fr - LB;
		delete[]step; delete[]x;
		delete[]Fvalues;
		*Frp = Fr;
		*LBp = LB;
	}
};

template <class T>
struct part {	//description of hyperinterval
	T *a = nullptr;
	T *b = nullptr;
	T LocUP, LocLO, deltaL;
};

template <class T>
struct Partitions {		//all hyperintervals obtained on current step of algorithm
	const int chunk = 16;
	int size, dim;			//like std::vector, but in C language with only needed methods
	int cur_alloc;		//memory already allocated for array of struct part
	struct part<T>* base;
	Partitions(int n) {
		dim = n;
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
	int add(const T* toa, const T *tob) {
		if (cur_alloc == 0) {
			base = (struct part<T>*)malloc(chunk * sizeof(struct part<T>));
			if (base) {
				cur_alloc = chunk;
			}
			else {
				std::cerr << "Error in allocation procedure" << std::endl;
				erase();
				return -1;
			}
		}
		if (size == cur_alloc) {
			base = (struct part<T>*)realloc(base, (cur_alloc + chunk) * sizeof(struct part<T>));
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
	part<T> & operator[](int n) {
		if (n > size - 1) {
			std::cerr << "Out of size" << std::endl;
		}
		return base[n];
	}
	Partitions<T> & operator=(Partitions<T>& P) {
		erase();
		for (int i = 0; i < P.size; i++) {
			add(P[i].a, P[i].b);
		}
		return (*this);
	}
};

#endif /* GRIDSOLVER_HPP */
