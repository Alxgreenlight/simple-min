
#ifndef GRIDSOLVER_HPP
#define GRIDSOLVER_HPP

#include <cmath>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <vector>
#include <fstream>
#include "../common/bbsolver.hpp"

template <class T>
class part;


template <class T>
class GridSolver : public BlackBoxSolver <T> {
public:
	/**
	* Constructor
	*/
	GridSolver() {
		errcode = 0;
		eps = 100;
		nodes = 2;
		param = false;
	}

	/**
	* Set params of solver
	* @param n number of nodes per dimension
	* @param e required accuracy
	*/
	virtual void setparams(int n, T e) {
		eps = e;
		nodes = n;
		param = true;
	}

	/**
	* Search with grid solver
	* @param n number of task dimensions
	* @param x coordinates of founded minimum (retvalue)
	* @param a,b left/right bounds of search region
	* @param f pointer to function for which search minimum
	*/
	virtual T search(int n, T* x, const T * const a, const T * const b, const std::function<T(const T * const)> &f) {
		/* reset variables */
		fevals = 0;
		iters = 0;
		/* create 2 vectors */
		/* P contains parts (hyperintervals on which search must be performed */
		/* P1 temporary */
		std::vector<part<T>> P, P1;
		/* Upper bound */
		UPB = std::numeric_limits<T>::max();
		/* check if function pointer provided */
		if (f == nullptr) {
			errcode = -1;
			return UPB;
		}

		/* Add first hyperinterval */

		part<T> pt(n, a, b);

		try {
			P.emplace_back(pt);
		}
		catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			errcode = -2;
			return UPB;
		}

		/* If hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates */
		/* xs temporary array for storage local min coordinates */
		T *a1, *b1, *xs;
		try {
			a1 = new T[n];
			b1 = new T[n];
			xs = new T[n];
		}
		catch (std::bad_alloc& ba) {
			std::cerr << ba.what() << std::endl;
			errcode = -2;
			return UPB;
		}

		/* Each hyperinterval can be subdivided or pruned (if non-promisable or fits accuracy) */
		while (!P.empty()) {
			/* number of iterations on this step (BFS) */
			unsigned int parts = P.size();
			iters += parts;

			/* For all hyperintervals on this step perform grid search */
			for (unsigned int i = 0; i < parts; i++) {
				/* local values of upper and lower bounds, value of delta*L (Lipshitz const) */
				T lUPB, lLOB, ldeltaL;
				T* ta = P[i].a, *tb = P[i].b;
				GridEvaluator(n, ta, tb, xs, &lUPB, &lLOB, &ldeltaL, f);
				P[i].LocLO = lLOB;
				P[i].LocUP = lUPB;
				P[i].deltaL = ldeltaL;
				/* remember new results if less then previous */
				UpdateRecords(lUPB, x, xs, n);
			}

			/* Choose which hyperintervals should be subdivided */
			for (unsigned int i = 0; i < parts; i++) {
				/* Subdivision criteria */
				if (P[i].LocLO < (UPB - eps)) {
					/* If subdivide, choose dimension (the longest side) */
					int choosen = ChooseDim(n, P[i].a, P[i].b);

					/* Make new edges for 2 new hyperintervals */
					for (int j = 0; j < n; j++) {
						if (j != choosen) {			/* [a .. b1] [a1 .. b] */
							a1[j] = P[i].a[j];		/* where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ] */
							b1[j] = P[i].b[j];		/* and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ] */
						}
						else {
							a1[j] = P[i].a[j] + fabs(P[i].b[j] - P[i].a[j]) / 2.0;
							b1[j] = a1[j];
						}
					}
					/* Add 2 new hyperintervals, parent HI no longer considered */
					try {
						pt.change(P[i].a, b1);
						P1.emplace_back(pt);
					}
					catch (std::exception& e) {
						std::cerr << e.what() << std::endl;
						errcode = -2;
						return UPB;
					}
					try {
						pt.change(a1, P[i].b);
						P1.emplace_back(pt);
					}
					catch (std::exception& e) {
						std::cerr << e.what() << std::endl;
						errcode = -2;
						return UPB;
					}
				}
			}

			P.clear();
			P.swap(P1);
		}
		delete[]a1; delete[]b1; delete[]xs;
		P.clear();
		return UPB;
	}

	/**
	* Check if there are errors during search process
	*/
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
		if (!param) {
			std::cerr << "Parameters have not been initialized, solved with default params:" << std::endl << \
				"        nodes = 2; eps = 100" << std::endl;
		}
	}

	/**
	* get number of consumed func evaluations and algorithm iterations after search completed
	* @param evs (retvalue) number of function evaluations
	* @param its (retvalue) number of algorithm iterations
	*/
	virtual void getinfo(unsigned long long int &evs, unsigned long int &its) {
		evs = fevals;
		its = iters;
	}

protected:

	unsigned long long fevals; /* number of function evaluations that search consumed */
	unsigned long iters; /* nember of algorithm itaretions that search consumed */
	T eps;	/* required accuracy */
	int errcode, nodes;	/* internal varibale for handlig errors and number of nodes per dimension */
	bool param;		/* indicate if nodes and eps was initialized */
	T UPB, LOB;		/* obtained upper bound and lower bound */

	/* Get R (reliable coefficient) for the corresponding step lenght*/
	virtual double getR(const T delta) {
		/* The value depends on step lenght (test formula, may be changed) */
		return static_cast<double>(exp(delta));
	}

	/* Update the current record and its coordinates in accordance with new results obtained on some hyperinterval*/
	virtual void UpdateRecords(const T LU, T* x, const T *xs, int n) {
		if (LU < UPB) {
			UPB = LU;
			for (int i = 0; i < n; i++) {
				x[i] = xs[i];
			}
		}
	}

	/* Select dimension for subdivide hyperinteral (choose the longest side) */
	virtual int ChooseDim(int dim, const T *a, const T *b) {
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
		try {
			/* step of grid in every dimension */
			step = new T[dim];
			x = new T[dim];
		}
		catch (std::bad_alloc& ba) {
			std::cerr << ba.what() << std::endl;
			errcode = -2;
			*Frp = Fr;
			*LBp = L;
			*dL = delta;
			return;
		}

		for (int i = 0; i < dim; i++) {
			step[i] = fabs(b[i] - a[i]) / (nodes - 1);
			delta = step[i] > delta ? step[i] : delta;
		}
		delta = delta * 0.5 * dim;
		R = getR(delta);
		int allnodes = static_cast<int>(pow(nodes, dim)), node = 0;
		try {
			Fvalues = new T[allnodes];
		}
		catch (std::bad_alloc& ba) {
			std::cerr << ba.what() << std::endl;
			errcode = -2;
			*Frp = Fr;
			*LBp = L;
			*dL = delta;
			return;
		}

		/* Calculate and cache the value of the function in all points of the grid */
		for (int j = 0; j < allnodes; j++) {
			int point = j;
			for (int k = dim - 1; k >= 0; k--) {
				int t = point % nodes;
				point = (int)(point / nodes);
				x[k] = a[k] + t * step[k];
			}
			T rs = compute(x);
			Fvalues[j] = rs;
			/* also remember minimum value across the grid */
			if (rs < Fr) {
				Fr = rs;
				node = j;
			}
		}

		fevals += allnodes;

		/* Calculate coordinates of obtained upper bound */

		for (int k = dim - 1; k >= 0; k--) {
			int t = node % nodes;
			node = (int)(node / nodes);
			xfound[k] = a[k] + t * step[k];
		}

		/* Calculate all estimations of Lipshitz constant and choose maximum estimation */
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

		/* final calculation */

		LB = R * L * delta;
		*dL = LB;
		LB = Fr - LB;
		delete[]step; delete[]x;
		delete[]Fvalues;
		*Frp = Fr;
		*LBp = LB;
	}
};

/* Search on grid with OpenMP parallelization */

template <class T> class GridSolverOMP : public GridSolver <T> {
public:

	/**
	* Constructor
	*/
	GridSolverOMP() {
		np = omp_get_num_procs();
		omp_set_dynamic(0);
		omp_set_num_threads(np);
	}

protected:

	/* number of available processors */
	int np;

	virtual void GridEvaluator(int dim, const T *a, const T *b, T* xfound, T *Frp, T *LBp, T *dL, const std::function<T(const T * const)> &compute) {
		T Fr = std::numeric_limits<T>::max();
		T L = std::numeric_limits<T>::min();
		T delta = std::numeric_limits<T>::min();
		T LB;
		double R;
		T *step, *x, *Fvalues, *Frs, *Ls;
		x = nullptr;
		int *pts;
		try {
			step = new T[dim];
			Frs = new T[np];
			Ls = new T[np];
			pts = new int[np];
		}
		catch (std::bad_alloc& ba) {
			std::cerr << ba.what() << std::endl;
			this->errcode = -2;
			*Frp = Fr;
			*LBp = L;
			*dL = delta;
			return;
		}

		for (int k = 0; k < np; k++) {
			Frs[k] = std::numeric_limits<T>::max();
			Ls[k] = std::numeric_limits<T>::min();
		}

		for (int i = 0; i < dim; i++) {
			step[i] = fabs(b[i] - a[i]) / (this->nodes - 1);
			delta = step[i] > delta ? step[i] : delta;
		}

		delta = delta * 0.5 * dim;
		R = this->getR(delta);
		int allnodes = static_cast<int>(pow(this->nodes, dim));
		int node;
		try {
			Fvalues = new T[allnodes];
		}
		catch (std::bad_alloc& ba) {
			std::cerr << ba.what() << std::endl;
			this->errcode = -2;
			*Frp = Fr;
			*LBp = L;
			*dL = delta;
			return;
		}

#pragma omp parallel private(x) shared(dim, step, a, Fvalues)
		{
				x = new T[dim];
#pragma omp for
				for (int j = 0; j < allnodes; j++) {
					int point = j;
					for (int k = dim - 1; k >= 0; k--) {
						int t = point % this->nodes;
						point = (int)(point / this->nodes);
						x[k] = a[k] + t * step[k];
					}
					T rs = compute(x);
					Fvalues[j] = rs;
					int nt = omp_get_thread_num();
					if (rs < Frs[nt]) {
						Frs[nt] = rs;
						pts[nt] = j;
					}
				}
				delete[]x;
		}

		this->fevals += allnodes;

#pragma omp parallel for shared(dim,allnodes,Fvalues)
		for (int j = 0; j < allnodes; j++) {
			int neighbour;
			T loc = std::numeric_limits<T>::min();
			for (int k = 0; k < dim; k++) {
				int board = (int)pow(this->nodes, k + 1);
				neighbour = j + board / this->nodes;
				if ((neighbour < allnodes) && ((j / board) == (neighbour / board))) {
					loc = static_cast<T>(fabs(Fvalues[j] - Fvalues[neighbour])) / step[dim - 1 - k];
				}
				int nt = omp_get_thread_num();
				if (loc > Ls[nt]) {
					Ls[nt] = loc;
				}
			}
		}
		for (int k = 0; k < np; k++) {
			if (Frs[k] < Fr) {
				Fr = Frs[k];
				node = pts[k];
			}
			L = Ls[k] > L ? Ls[k] : L;
		}
		for (int k = dim - 1; k >= 0; k--) {
			int t = node % this->nodes;
			node = (int)(node / this->nodes);
			xfound[k] = a[k] + t * step[k];
		}
		/* final calculation */
		LB = R * L * delta;
		*dL = LB;
		LB = Fr - LB;
		delete[]step;
		delete[]Fvalues;
		delete[]pts; delete[]Frs; delete[]Ls;
		*Frp = Fr;
		*LBp = LB;
	}

};

/* hyperinterval info */
template <class T>
class part {
	int size;
public:
	T * a = nullptr;
	T *b = nullptr;
	T LocUP, LocLO, deltaL;
	part(int n, const T* a, const T* b) {
		try {
			size = n;
			this->a = new T[size];
			this->b = new T[size];
		}
		catch (std::bad_alloc e) {
			throw e;
			if (this->a) delete[]this->a;
			if (this->b) delete[]this->b;
			return;
		}
		for (int i = 0; i < size; i++) {
			this->a[i] = a[i];
			this->b[i] = b[i];
		}
	}
	part(const part& p) {
		if (this->a) delete[]a;
		if (this->b) delete[]b;
		this->size = p.size;
		try {
			this->a = new T[this->size];
			this->b = new T[this->size];
		}
		catch (std::bad_alloc e) {
			throw e;
			if (this->a) delete[]a;
			if (this->b) delete[]b;
			return;
		}
		for (int i = 0; i < this->size; i++) {
			this->a[i] = p.a[i];
			this->b[i] = p.b[i];
		}
		this->deltaL = p.deltaL;
		this->LocLO = p.LocLO;
		this->LocUP = p.LocUP;
	}
	part & operator=(const part & p) {
		if (this->a) delete[]a;
		if (this->b) delete[]b;
		this->size = p.size;
		try {
			this->a = new T[this->size];
			this->b = new T[this->size];
		}
		catch (std::bad_alloc e) {
			throw e;
			if (this->a) delete[]a;
			if (this->b) delete[]b;
			return *this;
		}
		for (int i = 0; i < this->size; i++) {
			this->a[i] = p.a[i];
			this->b[i] = p.b[i];
		}
		this->deltaL = p.deltaL;
		this->LocLO = p.LocLO;
		this->LocUP = p.LocUP;
		return *this;
	}
	void change(const T*a, const T*b) {
		for (int i = 0; i < this->size; i++) {
			this->a[i] = a[i];
			this->b[i] = b[i];
		}
	}
	~part() {
		if (a) delete[]a;
		if (b) delete[]b;
	}
};


#endif /* GRIDSOLVER_HPP */
