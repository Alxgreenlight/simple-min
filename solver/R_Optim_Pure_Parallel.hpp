#ifndef R_OPTIM_HPP
#define R_OPTIM_HPP

#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <omp.h>
#include "../common/bbsolver.hpp"
#include <fstream>
#include <chrono>


template <class T>
class Box;

using us = std::chrono::duration_values<std::chrono::microseconds>;

template <class T>
class PPrOptimizer : public BlackBoxSolver <T> {
protected:
	int dim, np; //dimension and number of threads
	int all; //number of nodes in mesh
	T eps; //accuracy
	T UpperBound;
	T maxL; //maximum used value of Lipschitz constant per search
	unsigned long long FuncEvals; //number of evaluations per search
	unsigned long Iterations; //number of BnB iterations per search
	bool initialized; //some parameters must be provided before search
	const int l = 5, m = 9, h = 17; //number of nodes per dimension in 3 meshes
	T *grid = nullptr; //array of function values in nodes of mesh
	T *x = nullptr; //contains coordinates in which function value calculating
	T *step = nullptr; //step of mesh in every dimension
	T *ests = nullptr; //for parallel estimations of L calculating
	int *coords = nullptr; //need for nested iteration
	/*temporary*/
		//for time measurement
		std::chrono::steady_clock sc;
		std::chrono::microseconds gridtime = us::zero(), esttime = us::zero(), est1 = us::zero(), est2 = us::zero(), est3 = us::zero(), minsearch = us::zero();
	/*temporary*/
public:
	PPrOptimizer() { //throws
		if (!std::numeric_limits<T>::is_specialized) {
			throw std::runtime_error("Bad Type");
		}
		initialized = false;
		np = omp_get_num_procs();
		omp_set_dynamic(0);
		omp_set_num_threads(np);
	}

	~PPrOptimizer() {
		if (x != nullptr) delete[]x;
		if (step != nullptr) delete[]step;
		if (grid != nullptr) delete[]grid;
		if (coords != nullptr) delete[]coords;
		if (ests != nullptr) delete[]ests;
	}

	virtual void clear() {
		if (x != nullptr) delete[]x;
		if (step != nullptr) delete[]step;
		if (grid != nullptr) delete[]grid;
		if (coords != nullptr) delete[]coords;
		if (ests != nullptr) delete[]ests;
		x = nullptr;
		step = nullptr;
		grid = nullptr;
		coords = nullptr;
		ests = nullptr;
		initialized = false;
	}

	//set required parameters (dimension of task and required accuracy)
	virtual void init(const int d, const T e) {  //throws
		clear(); 
		if ((d >= 1) && (e > std::numeric_limits<T>::min())) {
			dim = d;
			eps = e;
			all = static_cast<int>(pow(h, dim));
		}
		else {
			throw std::runtime_error("Bad initialization");
		}
		try {
			grid = new T[all];
			x = new T[dim * np];
			step = new T[dim];
			coords = new int[np*(dim + 1)];
			ests = new T[np];
		}
		catch (std::bad_alloc& e) {
			throw e;
		}
		initialized = true;
	}

	virtual void getInfo(unsigned long long int &fevs, unsigned long int &iters, T &Lmax) {
		fevs = FuncEvals;
		iters = Iterations;
		Lmax = maxL;
	}

	virtual void debug_time(std::ofstream & fp) {
		fp << "Timing here:" << std::endl;
		fp << "Time in grid evaluation: " << gridtime.count() << " us" << std::endl;
		fp << "Time in estimations: " << esttime.count() << " us" << std::endl;
			/*
			fp << "Time in estimation 1: " << est1.count() << " us" << std::endl;
			fp << "Time in estimation 2: " << est2.count() << " us" << std::endl;
			fp << "Time in estimation 3: " << est3.count() << " us" << std::endl;
			*/
		fp << "Time in searching minimum: " << minsearch.count() << " us" << std::endl;
	}

	virtual T search(int n, T* xfound, const T * const a, const T * const b, const std::function<T(const T * const)> &f) {  //throws
		if ((!initialized) || (n != dim)) {
			throw std::runtime_error("Provide correct initialization first (from search func.)");
		}
		if (f == nullptr) {
			throw std::runtime_error("Pointer to calculation function is incorrect");
		}

		FuncEvals = 0;
		Iterations = 0;
		gridtime = us::zero(); esttime = us::zero(); est1 = us::zero(); est2 = us::zero(); est3 = us::zero(); minsearch = us::zero();
		UpperBound = std::numeric_limits<T>::max();
		maxL = std::numeric_limits<T>::min();

		std::vector<Box<T>> curBox, nextBox;
		try {
			curBox.emplace_back(static_cast<unsigned short>(dim), a, b, false);
		}
		catch (std::exception& e) {
			throw e;
		}

		/* If hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates */
		/* xs temporary array for storage local min coordinates */
		T *a1, *b1, *xs;
		try {
			a1 = new T[dim];
			b1 = new T[dim];
			xs = new T[dim];
		}
		catch (std::bad_alloc& ba) {
			throw ba;
		}

		/* Each hyperinterval can be subdivided or pruned (if non-promisable or fits accuracy) */
		while (!curBox.empty()) {
			unsigned int boxes = curBox.size(); /* number of iterations on this step */
			Iterations += boxes;
			
			/* For all hyperintervals on this step perform grid search */
			for (unsigned int i = 0; i < boxes; i++) {
				T lUPB;
				BoxOptimizator(curBox[i], xs, lUPB, f);
				/* remember new results if less then previous */
				UpdateRecords(lUPB, xfound, xs);
			}
			
			/* Choose which hyperintervals should be subdivided */
			for (unsigned int i = 0; i < boxes; i++) {
				/* Subdivision criteria */
				if ((!curBox[i].ready) || ((curBox[i].ready) &&  (curBox[i].LocLO < (UpperBound - eps)))) {
					/* If subdivide, choose dimension (the longest side) */
					
					int choosen = ChooseDim(curBox[i].a, curBox[i].b);
					/* Make new edges for 2 new hyperintervals */
					for (int j = 0; j < dim; j++) {
						if (j != choosen) {			/* [a .. b1] [a1 .. b] */
							a1[j] = curBox[i].a[j];		/* where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ] */
							b1[j] = curBox[i].b[j];		/* and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ] */
						}
						else {
							a1[j] =curBox[i].a[j] + fabs(curBox[i].b[j] - curBox[i].a[j]) / 2.0;
							b1[j] = a1[j];
						}
					}
					/* Add 2 new hyperintervals, parent HI no longer considered */
					try {
						nextBox.emplace_back(static_cast<unsigned short>(dim), curBox[i].a, b1, curBox[i].ready, curBox[i].L_estim, curBox[i].attempts);
						nextBox.emplace_back(static_cast<unsigned short>(dim), a1, curBox[i].b, curBox[i].ready, curBox[i].L_estim, curBox[i].attempts);
					}
					catch (std::exception& e) {
						throw e;
					}
				}
			}

			curBox.clear();
			curBox.swap(nextBox);
		}
		
		delete[]a1; delete[]b1; delete[]xs;
		curBox.clear();
		return UpperBound;
	}

protected:
	virtual void UpdateRecords(const T LU, T* x, const T *xs) {
		if (LU < UpperBound) {
			UpperBound = LU;
			for (int i = 0; i < dim; i++) {
				x[i] = xs[i];
			}
		}
	}

	virtual int ChooseDim(const T *a, const T *b) {
		T max = std::numeric_limits<T>::min(), diff;
		int i, maxI = 0;
		for (i = 0; i < dim; i++) {
			diff = fabs(b[i] - a[i]);
			if (diff > max) {
				max = diff;
				maxI = i;
			}
		}
		return maxI;
	}

	virtual int getIndex(int inc, int *coords, const int &nod) {
		int i;
		for (i = 0; i <= dim; i++) {
			if (coords[i] >= nod) {
				coords[i + 1] += static_cast<int>(coords[i] / nod);
				coords[i] = coords[i] % nod;
			}
			else {
				break;
			}
		}
		int j = 0;
		for (i = 0; i <= dim; i++) {
			j += coords[i] * inc* static_cast<int>(pow(h, i));
		}
		return j;
	}


	virtual T estimation(const int nod) {
		int inc = 1;
		T est = std::numeric_limits<T>::min();
		for (int i = 0; i < np; i++) {
			ests[i] = std::numeric_limits<T>::min();
		}
		if (nod == l) {
			inc = 4;
		}
		else if (nod == m) {
			inc = 2;
		}
		else if (nod == h) {
			inc = 1;
		}
		if (inc != 1) {
			for (int i = 0; i < np*(dim + 1); i++) {
				coords[i] = 0;
			}
		}
#pragma omp parallel
		{
			int nt = omp_get_thread_num();
			int j = nt;
			if (inc != 1) {
				coords[nt*(dim + 1)] = nt;
				j = getIndex(inc, coords + nt * (dim + 1), nod);
			}
			while (j < all) {
				int neighbour;
				if (inc != 1) {
					j = getIndex(inc, coords + nt * (dim + 1), nod);
					coords[nt*(dim + 1)] += np;
				}
				for (int k = 0; k < dim; k++) {
					int board = static_cast<int>(pow(h, k + 1));
					neighbour = j + inc * board / h;
					if ((neighbour < all) && ((j / board) == (neighbour / board))) {
						T loc = fabs(grid[j] - grid[neighbour]) / (inc * step[dim - 1 - k]);
						ests[nt] = loc > ests[nt] ? loc : ests[nt];
					}
				}
				if (inc == 1) {
					j += np;
				}
			}
		}
		for (int i = 0; i < np; i++) {
			est = ests[i] > est ? ests[i] : est;
		}
		return est;
	}

	virtual T clarification(const int param, const T le, const T me, const T he) {
		T delta = 0.0; 
		switch (param) {
		case 1:
			for (int i = 0; i < dim; i++) {
				delta += step[i] / 2.0;
			}
			return he * exp(delta);
			break;
		case 2:
			return he + fabs(he - le) / 2.0;
			break;
		default:
			return he;
			break;
		}
	}

	virtual bool isLstable(Box<T>& B, const std::function<T(const T * const)> &compute) {
		for (int i = 0; i < dim; i++) {
			step[i] = fabs(B.b[i] - B.a[i]) / (h - 1);
		}
		auto s = sc.now();
#pragma omp parallel for
		for (int i = 0; i < all; i++) {
			int nt = omp_get_thread_num();
			int point = i;
			for (int k = dim - 1; k >= 0; k--) {
				int t = point % h;
				point = point / h;
				x[nt*dim + k] = B.a[k] + t * step[k];
			}
			grid[i] = compute(x + nt * dim);
		}
		auto e = sc.now();
		gridtime += std::chrono::duration_cast<std::chrono::microseconds>(e - s);
		FuncEvals += all;
		T lEst, mEst, hEst;
		s = sc.now();
			lEst = estimation(l);
			mEst = estimation(m);
			hEst = estimation(h);
		e = sc.now();
		esttime += std::chrono::duration_cast<std::chrono::microseconds>(e - s);
		if ((fabs(mEst - lEst) > fabs(hEst - mEst)) && (fabs(hEst - mEst) < 0.05*hEst) && (fabs(mEst - lEst) < 0.1*mEst)) {
			T Local_estim = clarification(2, lEst, mEst, hEst);
			if ((B.ready == false) || (fabs(B.L_estim - Local_estim) > 0.1*[](T f, T s) { return f > s ? f : s; }(B.L_estim, Local_estim))) {
				B.ready = true;
				B.attempts = 0;
				B.L_estim = Local_estim;
				if (B.L_estim > maxL) {
					maxL = B.L_estim;
				}
			}
			return true;
		}
		else {
			B.ready = false;
			return false;
		}
	}

	virtual void BoxOptimizator(Box<T>& B, T* xfound, T& UpperBound, const std::function<T(const T * const)> &compute) {
		UpperBound = std::numeric_limits<T>::max();
		if (!isLstable(B, compute)) {
			B.attempts++;
			return;
		}
		int node = 0;
		T delta = 0.0;
		T Fr = std::numeric_limits<T>::max();
		for (int i = 0; i < dim; i++) {
			delta += step[i] / 2.0;
		}

		auto s = sc.now();
		for (int j = 0; j < all; j++) {
			if (grid[j] < Fr) {
				Fr = grid[j];
				node = j;
			}
		}
		auto e = sc.now();
		minsearch += std::chrono::duration_cast<std::chrono::microseconds>(e - s);

		for (int k = dim - 1; k >= 0; k--) {
			int t = node % h;
			node = node / h;
			xfound[k] = B.a[k] + t * step[k];
		}

		UpperBound = Fr;
		B.LocUP = Fr;
		B.LocLO = Fr - B.L_estim * delta;
	}
};

/* hyperinterval info */
template <class T>
class Box {
	unsigned short size; //dimension of task
public:
	T *a = nullptr; //left bounds of hyperinterval
	T *b = nullptr; //right bounds of hyperinterval
	T LocUP, LocLO; //upper and lower estimates of minimum on the hyperinterval
	T L_estim; //estimation for Lipshitz const (usable only if (ready))
	bool ready = false; //Lipshitz constant estimated reliably
	unsigned int attempts = 0; //how many times L had tried to be estimated, but unsuccessful
	Box(const unsigned short n, const T* a, const T* b, const bool r, const double L = 0.0, const unsigned short at = 0) {
		try {
			this->size = n;
			this->a = new T[size];
			this->b = new T[size];
			this->ready = r;
			this->L_estim = L;
			this->attempts = at;
		}
		catch (const std::bad_alloc& e) {
			if (this->a) delete[]this->a;
			if (this->b) delete[]this->b;
			throw e;
			return;
		}
		for (int i = 0; i < size; i++) {
			this->a[i] = a[i];
			this->b[i] = b[i];
		}
	}
	Box(const Box& B) {
		if (this->a) delete[]a;
		if (this->b) delete[]b;
		this->size = B.size;
		this->ready = B.ready;
		try {
			this->a = new T[this->size];
			this->b = new T[this->size];
		}
		catch (const std::bad_alloc& e) {
			if (this->a) delete[]this->a;
			if (this->b) delete[]this->b;
			throw e;
			return;
		}
		for (int i = 0; i < this->size; i++) {
			this->a[i] = B.a[i];
			this->b[i] = B.b[i];
		}
		this->L_estim = B.L_estim;
		this->attempts = B.attempts;
		this->LocLO = B.LocLO;
		this->LocUP = B.LocUP;
	}
	Box & operator=(const Box & B) {
		if (this->a) delete[]this->a;
		if (this->b) delete[]this->b;
		this->size = B.size;
		try {
			this->a = new T[this->size];
			this->b = new T[this->size];
		}
		catch (const std::bad_alloc& e) {
			if (this->a) delete[]a;
			if (this->b) delete[]b;
			throw e;
			return *this;
		}
		for (int i = 0; i < this->size; i++) {
			this->a[i] = B.a[i];
			this->b[i] = B.b[i];
		}  
		this->L_estim = B.L_estim;
		this->attempts = B.attempts;
		this->LocLO = B.LocLO;
		this->LocUP = B.LocUP;
		return *this;
	}
	~Box() {
		if (this->a) delete[]this->a;
		if (this->b) delete[]this->b;
	}
};
#endif
