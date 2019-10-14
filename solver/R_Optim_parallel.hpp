#ifndef R_OPTIM_PARALLEL_HPP
#define R_OPTIM_PARALLEL_HPP

#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <omp.h>
#include "R_Optim.hpp"
/*for debugging*/ #include <iostream>


template <class T>
class PrOptimizer_v1 : public rOptimizer <T> {
protected:
	int np;
public:
	PrOptimizer_v1() try : rOptimizer<T>() {//throws
		np = omp_get_num_procs();
		omp_set_dynamic(0);
		omp_set_num_threads(np);
	}
	catch (std::exception& e) {
		throw e;
	}

	~PrOptimizer_v1() {
	}

	virtual void init(const int d, const T e) {  //throws
		if ((d >= 1) && (e > std::numeric_limits<T>::min())) {
			this->dim = d;
			this->eps = e;
			this->all = static_cast<int>(pow(this->h, this->dim));
			this->initialized = true;
		}
		else {
			this->initialized = false;
			throw std::runtime_error("Bad initialization");
		}
		try {
			this->grid = new T[static_cast<int>(pow(this->h, this->dim))];
			this->x = new T[this->dim * np];
			this->step = new T[this->dim];
			this->coords = new int[this->dim + 1];
		}
		catch (std::bad_alloc& e) {
			this->initialized = false;
			throw e;
		}
	}


protected:

	virtual int incr(int inc) {
		int i;
		for (i = 0; i <= this->dim; i++) {
			this->coords[i] += inc;
			if (this->coords[i] > this->h) {
				this->coords[i] = 0;
			}
			else {
				break;
			}
		}
		int j = 0;
		for (i = 0; i <= this->dim; i++) {
			j += this->coords[i] * static_cast<int>(pow(this->h, i));
		}
		return j;
	}


	virtual T estimation(const int nod) {
		int inc;
		T est = std::numeric_limits<T>::min();
		if (nod == this->l) {
			inc = 4;
		}
		else if (nod == this->m) {
			inc = 2;
		}
		else if (nod == this->h) {
			inc = 1;
		}
		if (inc != 1) {
			for (int i = 0; i <= this->dim; i++) {
				this->coords[i] = 0;
			}
		}
		int j = 0;
		while (j < this->all) {
			int neighbour;
			for (int k = 0; k < this->dim; k++) {
				int board = static_cast<int>(pow(this->h, k + 1));
				neighbour = j + inc * board / this->h;
				if ((neighbour < this->all) && ((j / board) == (neighbour / board))) {
					T loc = fabs(this->grid[j] - this->grid[neighbour]) / (inc * this->step[this->dim - 1 - k]);
					est = loc > est ? loc : est;
				}
			}
			if (inc == 1) {
				j++;
			}
			else {
				j = incr(inc);
			}
		}
		return est;
	}

	virtual bool isLstable(Box<T>& B, const std::function<T(const T * const)> &compute) {
		for (int i = 0; i < this->dim; i++) {
			this->step[i] = fabs(B.b[i] - B.a[i]) / (this->h - 1);
		}
#pragma omp parallel for
		for (int i = 0; i < this->all; i++) {
			int nt = omp_get_thread_num();
			int point = i;
			for (int k = this->dim - 1; k >= 0; k--) {
				int t = point % this->h;
				point = point / this->h;
				this->x[nt*this->dim + k] = B.a[k] + t * this->step[k];
			}
			this->grid[i] = compute(this->x + nt * this->dim);
		}
		this->FuncEvals += this->all;
		T lEst, mEst, hEst;
#pragma omp sections
		{
#pragma omp section
			{
				lEst = estimation(this->l);
				mEst = estimation(this->m);
			}
#pragma omp section
			{
				hEst = estimation(this->h);
			}
		}
		if ((fabs(mEst - lEst) > fabs(hEst - mEst)) && (fabs(hEst - mEst) < 0.05*hEst) && (fabs(mEst - lEst) < 0.1*mEst)) {
			T Local_estim = this->clarification(2, lEst, mEst, hEst);
			if ((B.ready == false) || (fabs(B.L_estim - Local_estim) > 0.1*this->mymax(B.L_estim, Local_estim))) {
				B.ready = true;
				B.L_estim = Local_estim;
				if (B.L_estim > this->maxL) {
					this->maxL = B.L_estim;
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
		int node;
		T delta = 0.0;
		T Fr = std::numeric_limits<T>::max();
		for (int i = 0; i < this->dim; i++) {
			delta += this->step[i] / 2.0;
		}
		for (int j = 0; j < this->all; j++) {
			if (this->grid[j] < Fr) {
				Fr = this->grid[j];
				node = j;
			}
		}



		for (int k = this->dim - 1; k >= 0; k--) {
			int t = node % this->h;
			node = node / this->h;
			xfound[k] = B.a[k] + t * this->step[k];
		}

		UpperBound = Fr;
		B.LocUP = Fr;
		B.LocLO = Fr - B.L_estim * delta;
	}
};

template <class T>
class PrOptimizer_v2 : public rOptimizer <T> {
protected:
	int np;
	T* UPBs, *XFs;
public:
	PrOptimizer_v2() try : rOptimizer<T>() { //throws
		np = omp_get_num_procs();
		omp_set_dynamic(0);
		omp_set_num_threads(np);
	}
	catch (std::exception &e) {
		throw e;
	}

	~PrOptimizer_v2() {
		if (UPBs != nullptr) delete[]UPBs;
		if (XFs != nullptr) delete[]XFs;
	}

	virtual void clear() {
		if (this->x != nullptr) delete[]this->x;
		if (this->step != nullptr) delete[]this->step;
		if (this->grid != nullptr) delete[]this->grid;
		if (this->coords != nullptr) delete[]this->coords;
		if (UPBs != nullptr) delete[]UPBs;
		if (XFs != nullptr) delete[]XFs;
		this->x = nullptr;
		this->step = nullptr;
		this->grid = nullptr;
		this->coords = nullptr;
		UPBs = nullptr;
		XFs = nullptr;
		this->initialized = false;
	}

	virtual void init(const int d, const T e) {  //throws
		if ((d >= 1) && (e > std::numeric_limits<T>::min())) {
			this->dim = d;
			this->eps = e;
			this->all = static_cast<int>(pow(this->h, this->dim));
			this->initialized = true;
		}
		else {
			this->initialized = false;
			throw std::runtime_error("Bad initialization");
		}
		try {
			this->grid = new T[np*this->all];
			this->x = new T[np*this->dim];
			this->step = new T[np*this->dim];
			this->coords = new int[np*(this->dim + 1)];
			UPBs = new T[np];
			XFs = new T[this->dim*np];
		}
		catch (std::bad_alloc& e) {
			this->initialized = false;
			throw e;
		}
	}

	virtual T search(int n, T* xfound, const T * const a, const T * const b, const std::function<T(const T * const)> &f) {  //throws
		if ((!this->initialized) || (n != this->dim)) {
			throw std::runtime_error("Provide correct initialization first (from search func.)");
		}
		this->FuncEvals = 0;
		this->Iterations = 0;
		std::vector<Box<T>> curBox, nextBox;
		this->UpperBound = std::numeric_limits<T>::max();
		for (int i = 0; i < np; i++) {
			UPBs[i] = std::numeric_limits<T>::max();
		}
		if (f == nullptr) {
			throw std::runtime_error("Pointer to calculation function is incorrect");
		}
		try {
			curBox.emplace_back(static_cast<unsigned short>(this->dim), a, b, false);
		}
		catch (std::exception& e) {
			throw e;
		}

		/* If hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates */
		/* xs temporary array for storage local min coordinates */
		T *a1, *b1, *xs;
		try {
			a1 = new T[this->dim];
			b1 = new T[this->dim];
			xs = new T[np*this->dim];
		}
		catch (std::bad_alloc& ba) {
			throw ba;
		}

		/* Each hyperinterval can be subdivided or pruned (if non-promisable or fits accuracy) */
		while (!curBox.empty()) {
			/* number of iterations on this step (BFS) */
			int boxes = curBox.size();
			this->Iterations += boxes;

			/* For all hyperintervals on this step perform grid search */
#pragma omp parallel for
			for (int i = 0; i < boxes; i++) {
				T lUPB;
				int nt = omp_get_thread_num();
				try {
					this->BoxOptimizator(curBox[i], xs, lUPB, f, nt);
				}
				catch (std::bad_alloc &ba) {
					throw ba;
				}
				/* remember new results if less then previous */
				UpdateRecords(lUPB, nt, XFs, xs);
			}

			FinalUpdate(xfound);

			/* Choose which hyperintervals should be subdivided */
			for (int i = 0; i < boxes; i++) {
				/* Subdivision criteria */
				if ((!curBox[i].ready) || ((curBox[i].ready) && (curBox[i].LocLO < (this->UpperBound - this->eps)))) {
					/* If subdivide, choose dimension (the longest side) */

					int choosen = this->ChooseDim(curBox[i].a, curBox[i].b);
					/* Make new edges for 2 new hyperintervals */
					for (int j = 0; j < this->dim; j++) {
						if (j != choosen) {			/* [a .. b1] [a1 .. b] */
							a1[j] = curBox[i].a[j];		/* where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ] */
							b1[j] = curBox[i].b[j];		/* and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ] */
						}
						else {
							a1[j] = curBox[i].a[j] + fabs(curBox[i].b[j] - curBox[i].a[j]) / 2.0;
							b1[j] = a1[j];
						}
					}
					/* Add 2 new hyperintervals, parent HI no longer considered */
					try {
						nextBox.emplace_back(static_cast<unsigned short>(this->dim), curBox[i].a, b1, curBox[i].ready, curBox[i].L_estim, curBox[i].attempts);
						nextBox.emplace_back(static_cast<unsigned short>(this->dim), a1, curBox[i].b, curBox[i].ready, curBox[i].L_estim, curBox[i].attempts);
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
		return this->UpperBound;
	}

	void printx() {
		for (int i = 0; i < np; i++) {
			for (int j = 0; j < this->dim; j++) {
				std::cout << XFs[i*this->dim + j] << std::endl;
			}
		}
	}

protected:
	virtual void FinalUpdate(T *x_dest) {
		for (int i = 0; i < np; i++) {
			if (UPBs[i] < this->UpperBound) {
				this->UpperBound = UPBs[i];
				for (int j = 0; j < this->dim; j++) {
					x_dest[j] = XFs[i*this->dim + j];
				}
			}
		}
	}

	virtual void UpdateRecords(const T LU, const int nt, T* x_dest, const T *x_source) {
		if (LU < UPBs[nt]) {
			UPBs[nt] = LU;
			for (int i = 0; i < this->dim; i++) {
				x_dest[nt*this->dim + i] = x_source[nt*this->dim + i];
			}
		}
	}

	virtual int increment(int inc, const int nt) {
		int i;
		for (i = 0; i <= this->dim; i++) {
			this->coords[(this->dim+1)*nt + i] += inc;
			if (this->coords[(this->dim+1)*nt + i] > this->h) {
				this->coords[(this->dim+1)*nt + i] = 0;
			}
			else {
				break;
			}
		}
		int j = 0;
		for (i = 0; i <= this->dim; i++) {
			j += this->coords[(this->dim+1)*nt + i] * static_cast<int>(pow(this->h, i));
		}
		return j;
	}


	virtual T estimation(const int nod, const int nt) {
		int inc;
		T est = std::numeric_limits<T>::min();
		if (nod == this->l) {
			inc = 4;
		}
		else if (nod == this->m) {
			inc = 2;
		}
		else if (nod == this->h) {
			inc = 1;
		}
		if (inc != 1) {
			for (int i = 0; i <= this->dim; i++) {
				this->coords[nt*(this->dim + 1) + i] = 0;
			}
		}
		int j = 0;
		while (j < this->all) {
			int neighbour;
			for (int k = 0; k < this->dim; k++) {
				int board = static_cast<int>(pow(this->h, k + 1));
				neighbour = j + inc * board / this->h;
				if ((neighbour < this->all) && ((j / board) == (neighbour / board))) {
					T loc = fabs(this->grid[this->all*nt + j] - this->grid[this->all*nt + neighbour]) / (inc * this->step[this->dim*nt + (this->dim - 1 - k)]);
					est = loc > est ? loc : est;
				}
			}
			if (inc == 1) {
				j++;
			}
			else {
				j = this->increment(inc, nt);
			}
		}
		return est;
	}

	virtual T clarification(const int nt, const int param, const T le, const T me, const T he) {
		T delta = 0.0;
		switch (param) {
		case 1:
			for (int i = 0; i < this->dim; i++) {
				delta += this->step[nt*this->dim + i] / 2.0;
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


	virtual bool isLstable(Box<T>& B, const std::function<T(const T * const)> &compute, const int nt) {
		for (int i = 0; i < this->dim; i++) {
			this->step[this->dim*nt + i] = fabs(B.b[i] - B.a[i]) / (this->h - 1);
		}
		for (int i = 0; i < this->all; i++) {
			int point = i;
			for (int k = this->dim - 1; k >= 0; k--) {
				int t = point % this->h;
				point = point / this->h;
				this->x[nt*this->dim + k] = B.a[k] + t * this->step[nt*this->dim + k];
			}
			this->grid[this->all*nt + i] = compute(this->x + nt * this->dim);
		}
//#pragma omp atomic
		this->FuncEvals += this->all;
		T lEst, mEst, hEst;
		lEst = this->estimation(this->l, nt);
		mEst = this->estimation(this->m, nt);
		hEst = this->estimation(this->h, nt);
		if ((fabs(mEst - lEst) > fabs(hEst - mEst)) && (fabs(hEst - mEst) < 0.05*hEst) && (fabs(mEst - lEst) < 0.1*mEst)) {
			T Local_estim = this->clarification(nt, 2, lEst, mEst, hEst);
			if ((B.ready == false) || (fabs(B.L_estim - Local_estim) > 0.1*this->mymax(B.L_estim, Local_estim))) {
				B.ready = true;
				B.L_estim = Local_estim;
				if (B.L_estim > this->maxL) {
					this->maxL = B.L_estim;
				}
			}
			return true;
		}
		else {
			B.ready = false;
			return false;
		}
	}

	virtual void BoxOptimizator(Box<T>& B, T* xfound, T& UpperBound, const std::function<T(const T * const)> &compute, const int nt) {
		UpperBound = std::numeric_limits<T>::max();
		if (!isLstable(B, compute, nt)) {
			B.attempts++;
			return;
		}
		int node;
		T delta = 0.0;
		T Fr = std::numeric_limits<T>::max();
		for (int i = 0; i < this->dim; i++) {
			delta += this->step[nt*this->dim + i] / 2.0;
		}
		for (int j = 0; j < this->all; j++) {
			if (this->grid[this->all*nt + j] < Fr) {
				Fr = this->grid[this->all*nt + j];
				node = j;
			}
		}

		for (int k = this->dim - 1; k >= 0; k--) {
			int t = node % this->h;
			node = node / this->h;
			xfound[nt*this->dim + k] = B.a[k] + t * this->step[nt*this->dim + k];
		}

		UpperBound = Fr;
		B.LocUP = Fr;
		B.LocLO = Fr - B.L_estim * delta;
	}
};













//virtual void set_parallel() {
//	int s = 0;
//	hp = 2; lp = 2;
//	int edge = this->dim * 4;
//	while ((hp*lp < np) && (hp*lp < edge)) {
//		if (s < 2) {
//			hp++;
//		}
//		else {
//			lp++;
//		}
//		if (s == 3) {
//			s = 0;
//		}
//		else {
//			s++;
//		}
//	}
//}


#endif
