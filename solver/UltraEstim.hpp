#ifndef ULTRAESTIM_HPP
#define ULTRAESTIM_HPP

#include <omp.h>
#include <vector>
#include <functional>



template<class T> 
class L_accurate {
protected:
	T* grid = nullptr, *step = nullptr, *xp = nullptr, *Frs = nullptr, *Ls = nullptr;
	int dimension, all, nodes;
	bool fixed;
	int np, *pts = nullptr;
public:
	L_accurate() { //throws
		fixed = false;
		np = omp_get_num_procs();
		omp_set_dynamic(0);
		omp_set_num_threads(np);
		try {
			Frs = new T[np];
			Ls = new T[np];
			pts = new int[np];
		}
		catch (std::bad_alloc& e) {
			if (Frs != nullptr) delete[]Frs;
			if (Ls != nullptr) delete[]Ls;
			throw e;
		}
	}
	~L_accurate() {
		if (grid != nullptr) delete[]grid;
		if (xp != nullptr) delete[]xp;
		if (step != nullptr) delete[]step;
		if (Frs != nullptr) delete[]Frs;
		if (Ls != nullptr) delete[]Ls;
		if (pts != nullptr) delete[]pts;
	}
	virtual void stay_fixed(int d, int n) {
		if ((d > 0) && (n >= 2)) {
			fixed = true;
			dimension = d;
			nodes = n;
			all = static_cast<int>(pow(nodes, dimension));
			try {
				grid = new T[all];
				step = new T[dimension];
				xp = new T[np*dimension];
			}
			catch (std::bad_alloc &e) {
				if (grid != nullptr) delete[]grid;
				if (step != nullptr) delete[]step;
				fixed = false;
				throw e;
			}
		}
		else {
			fixed = false;
		}
	}
	virtual void unfix() {
		fixed = false;
		delete[]grid;
		delete[]step;
		delete[]xp;
		grid = nullptr;
		step = nullptr;
		xp = nullptr;
	}

	/*
		nodes - number of nodes in grid
		a - left bound of search region
		b - right bound of search region
		xfound - coordinate in which the global minimum found
		Frp - the global minimum
		compute - pointer to function, used for calculation
	*/
	/* simple function for calculation accurate estimation of L (implies a grid with a large number of nodes) */
	/* for case of one-dimensional problem, includes code-simplification and increased usability */
	T ultraoptimizer_1(const int nodes, const T a, const T b, T &xfound, T &Frp, const std::function<T(const T * const)> &compute) {
		T *Fvalues = nullptr;
		T step;

		try {
			Fvalues = new T[nodes];
			if (xp != nullptr) delete[]xp;
			xp = new T[np];
		}
		catch (std::bad_alloc& e) {
			if (Fvalues != nullptr) delete[]Fvalues;
			throw e;
		}

		T Fr = std::numeric_limits<T>::max();
		T L = std::numeric_limits<T>::min();
		int node;

		for (int k = 0; k < np; k++) {
			Frs[k] = std::numeric_limits<T>::max();
			Ls[k] = std::numeric_limits<T>::min();
		}

		step = fabs(b - a) / (nodes - 1);

#pragma omp parallel for
		for (int j = 0; j < nodes; j++) {
			int nt = omp_get_thread_num();
			int point = j;
			int t = point % nodes;
			point = (int)(point / nodes);
			xp[nt] = a + t * step;
			T rs = compute((xp + nt));
			Fvalues[j] = rs;
			if (rs < Frs[nt]) {
				Frs[nt] = rs;
				pts[nt] = j;
			}
		}


#pragma omp parallel for
		for (int j = 0; j < nodes; j++) {
			int neighbour;
			T loc = std::numeric_limits<T>::min();
			neighbour = j + 1;
			if (neighbour < nodes) {
				loc = static_cast<T>(fabs(Fvalues[j] - Fvalues[neighbour])) / step;
			}
			int nt = omp_get_thread_num();
			if (loc > Ls[nt]) {
				Ls[nt] = loc;
			}
		}
		for (int k = 0; k < np; k++) {
			if (Frs[k] < Fr) {
				Fr = Frs[k];
				node = pts[k];
			}
			L = Ls[k] > L ? Ls[k] : L;
		}
		xfound = a + node * step;

		Frp = Fr;
		if (Fvalues != nullptr) delete[]Fvalues;
		if (xp != nullptr) delete[]xp;
		xp = nullptr;
		return L;
	}

	/* Standart function for accurate estimating L with possible large grid
	   In case of multidimensional problem */
	/*
		Additional parameters:
		dim - dimensionality of a problem
		Also because the problem multidimensional there are arrays
		for left and right bound, such for coordinates of global minimum too
	*/
	T ultraoptimizer_m(const int dim, const int nodes, const T *a, const T *b, T* xfound, T &Frp, const std::function<T(const T * const)> &compute) {
		T *Fvalues = nullptr, *step = nullptr;

		T Fr = std::numeric_limits<T>::max();
		T L = std::numeric_limits<T>::min();
		T delta = 0.0;
		int node;
		int allnodes = static_cast<int>(pow(nodes, dim));

		try {
			step = new T[dim];
			if (xp != nullptr) delete[]xp;
			xp = new T[np*dim];
			Fvalues = new T[allnodes];
		}
		catch (std::bad_alloc& e) {
			if (step != nullptr) delete[]step;
			if (xp != nullptr) delete[]xp;
			throw e;
		}
		for (int k = 0; k < np; k++) {
			Frs[k] = std::numeric_limits<T>::max();
			Ls[k] = std::numeric_limits<T>::min();
		}

		for (int i = 0; i < dim; i++) {
			step[i] = fabs(b[i] - a[i]) / (nodes - 1);
			delta += step[i] / 2.0;
		}


#pragma omp parallel for
		for (int j = 0; j < allnodes; j++) {
			int nt = omp_get_thread_num();
			int point = j;
			for (int k = dim - 1; k >= 0; k--) {
				int t = point % nodes;
				point = (int)(point / nodes);
				xp[dim*nt + k] = a[k] + t * step[k];
			}
			T rs = compute((xp + dim * nt));
			Fvalues[j] = rs;
			if (rs < Frs[nt]) {
				Frs[nt] = rs;
				pts[nt] = j;
			}
		}

#pragma omp parallel for
		for (int j = 0; j < allnodes; j++) {
			int neighbour;
			T loc = std::numeric_limits<T>::min();
			for (int k = 0; k < dim; k++) {
				int board = static_cast<int>(pow(nodes, k + 1));
				neighbour = j + board / nodes;
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
			int t = node % nodes;
			node = static_cast<int>(node / nodes);
			xfound[k] = a[k] + t * step[k];
		}

		Frp = Fr;
		if (Fvalues != nullptr) delete[]Fvalues;
		if (xp != nullptr) delete[]xp;
		xp = nullptr;
		if (step != nullptr) delete[]step;
		return L;

	}

	T ultraoptimizer_f(const T *a, const T *b, T* xfound, T &Frp, const std::function<T(const T * const)> &compute) {
		T Fr = std::numeric_limits<T>::max();
		T L = std::numeric_limits<T>::min();
		T delta = static_cast<T>(0);
		if (this->fixed == false) {
			return L;
		}
		int node;

		for (int k = 0; k < np; k++) {
			Frs[k] = std::numeric_limits<T>::max();
			Ls[k] = std::numeric_limits<T>::min();
		}

		for (int i = 0; i < dimension; i++) {
			step[i] = fabs(b[i] - a[i]) / (nodes - 1);
			delta += step[i] / 2.0;
		}


#pragma omp parallel for
		for (int j = 0; j < all; j++) {
			int nt = omp_get_thread_num();
			int point = j;
			for (int k = dimension - 1; k >= 0; k--) {
				int t = point % nodes;
				point = (int)(point / nodes);
				xp[dimension*nt + k] = a[k] + t * step[k];
			}
			T rs = compute((xp + dimension * nt));
			grid[j] = rs;
			if (rs < Frs[nt]) {
				Frs[nt] = rs;
				pts[nt] = j;
			}
		}

#pragma omp parallel for
		for (int j = 0; j < all; j++) {
			int neighbour;
			T loc = std::numeric_limits<T>::min();
			for (int k = 0; k < dimension; k++) {
				int board = static_cast<int>(pow(nodes, k + 1));
				neighbour = j + board / nodes;
				if ((neighbour < all) && ((j / board) == (neighbour / board))) {
					loc = static_cast<T>(fabs(grid[j] - grid[neighbour])) / step[dimension - 1 - k];
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
		for (int k = dimension - 1; k >= 0; k--) {
			int t = node % this->nodes;
			node = static_cast<int>(node / nodes);
			xfound[k] = a[k] + t * step[k];
		}

		Frp = Fr;
		return L;

	}
};

#endif