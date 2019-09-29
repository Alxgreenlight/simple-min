#ifndef ULTRAONEDIM_HPP
#define ULTRAONEDIM_HPP

#include <omp.h>
#include <functional>

template<class T>
T ultraoptimizer(const int nodes, const T *a, const T *b, T* xfound, T *Frp, const std::function<T(const T * const)> &compute) {
	int np, xpp;
	T *Fvalues = nullptr; 
	T *Frs = nullptr, *Ls = nullptr, *xp = nullptr;
	int *pts = nullptr;
	T step;

	np = omp_get_num_procs();
	omp_set_dynamic(0);
	omp_set_num_threads(np);
	Fvalues = new T[nodes];
	Frs = new T[np];
	Ls = new T[np];
	pts = new int[np];
	xp = new T[np];
	xpp = 1;

	T Fr = std::numeric_limits<T>::max();
	T L = std::numeric_limits<T>::min();
	int node;

	for (int k = 0; k < np; k++) {
		Frs[k] = std::numeric_limits<T>::max();
		Ls[k] = std::numeric_limits<T>::min();
	}

	step = fabs(b[0] - a[0]) / (nodes - 1);

#pragma omp parallel for
	for (int j = 0; j < nodes; j++) {
		int nt = omp_get_thread_num();
		int point = j;
		int t = point % nodes;
		point = (int)(point / nodes);
		xp[nt] = a[0] + t * step;
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
	xfound[0] = a[0] + node * step;
	/* final calculation */
	*Frp = Fr;
	if (Fvalues != nullptr) delete[]Fvalues;
	if (Frs != nullptr) delete[]Frs;
	if (Ls != nullptr) delete[]Ls;
	if (pts != nullptr) delete[]pts;
	if (xp != nullptr) delete[]xp;
	return L;
}

#endif