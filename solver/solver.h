#ifndef SOLVER_H
#define SOLVER_H

namespace solver {
	extern int nodes, dim;
	extern double eps;
	extern double UPB, LOB, deltaL, glob;
	extern double(*compute)(double *x);
	void checkErrors(int code); //check errors after call min_search and print error info
	int min_search(double*a, double *b);
}
/*find upper and lower bound for multidimensional [a,b]
(pointer *compute, amount of nodes per dimension - nodes, dimension of problem - dim and accuracy eps
should be set before call this function
one can set this values using solver namespace in main code*/

#endif
