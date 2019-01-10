#ifndef SOLVER_H
#define SOLVER_H

void checkErrors(int code); //check errors after call min_search and print error info
int min_search(double*a, double *b); 
/*find upper and lower bound for multidimensinal [a,b]
(pointer *compute, amount of nodes per dimension - nodes, dimension of problem - dim and accuracy eps
should be set before call this function
one can set this values using externs of this global values in main code*/

#endif
