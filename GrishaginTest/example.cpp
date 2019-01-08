#include <cstdio>
#include <chrono>
#include <iostream>
#include <vector>
#include "../solver/solver.h"
extern "C" {
#include "vagris.h"
}

//externs from solver.cpp
extern int nodes, dim;
extern double eps, UPB, LOB, deltaL, glob;
extern double(*compute)(double x1, double x2);

double *a, *b; //sets the search area for task

std::chrono::steady_clock sc; //for runtime measurement

int main()
{
 int nf;
 char filename[12];
 FILE *fp;
 sprintf(filename,"results.txt");

 fp = fopen(filename, "wt");
 if (fp == NULL) {
	 printf("Problem with file...\n");
	 std::cin.get();
	 return 1;
 }
 
 printf("Set number of nodes per dimension\nIt can significantly affect on performance!\n");
 while (scanf("%d", &nodes) != 1) {
	 printf("Please, repeat input\n");
 }

 printf("Set accuracy\nIt can significantly affect on performance!\n");
 while (scanf("%lf", &eps) != 1) {
	 printf("Please, repeat input\n");
 }

 a = (double*)malloc(dim * sizeof(double));
 b = (double*)malloc(dim * sizeof(double));

 for (int i = 0; i < dim; i++) {
	 a[i] = 0.0;
	 b[i] = 1.0;
 } //set search area for all tests to [0..1] [0..1] (like in tutorial for this tests)

 compute = random_func; //random_func returns value of function in x1,x2 (from vagris.h)

 auto astart = sc.now(); //start time for full set of tests
for (nf = 1; nf <= 100; nf++)
	 {
		 set_random(nf);

		 fprintf(fp, "\nFunction number %d", nf);
		 fprintf(fp, "\nIts global minimizer is (%lf, %lf)\n", rand_minimums[(nf - 1) * 2], rand_minimums[(nf - 1) * 2 + 1]);
			 glob = random_func(rand_minimums[(nf - 1) * 2], rand_minimums[(nf - 1) * 2 + 1]);
		 printf("\nComputing function %d...\n", nf);

		 auto start = sc.now(); //start time for one function evaluating

		 int r = min_search(a, b);

		 if (r) {
			 checkErrors(r);
			 fclose(fp);
			 free(a); free(b);
			 std::cin.get();
			 return -1;
		 }

		 auto end = sc.now();
		 auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		 fprintf(fp, "Search complete:\nUpper bound: %lf\nReal global minimum: %lf\nLower bound: %lf\nDiff: %lf\n", \
			 UPB, glob, LOB, 1.0*fabs(glob - UPB));
		 fprintf(fp, "Evaluation time: %I64u ms\n", time_span.count());
		 fprintf(fp, "\n\n");
} /* for nf*/

 printf("\n");

 auto aend = sc.now();
 auto atime_span = std::chrono::duration_cast<std::chrono::seconds>(aend - astart);

 fprintf(fp, "Total evaluation time: %I64u s\n", atime_span.count());
 printf("Total evaluation time: %I64u s\n", atime_span.count());

 fclose(fp);
 free(a); free(b);
 std::cin.get();
 return 0;
}






