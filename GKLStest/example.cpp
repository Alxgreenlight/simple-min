#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <Windows.h>
#include "solver.h"

extern "C" {
#include "gkls.h"
#include "rnd_gen.h"
}

//externs from solver.cpp
extern int nodes, dim;
extern double eps, UPB, LOB, deltaL, glob;
extern double(*compute)(double* x);

double *a, *b; //sets the search area for task

void print_error_msg  (int);

std::chrono::steady_clock sc; //for runtime measurement


int main()
{
 int error_code;    /* error codes variable */
 int func_num;      /* test function number within a class     */
 FILE *fp;
 char filename[12]; /* name of files */

printf("GKLS-Generator of Classes of ND, D, and D2 Test Functions");
printf("\nfor Global Optimization,");
printf("\n(C) 2002-2003, M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev\n\n");


/* Set the input parameters */
/*if ((error_code=GKLS_set_default()) != GKLS_OK) {
	print_error_msg(error_code);
	return error_code;
}*/
 GKLS_dim = 2;                                             
 GKLS_num_minima = 10;                                     
 if ((error_code = GKLS_domain_alloc()) != GKLS_OK)        
    return error_code;                                     
 GKLS_global_dist = 2.0/3.0;                               
 GKLS_global_radius = 0.5*GKLS_global_dist;                
 GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;                
 if ((error_code = GKLS_parameters_check()) != GKLS_OK)    
    return error_code;                                     

dim = GKLS_dim;

sprintf(filename, "results.txt");

fp = fopen(filename, "wt");
if (fp == NULL) {
	printf("Problem with file...\n");
	system("pause");
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

/* Generate the class of 100 D-type functions */
auto astart = sc.now(); //start time for full set of tests

for (func_num = 1; func_num <= 100; func_num++)
{
  if((error_code=GKLS_arg_generate (func_num)) != GKLS_OK) {
	print_error_msg(error_code);
	return error_code;
  }

  for (int i = 0; i < dim; i++) {
	  a[i] = GKLS_domain_left[i];
	  b[i] = GKLS_domain_right[i];
  } //set search area for test

  compute = GKLS_D_func;
  /* z=GKLS_ND_func(xx); -- for ND-type test function */
  /* z=GKLS_D2_func(xx); -- for D2-type test function */

  printf("\nGenerating the function number %d\n", func_num);

  fprintf(fp,"D-type function number %d", func_num);
  fprintf(fp,"\nof the class with the following parameters:");
  fprintf(fp,"\n    global minimum value   = %f;",GKLS_global_value);

  glob = GKLS_global_value;


  /* Information about global minimizers */
  if (GKLS_glob.gm_index == 0)
	fprintf(fp,"\nAn error during the global minimum searching was occurred!");
  else {
    if (GKLS_glob.num_global_minima == 1) {
	  fprintf(fp,"\n\nThere is one global minimizer.");
  	  fprintf(fp,"\nThe global minimizer coordinates are: ");
	}
	else {
	  fprintf(fp,"\n\nThere are %u global minimizers.",GKLS_glob.num_global_minima);
  	  fprintf(fp,"\nThe global minimizers coordinates are: ");
	}
    for (unsigned int i = 0; i < GKLS_glob.num_global_minima; i++) {
		fprintf(fp, "%d:	", i+1);
	  double *x = GKLS_minima.local_min[GKLS_glob.gm_index[i]];
	  fprintf(fp, "[");
	  for (int k = 0; k < dim; k++) {
		  if (k != dim - 1)
			  fprintf(fp, "%.3lf, ", x[k]);
		  else
			  fprintf(fp, "%.3lf]\n", x[k]);
	  }
	}
  }

  /* Function evaluating */
  auto start = sc.now(); //start time for one function evaluating

  int r = min_search(a, b);

  if (r) {
	  checkErrors(r);
	  fclose(fp);
	  free(a); free(b);
	  system("pause");
	  GKLS_free();
	  GKLS_domain_free();
	  return -1;
  }

  auto end = sc.now();
  auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  fprintf(fp, "\nSearch complete:\nUpper bound: %lf\nReal global minimum: %lf\nLower bound: %lf\nDiff: %lf\n", \
	  UPB, glob, LOB, 1.0*fabs(glob - UPB));
  fprintf(fp, "Evaluation time: %I64u ms\n", time_span.count());
  fprintf(fp, "\n\n");

  
  /* Deallocate memory */
  GKLS_free();
} /* for func_num*/

/* Close files */
fclose(fp);


 /* Deallocate the boundary vectors */
 GKLS_domain_free();
 system("pause");
 return 0;

} /* example.c */


  /* Print an error message */
void print_error_msg(int error_code)
{
	switch (error_code)
	{
	case GKLS_OK:
		printf("\nGKLS_OK: There is no error.");
		break;
	case GKLS_DIM_ERROR:
		printf("\nGKLS_DIM_ERROR: The problem dimension is out of the valid range [1,%u].",
			(unsigned int)NUM_RND);
		break;
	case GKLS_NUM_MINIMA_ERROR:
		printf("\nGKLS_NUM_MINIMA_ERROR: The number of local minima must be greater than 1.");
		break;
	case GKLS_FUNC_NUMBER_ERROR:
		printf("\nGKLS_FUNC_NUMBER_ERROR: The number of the test function to be generated is out of the range [1,100].");
		break;
	case GKLS_BOUNDARY_ERROR:
		printf("\nGKLS_BOUNDARY_ERROR: The admissible region boundary vectors are not defined or ill-defined.");
		break;
	case GKLS_GLOBAL_MIN_VALUE_ERROR:
		printf("\nGKLS_GLOBAL_MIN_VALUE_ERROR: The global minimum value must be greater than %f.",
			(double)GKLS_PARABOLOID_MIN);
		break;
	case GKLS_GLOBAL_DIST_ERROR:
		printf("\nGKLS_GLOBAL_DIST_ERROR: The distance from the paraboloid vertex to the global minimizer is too great.");
		break;
	case GKLS_GLOBAL_RADIUS_ERROR:
		printf("\nGKLS_GLOBAL_RADIUS_ERROR: The radius of the attraction region of the global minimizer is too high.");
		break;
	case GKLS_MEMORY_ERROR:
		printf("\nGKLS_MEMORY_ERROR: There is not enough memory to allocate.");
		break;
	case GKLS_DERIV_EVAL_ERROR:
		printf("\nGKLS_DERIV_EVAL_ERROR: An error occurs during derivative evaluation.");
		break;
	case GKLS_FLOATING_POINT_ERROR:
	default:
		printf("\nUnknown error.");
	}
	printf("\n");
} /* print_error_msg() */