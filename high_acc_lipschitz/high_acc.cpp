#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <exception>
#include "solver/UltraEstim.hpp"
#include "util/helper.hpp"
extern "C"
{
#include "GKLStest/gkls.h"
#include "GKLStest/rnd_gen.h"
}

#define NUMNOD 1000000

void doGKLS(int dim, int nod);
void doMatx(int nod);

double *a = nullptr, *b = nullptr, *x = nullptr, L;
unsigned long long int fevals;
unsigned long int iters;
std::chrono::steady_clock sc;
std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();
std::ofstream fp;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        //helper::help(ltool);
        return -1;
    }
    std::string param(argv[1]);
    if (param == "g2")
    {
        try
        {
            doGKLS(2, NUMNOD);
        }
        catch (std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
    else if (param == "g3")
    {
        try
        {
            doGKLS(3, NUMNOD);
        }
        catch (std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
    else if (param == "ml")
    {
        doMatx(NUMNOD);
    }
    else
    {
        std::cout << "Parameter is incorrect" << std::endl;
        //helper::help(ltool);
    }

    return 0;
}

std::string print_error_msg(int);

double gFunc(const double *x)
{
    double *xc;
    xc = new double[GKLS_dim];
    for (unsigned int i = 0; i < GKLS_dim; i++)
    {
        xc[i] = x[i];
    }
    double r = GKLS_D_func(xc);
    delete[] xc;
    return r;
}

void doGKLS(int dim, int nod)
{
    std::cout << "GKLS functions, dimension = " << dim << ", with nodes per dimension: " << nod << std::endl;
    int error_code;
    int func_num;
    double UPB;

    GKLS_dim = dim;
    GKLS_num_minima = 10;
    if ((error_code = GKLS_domain_alloc()) != GKLS_OK)
    {
        throw(std::runtime_error(print_error_msg(error_code)));
    }
    GKLS_global_dist = 2.0 / 3.0;
    GKLS_global_radius = 0.5 * GKLS_global_dist;
    GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
    if ((error_code = GKLS_parameters_check()) != GKLS_OK)
    {
        throw(std::runtime_error(print_error_msg(error_code)));
    }

    fp.open("results_g" + std::to_string(dim) + ".txt", std::ios::out);
    if (!fp.is_open())
    {
        throw(std::runtime_error("Problem when trying to open results file"));
    }

    try
    {
        a = new double[GKLS_dim];
        b = new double[GKLS_dim];
        x = new double[GKLS_dim];

        L_accurate<double> La;
        La.stay_fixed(GKLS_dim, nod);

        helper::progress_bar(0, 100);
        for (func_num = 1; func_num <= 100; func_num++)
        {

            if ((error_code = GKLS_arg_generate(func_num)) != GKLS_OK)
            {
                throw(std::runtime_error(print_error_msg(error_code)));
            }

            for (unsigned int i = 0; i < GKLS_dim; i++)
            {
                a[i] = GKLS_domain_left[i];
                b[i] = GKLS_domain_right[i];
            }

            fp << "D-type function number " << func_num << std::endl;
            auto start = sc.now();
            fp << "Real accurate estimated L: " << La.ultraoptimizer_f(a, b, x, UPB, gFunc) << std::endl;
            auto end = sc.now();
            helper::progress_bar(func_num, 100);
            auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            fp << "Evaluation time: " << time_span.count() << " ms" << std::endl;
            fp << std::endl
               << std::endl;
            GKLS_free();
        }
        La.unfix();
    }
    catch (std::exception &e)
    {
        GKLS_free();
        GKLS_domain_free();
        delete[] a;
        delete[] b;
        delete[] x;
        fp.close();
        throw e;
    }
    fp.close();
    delete[] a;
    delete[] b;
    delete[] x;
    GKLS_domain_free();
    return;
}

void doMatx(int nod)
{
    std::cout << "Matexplib functions with nodes per dimension: " << nod << std::endl;
    return;
}

std::string print_error_msg(int error_code)
{
    std::string output;
    switch (error_code)
    {
    case GKLS_OK:
        output = "GKLS_OK: There is no error.";
        break;
    case GKLS_DIM_ERROR:
        output = "GKLS_DIM_ERROR: The problem dimension is out of the valid range [1," + std::to_string((unsigned int)NUM_RND) + "].";
        break;
    case GKLS_NUM_MINIMA_ERROR:
        output = "GKLS_NUM_MINIMA_ERROR: The number of local minima must be greater than 1.";
        break;
    case GKLS_FUNC_NUMBER_ERROR:
        output = "GKLS_FUNC_NUMBER_ERROR: The number of the test function to be generated is out of the range [1,100].";
        break;
    case GKLS_BOUNDARY_ERROR:
        output = "GKLS_BOUNDARY_ERROR: The admissible region boundary vectors are not defined or ill-defined.";
        break;
    case GKLS_GLOBAL_MIN_VALUE_ERROR:
        output = "GKLS_GLOBAL_MIN_VALUE_ERROR: The global minimum value must be greater than " + std::to_string((double)GKLS_PARABOLOID_MIN);
        break;
    case GKLS_GLOBAL_DIST_ERROR:
        output = "GKLS_GLOBAL_DIST_ERROR: The distance from the paraboloid vertex to the global minimizer is too great.";
        break;
    case GKLS_GLOBAL_RADIUS_ERROR:
        output = "GKLS_GLOBAL_RADIUS_ERROR: The radius of the attraction region of the global minimizer is too high.";
        break;
    case GKLS_MEMORY_ERROR:
        output = "GKLS_MEMORY_ERROR: There is not enough memory to allocate.";
        break;
    case GKLS_DERIV_EVAL_ERROR:
        output = "GKLS_DERIV_EVAL_ERROR: An error occurs during derivative evaluation.";
        break;
    case GKLS_FLOATING_POINT_ERROR:
        output = "GKLS_FLOATING_POINT_ERROR: An error occurs in floating point operation";
        break;
    default:
        output = "Unknown error.";
    }
    return output;
}