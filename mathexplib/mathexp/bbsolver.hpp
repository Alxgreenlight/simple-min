/* 
 * File:   bbsolver.hpp
 * Author: mikhail
 *
 * Created on February 22, 2019, 9:58 PM
 */

#ifndef BBSOLVER_HPP
#define BBSOLVER_HPP

#include <functional>

/**
 * Generic black box solver interface
 */
template <class T> class BlackBoxSolver {
    public:
        /**
         * The method that searches for a minimum of a function with interval constraints
         * @param n the number of parameters
         * @param x starting point on entry, result on exit (for some methods may be arbitrary on entry)
         * @param a lower (interval) bounds on variables
         * @param b upper (interval) bounds on variables
         * @param f the objective function
         * @return the found value
         */
        virtual T search(int n, T* x, const T* a, const T* b, const std::function<T ( const T* )> &f) = 0;
    
};

#endif /* BBSOLVER_HPP */

