#ifndef _VEC_HPP_
#define _VEC_HPP_
/**
 * Subroutines for frequently used vector operations
 *
 * @author Mikhail Posypkin, ISA RAS, posypkin@isa.ru
 *
 * @warning Using the code below requires an explicit permition from the author
 * @warning Any other use is prohibited
 *
 * @date Copyright by Mikhail Posypkin 2005-2015
 * @file vec.hpp
 */

#include <math.h>
#include <iostream>
#include <sstream>
#include "utilmacro.hpp"

namespace snowgoose {

    class VecUtils {
    public:

        /**
         * Initialize all elements in a vector with scalar values
         * @param x scalar value
         * @param y vector
         */
        template<class T> static void vecSet(int n, T x, T* y) {
            for (int i = 0; i < n; i++)
                y[i] = x;
        }

        /**
         * Adds scalar value to elements of a vector
         * @param x source vector
         * @param y destination vector (can conincide with x)
         * @param a scalar
         */
        template<class T> static void vecAddScalar(int n, T a, T* x, T* y) {
            for (int i = 0; i < n; i++)
                y[i] = x[i] + a;
        }

        /**
         * Compute sum of vector elements
         * @param n dimension
         * @param x vector
         * @return sum of element
         */
        template<class T> static T vecSum(int n, const T* x) {
            T v = 0.;
            for (int i = 0; i < n; i++)
                v += x[i];
            return v;
        }

        /**
         * Compute "first" vector norm 
         * @param n dimension
         * @param x vector
         * @return vector norm
         */
        template<class T> static T vecNormOne(int n, const T* x) {
            T v = 0.;
            for (int i = 0; i < n; i++)
                v = SGMAX(SGABS(x[i]), v);
            return v;
        }

        /**
         * Compute "second" vector norm in square
         * @param n dimension
         * @param x vector
         * @return vector norm
         */
        template<class T> static T vecNormTwoSqr(int n, const T* x) {
            T v = 0.;
            for (int i = 0; i < n; i++)
                v += SGSQR(x[i]);
            return v;
        }

        /**
         * Compute "second" vector norm 
         * @param n dimension
         * @param x vector
         * @return vector norm
         */
        template<class T> static T vecNormTwo(int n, const T* x) {
            return sqrt(vecNormTwoSqr(n, x));
        }

        /**
         * Compute scalar multiple of two vectors
         * @param n dimension
         * @param x first vector                
         * @param y second vector
         * @return scalar multiple
         */
        template<class T> static T vecScalarMult(int n, const T* x, const T* y) {
            T v;
            v = 0;
            for (int i = 0; i < n; i++)
                v += x[i] * y[i];
            return v;
        }

        /**
         * Compute the square distance between two vectors 
         * @param n dimension
         * @param x first vector
         * @param y second vector
         * @return distance
         */
        template<class T> static T vecDist(int n, const T* x, const T* y) {
            T v = 0.;
            for (int i = 0; i < n; i++)
                v += SGSQR(y[i] - x[i]);
            return sqrt(v);
        }

        /**
         * Compute the absolute distance between two vectors
         * @param n dimension
         * @param x first vector
         * @param y second vector
         * @return distance
         */
        template<class T> static T vecDistAbs(int n, const T* x, const T* y) {
            T v = 0.;
            for (int i = 0; i < n; i++) {
                T u = SGABS(y[i] - x[i]);
                v = SGMAX(u, v);
            }
            return v;
        }

        /**
         * Computes the central point between two vectors
         * @param n dimension
         * @param x first vector
         * @param y second vector
         * @param z central point
         */
        template<class T> static T vecCenter(int n, const T* x, const T* y, T* z) {
            for (int i = 0; i < n; i++)
                z[i] = 0.5 * (x[i] + y[i]);
        }

        /**
         * Copies one vector to another
         * @param n dimension
         * @param x source vector
         * @param y destination vector
         */
        template <class T> static void vecCopy(int n, const T * x, T* y) {
            for (int i = 0; i < n; i++)
                y[i] = x[i];
        }

        /**
         * Reverts the vector
         * @param n dimension 
         * @param x vector to revert
         */
        template <class T> static void revert(int n, T * x) {
            for (int i = 0; i < n; i++)
                x[i] *= -1;
        }

        /**
         * Multiplies vector by a scalar y = alpha * x
         * @param n dimension
         * @param x source vector
         * @param alhpa multiple
         * @param y destination vector
         */
        template <class T> static void vecMult(int n, const T * x, T alpha, T* y) {
            for (int i = 0; i < n; i++)
                y[i] = alpha * x[i];
        }

        /**
         * Multiplies vector by a scalar and adds another vector z = x + alpha * y
         * @param n dimension
         * @param x source vector
         * @param y source vector
         * @param alhpa multiple
         * @param z destination vector
         */
        template <class T> static void vecSaxpy(int n, const T * x, const T *y, T alpha, T* z) {
            for (int i = 0; i < n; i++)
                z[i] = x[i] + y[i] * alpha;
        }

       /**
        * Component-wise multiplication of two vectors Multiply 
        * @param n
        * @param x
        * @param y
        * @param z
        */
        template <class T> static void vecMultVect(int n, const T * x, const T *y, T* z) {
            for (int i = 0; i < n; i++)
                z[i] = x[i] * y[i];
        }

        /**
         * Prints the vector to string
         * @param n dimension 
         * @param x vector to print
         * @param prec precision
         * @return resulting string
         */
        template<class T> static std::string vecPrint(int n, const T* x, int prec = 0) {
            std::ostringstream os;
            if (prec) {
                os.precision(prec);
                os.setf(std::ios::fixed, std::ios::floatfield);
            }
            os << "[";
            for (int i = 0; i < n; i++) {
                os << x[i];
                if (i != (n - 1))
                    os << ",";
            }
            os << "]";
            return os.str();
        }

        /**
         * Reads the vector from string
         * @param s string
         * @param n dimension 
         * @param x vector to read
         */
        template<class T> static void vecRead(const std::string s, int n, T* x) {
            std::istringstream is(s);
            for (int i = 0; i < n; i++) {
                is >> x[i];
            }
        }

        /**
         * Computes the maximal value an its position in a vector
         * @param n dimension
         * @param x vector 
         * @param pos address of the maximal element
         * @return the maximal value
         */
        template<class T> static T max(int n, const T *x, int *pos = nullptr) {
            T mv = x[0];
            int p = 0;
            for (int i = 0; i < n; i++) {
                if (x[i] > mv) {
                    p = i;
                    mv = x[i];
                }
            }
            if (pos != nullptr)
                *pos = p;
            return mv;
        }

        /**
         * Computes the mininal element and its position in a vector
         * @param n dimension
         * @param x vector 
         * @param pos address of the minimal element
         * @return the minimal value
         */
        template<class T> static T min(int n, const T *x, int *pos = nullptr) {
            T mv = x[0];
            int p = 0;
            for (int i = 0; i < n; i++) {
                if (x[i] < mv) {
                    p = i;
                    mv = x[i];
                }
            }
            if (pos != nullptr)
                *pos = p;
            return mv;
        }

        /**
         * Computes the element with the maximal absolute value and its position in a vector
         * @param n dimension
         * @param x vector 
         * @param pos address of the element with the maximal absolute value
         * @return the maximal absolute value
         */
        template<class T> static T maxAbs(int n, const T *x, int *pos = nullptr) {
            T mv = SGABS(x[0]);
            int p = 0;
            for (int i = 0; i < n; i++) {
                if (SGABS(x[i]) > mv) {
                    p = i;
                    mv = SGABS(x[i]);
                }
            }
            if (pos != nullptr)
                *pos = p;
            return mv;
        }

        /**
         * Computes the element with the minimal absolute value and its position in a vector
         * @param n dimension
         * @param x vector 
         * @param pos address of the element with the minimal absolute value
         * @return the minimal absolute value
         */
        template<class T> static T minAbs(int n, const T *x, int *pos = nullptr) {
            T mv = SGABS(x[0]);
            int p = 0;
            for (int i = 0; i < n; i++) {
                if (SGABS(x[i]) < mv) {
                    *pos = i;
                    mv = SGABS(x[i]);
                }
            }
            if (pos != nullptr)
                *pos = p;
            return mv;
        }

        /**
         * Projector to the hyperplane
         * @param n dimension
         * @param nplane plane normal vector
         * @param sourcev vector to project
         * @param resultv projected vector     
         */
        template<class T> static void project(int n, const T* nplane, const T* sourcev, T* resultv) {
            T s = vecDotProd(n, nplane, sourcev);
            T csq = vecNormTwoSqr(n, nplane);
            T a = -s / csq;
            vecSaxpy(n, sourcev, nplane, a, resultv);
        }


    };
}
#endif
