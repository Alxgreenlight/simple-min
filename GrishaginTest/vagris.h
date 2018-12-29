#ifndef __VAGRIS_H
#define __VAGRIS_H

#ifndef M_PI
#define M_PI 3.141592654
#endif

/* Two-dimensional multiextremal test function of Vladimir A. Grishagin

   See:  Grishagin, V.A., Operating Characteristics of Some Global Search Algorithms.
         in: Problems in Random Search, Riga: Zinatne, 1978, issue 7, pp. 198--206;
   also: Strongin, R.G., and Sergeyev, Ya.D., Global Optimization with Non-Convex
         Constraints: Sequential and Parallel Algorithms. Dordrecht: Kluwer, 2000.
*/

/* First, subroutine set_random() has to be called (only once for a function nf);
   then,  function nf or its derivatives can be evaluated in a point (x,y)
          by using subroutines random_func() or random_pfuncx(), random_pfuncy().
*/



extern double rand_minimums[];            /* global minimum points                  */

void   set_random(int nf);                /* set function number nf, 1<=nf<=100     */
double random_func(double x, double y);   /* evaluate the current function in (x,y) */

double random_pfuncx(double x, double y); /* evaluate the first derivative on x     */
double random_pfuncy(double x, double y); /* evaluate the first derivative on y     */


#endif