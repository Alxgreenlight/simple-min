# simple-min

About this suite:
This kit contains a multidimensional multi-criteria global Lipschitz optimization method.
The method is presented in three versions: 
sequential algorithm and two different variants of the parallel algorithm.
Also in the set there are tests with which you can make sure that the method is operational.
The attached makefiles will help to understand the process of building applications using the developed method.
All test cases contain comments explaining how the solver can be used.

Necessary files:
Solver is lovated in solver/gridsolver.hpp or solver/gridsolver_hlp.hpp
First one contains serial and one of the parallel algorithms.
Second contains serial and both of parallel algorithms. Versions of
duplicate algorithms in both files are the same.
Solver files depende on common/bbsolver.hpp because it support interface
from https://github.com/mposypkin/blackbox

Now there are some usable directories:
solver/
	and common/ on which solvers depend
GKLStest/ contains example for running this method for solving GKLS tasks (http://wwwinfo.deis.unical.it/~yaro/GKLS.html)
GKLS_test_tool/ intended to run several experiments with different dimesions using GKLS suite (for researches)
GrishaginTest/ contains example for running this method with a class of two-dimensional test functions developed by V. Grishagin
(http://wwwinfo.deis.unical.it/yaro/Grishagin_web.zip)
mathexplib/ contains example for runnig this method with well-known real optimization problems, collected in the test suite
(https://www.degruyter.com/view/j/eng.2017.7.issue-1/eng-2017-0050/eng-2017-0050.xml)
