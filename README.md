# simple-min

About this suite:
This kit contains a multidimensional multi-criteria global Lipschitz optimization method.

Also in the set there are tests with which you can make sure that the method is operational.
The attached makefiles will help to understand the process of building applications using the developed method.
All test cases contain comments explaining how the solver can be used.

Necessary files:
Solver is located in solver folder. Currently used version is: R_Optim_Pure_Parallel.hpp
Solver files depende on common/bbsolver.hpp because it support interface
from https://github.com/mposypkin/blackbox

Serial solver available at solver/R_Optim.hpp
And parallel one at solver/R_Optim_Pure_Parallel.hpp

Now there are some usable directories:
solver/
	and common/ on which solvers depend
GKLStest/ contains example for running this method for solving GKLS tasks (http://wwwinfo.deis.unical.it/~yaro/GKLS.html)
GKLS_test_tool/ intended to run several experiments with different dimesions using GKLS suite (for researches)
GrishaginTest/ contains example for running this method with a class of two-dimensional test functions developed by V. Grishagin
(http://wwwinfo.deis.unical.it/yaro/Grishagin_web.zip)
mathexplib/ contains example for running this method with well-known real optimization problems, collected in the test suite
(https://www.degruyter.com/view/j/eng.2017.7.issue-1/eng-2017-0050/eng-2017-0050.xml)
P.S. Please take a look to .cpp main files used in Makefiles, some functionality may be turned off
