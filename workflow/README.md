This directory contain the files needed for the reproduction of all the results given in the paper, including tables and figures. The workflow is splited in two main parts, simulation examples and real data examples.

The following comments are of great importance for the undertanding of the results:

  - Many algorithms depend on random numbers of initializations nstart. This number is set with a low value in some examples, in order to reduce the computational time.
  In other cases, the number of iterations has also been decreased, with the same objective.
  Because of that, the optimal solution presented in the article might not be found, but increasing this value should allow you to achieve the best solution. 
  The examples where this occurs are commented in the code.  

 - Even if we take the same number of random initializations and iterations, techniques used for model validation (both in simulation examples
 and real data examples) are based on random methods, such as Monte Carlo and resampling.
 Results given in that case should slightly vary from the results given in the paper.
