This directory contain the files needed for the reproduction of the results given in the paper, including tables and figures. The workflow is splitted into two main parts, simulation examples and real data examples.

- For the simulation examples, the file simulationData.RData must be loaded.

- For the real data examples, the libraries MASS and pdfCluster must be loaded, and also the Cancer RNA-Seq Data Set, which can be downloaded from
https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq#

The following comments are of great importance for the undertanding of the results:

  - Many algorithms depend on a random numbers of initializations nstart. This number is set to a low value in some examples, in order to decrease the computational time.  In other cases, the number of iterations of various procedures has also been reduced, with the same objective.
  Because of that, the optimal solution presented in the article might not be found, but increasing this value should allow you to achieve the best solution. 
  The examples where this occurs are commented in the code.  
  
  
 - Even if we take the same number of random initializations and iterations, techniques used for model validation (both in simulation examples
 and real data examples) are based on random methods, such as Monte Carlo and resampling.
 Results given in that case would slightly vary from the results given in the paper.
