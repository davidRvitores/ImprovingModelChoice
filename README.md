## Improving Model Choice in Classification: An Approach Based on Clustering of Covariance Matrices.

### Reproducibility

This repository includes the needed material to reproduce the results in the paper "Improving Model Choice in Classification: An Approach Based on Clustering of Covariance Matrices". 

The code in this repository is a first version of the algorithms, created in the software R with the main objective of illustrating the potential of the methodology proposed in simulation and real data examples.  Not all technical and computational details have been covered in this version, and as a consequence a lot of computational time would be required to reproduce all the results of the article. In this repository, we include the code with the algorithms and the workflow needed to reproduce all the results presented in the paper. However, in order to enhance the understanding of the methods, in some examples in the workflow, the number of random starts or iterations has been significantly reduced, to avoid expensive computational times. Appropiate comments are included in these situations in the corresponding code.

### Instructions for use

To reproduce the workflow, the following steps must be taken:

- Run algorithms.R, in the /code directory, with all the libraries in the header. 

For the simulation examples:

-  Load the file simulationData.RData, in the /data directory.

- Run the archive: simulationWorkflow.R

For the real data examples:

-  Load the libraries MASS and pdfCluster, and also the Cancer RNA-Seq Data Set, which can be downloaded from
https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq#

- Run the archive: realDataWorkflow.R







