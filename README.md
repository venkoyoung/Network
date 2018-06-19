# Network scripts

The idea of this folder is contain the script for 3 different steps regarding the network analysis:

1.- Estimation of Rho parameter
Input file: data filtered and prepared. 
  -Not NAS, they were deleted or imputed previously
  -Filtered according variance or something
 Input parameter: 
 Number of random matrixes to be computed for the estimation of FDR
 range and interval of RHO scanning. (Rho goes from 0 to 1). 
 I suggest to start from 0.4 to 0.7, with intervals of 0.5. Then to reduce the range and interval (to 0.05 for example)
 
Given a matrix, this functions are useful for randomize a matrix by rows and compute true and random edges  upon a network reconstruction


Output plots:
With this data we analyze:
-Number of real edges versus random edges
-Number of real positive values versus random positive values
 

2.- Network reconstruction.
Upon definition of RHO parameter for a fdr of 10 or 5 %, build (the correlation) network 
Input file: data filtered and prepared. (Same as 1)
Output files: 
    -Network plot
    -Table with topologies description
    -Edge list


3.- Network comparison
Given 2 networks, extract relevant information about their overlap
