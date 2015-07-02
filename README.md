# BikeClusters

A collection of programs to analyse clusters within urban areas as produced by
usage patterns of hire bicycle systems. Requires the output of
`bike-correlations', in particular the correlation (R2) and distance matrices.
Currently able to analyse data from London, NYC, Boston, Chicago, and Washington
DC.

# Building and Usage

1.  Use makefile to produce the two executable files ClustersNeutral and
    ClustersActual 
2.  Run R routine `get-clusters` to produce results files of cluster memberships
    using preferred clustering algorithm.
3.  Run R routine `calc.pnc` (can be run independently) to generate
    randomly-distributed probabilities of peak heights in relationship between
    numbers of clusters and their T-values above neutral distributions (in file
    `results_prob_m.txt`).
4.  Run `>./ClustersNeutral` to generate random (neutral) values for distances
    ridden between randomly-generated clusters. 
5.  Run `>./ClustersActual` to calculate observed distances ridden between clusters as
    allocated according to specified clustering algorithm 
6.  Run R routine `num.clusts` to generate observed statistics for peak heights to
    be compared with neutral values.
7.  Run R routine `clust.sig` to produce final output of all statistics and
    graphics.

