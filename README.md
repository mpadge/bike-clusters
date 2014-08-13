# BikeClusters

A collection of programs to analyse clusters within urban areas as produced by
usage patterns of hire bicycle systems. Set up at present to analyse data from
London's system, with sample data provided.

# Building and Usage

1.  Use makefile to produce the three executable files:
  1. `make all`
  2. `make getr2`
  3. `make clusters`
2. Ensure that `./results` directory exists, and that `./data` holds `raw_data.txt` and
   `station_latlons.txt`
3. Run `>./getr2` to compute pairwise inter-station correlations from observed ride
   data. Also makes the following `./results` files:
  1. `ntrips` = Matrix of (double) numbers of trips per day between each pair of stations
  2. `results_list` = Table with 1 row for each station of all summary data
  3. `results_st_latlons` = Just the lat-lon coordinates of each station
  4. `results_r2from/to` = Correlation matrices between all pairs of stations
  5. `results_dists` = Matrix of pairwise distances between all stations
4.  Run R routine `get.clusters` to produce results files of cluster memberships
    using preferred clustering algorithm.
5.  Run R routine `calc.pnc` (can be run independently) to generate
    randomly-distributed probabilities of peak heights in relationship between
    numbers of clusters and their T-values above neutral distributions (in file
    `results_prob_m.txt`).
6.  Run `>./NeutralClusters` to generate random (neutral) values for distances ridden
    between randomly-generated clusters. Produces the file `results_neutral.txt`
7.  Run `>./Clusters` to calculate observed distances ridden between clusters as
    allocated according to specified clustering algorithm (stored in various
    `results_actual_...txt` files).
8.  Run R routine `num.clusts` to generate observed statistics for peak heights to
    be compared with neutral values.
9.  Run R routine `clust.sig` to produce final output of all statistics and
    graphics.

