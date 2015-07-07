# BikeClusters

A collection of programs to analyse clusters within urban areas as produced by
usage patterns of hire bicycle systems. Requires the output of
`bike-correlations', in particular the correlation (R2) and distance matrices.
Currently able to analyse data from London, NYC, Boston, Chicago, and Washington
DC, using cluster methods of `ward', `complete', `k-means', and `skater'.

Analyses are ultimately based on comparisons of the total distances ridden
within clusters to equivalent distances ridden between them.

Use makefile to build. Note that the C++ routines can take a long time to
execute, and so are built as stand-alones, rather than being integrated into R.
The two routines are `ClustersNeutral', which generates neutrally expected
values of inter- and intra-cluster distance, and `Clusters Actual', which does
the corresponding calculations for the observed rides.

All calculations are repeated for a range of numbers of clusters (up to 100).

# Usage for a given city and clustering method:

1. `>./ClustersNeutral city'
2. `R> get.clusters (city method)'
3. `R> get.skater.groups (city, method)'
4. `>./ClustersActual city'
5. 'R> calc.pnc (city, method'
6. 'R> clust.sig (city, method)'

See `aaaread-this' for further details.
