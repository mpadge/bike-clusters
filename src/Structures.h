/*
 * Structures.h
 */

const int MAX_CLUST_SIZE = 100, NUM_REPEATS = 100;
// MAX_CLUST_SIZE must be the same as the value used in the R routine "get.clusters"
// that is used to generate initial cluster memberships.
// NUM_REPEATS is used for NeutralClusters, to determine how many sets of cluster
// memberships are used to generate statistics.

struct stnData {
	int nstations_full, nstations;
	double lon_span, lat_span, lat_range [2], lon_range [2], mean_minDist, 
               max_minDist;
};

struct distStats {
    double meanProp, sdProp, d_in, d_out, d_total;
};

struct vectorData {
    std::vector <int> ivec, jvec;
    std::vector <double> lonveci, lonvecj, latveci, latvecj, distvec, r2from, r2to;
    dmat r2from_mat, r2to_mat;
};
