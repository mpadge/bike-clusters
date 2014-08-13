/*
 * Calculations.h
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include "Utils.h"

#endif

#ifndef STRUCTS_H
#define STRUCTS_H

#include "Structures.h"

#endif


ivec allocateClusters (int num_clusters, dmat* st_dists, 
    base_generator_type* generator);
distStats calcClusterDists (int nc, ivec cluster_ids, dmat* st_dists, dmat* ntrips);
double vectorAnalyses (vectorData *vec_data, dmat ntrips, stnData station_data, 
        dvec *lons, dvec *lats, dmat *st_dists);
