/*
 * DataProcessing.h
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


void getAllDists (dvec* lons, dvec* lats, dmat* st_dists,
         dmat* st_dists_x, dmat* st_dists_y, dmat* st_angles);
void getGridData (dvec* lons, dvec* lats, dmat* st_dists, stnData* station_data);
void resizeLatLonData (dvec* lons, dvec* lats, bvec* has_data, stnData station_data);
void resizeMatrices (dmat* st_dists, dmat* st_dists_x, dmat* st_dists_y,
        dmat* st_angles, dmat* ntrips, bvec* has_data, stnData station_data);
