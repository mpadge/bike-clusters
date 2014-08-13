/*
 * DataProcessing.cc
 *
 */

#include "DataProcessing.h"



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          GETALLDISTS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void getAllDists (dvec* lons, dvec* lats, dmat* st_dists, dmat* st_dists_x, 
        dmat* st_dists_y, dmat* st_angles)
{
    int nstations = (*lons).size ();
    DistStruct dists;

    for (int i=0; i<nstations; i++) {
        (*st_dists) (i, i) = 0.0;
        (*st_dists_x) (i, i) = 0.0;
        (*st_dists_y) (i, i) = 0.0;
        (*st_angles) (i, i) = NAN;
    }
    for (int i=0; i<(nstations - 1); i++) {
        for (int j=(i+1); j<nstations; j++) {
            if (!isnan ((*lons) (i)) && (!isnan ((*lons) (j)))) {
                dists = getdists ((*lons) (i), (*lats) (i), 
                        (*lons) (j), (*lats) (j));
                (*st_dists) (i, j) = (*st_dists) (j, i) = dists.d;
                (*st_dists_x) (i, j) = dists.dx;
                (*st_dists_x) (j, i) = -dists.dx;
                (*st_dists_y) (i, j) = dists.dy;
                (*st_dists_y) (j, i) = -dists.dy;
                (*st_angles) (i, j) = calc_angle ((*st_dists_y) (i, j), 
                        (*st_dists_x) (i, j));
                (*st_angles) (j, i) = calc_angle ((*st_dists_y) (j, i), 
                        (*st_dists_x) (j, i));
            } // end if
            else {
                (*st_dists) (i, j) = (*st_dists) (j, i) = NAN;
                (*st_dists_x) (i, j) = (*st_dists_x) (j, i) = NAN;
                (*st_dists_y) (i, j) = (*st_dists_y) (j, i) = NAN;
                (*st_angles) (i, j) = (*st_angles) (j, i) = NAN;
            }
        } // end for j
    } // end for i
} // end function getAllDists


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         GETGRIDDATA                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void getGridData (dvec* lons, dvec* lats, dmat* st_dists, stnData* station_data)
{
    // Measures how far the most isolated station is from its nearest neighbour,
    // and also measures lon & lat ranges of entire system. The nearest neighbour is
    // the minimal inter-station distance, so the desired value is the maximum over
    // all stations of this minimal distance, hence it's called max_minDist.
    int count = 0, nstations = (*lons).size ();
    double dmin;
    (*station_data).max_minDist = -DOUBLE_MAX;
    (*station_data).mean_minDist = 0.0;
    (*station_data).lon_range [0] = (*station_data).lat_range [0] = DOUBLE_MAX;
    (*station_data).lon_range [1] = (*station_data).lat_range [1] = -DOUBLE_MAX;
    for (int i=0; i<nstations; i++) {
        if (!isnan ((*lons) (i)) && (*lons) (i) < (*station_data).lon_range [0]) {
            (*station_data).lon_range [0] = (*lons) (i);
        }
        else if ((*lons) (i) > (*station_data).lon_range [1]) {
            (*station_data).lon_range [1] = (*lons) (i);
        }
        if (!isnan ((*lats) (i)) && (*lats) (i) < (*station_data).lat_range [0]) {
            (*station_data).lat_range [0] = (*lats) (i);
        }
        else if ((*lats) (i) > (*station_data).lat_range [1]) {
            (*station_data).lat_range [1] = (*lats) (i);
        }
        dmin = DOUBLE_MAX;
        for (int j=0; j<nstations; j++) {
            if (j != i && (*st_dists) (i, j) > 0.0 && (*st_dists) (i, j) < dmin) {
                dmin = (*st_dists) (i, j);
            }
        } // end for j
        if (dmin < DOUBLE_MAX) {
            count++;
            (*station_data).mean_minDist += dmin;
            if (dmin > (*station_data).max_minDist) {
                (*station_data).max_minDist = dmin;
            }
        }
    }
    (*station_data).mean_minDist = (*station_data).mean_minDist / (double) count;
    (*station_data).nstations_full = nstations;
} // end function getGridData


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       RESIZELATLONDATA                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void resizeLatLonData (dvec* lons, dvec* lats, bvec* has_data, stnData station_data)
{
    std::vector <double> tempvec1, tempvec2;
    tempvec1.resize (0);
    tempvec2.resize (0);
    for (int i=0; i<station_data.nstations_full; i++) {
        if ((*has_data) (i)) {
            tempvec1.push_back ((*lons) (i));
            tempvec2.push_back ((*lats) (i));
        }
    }
    (*lons).resize (station_data.nstations);
    (*lats).resize (station_data.nstations);
    for (int i=0; i<station_data.nstations; i++) {
        (*lons) (i) = tempvec1 [i];
        (*lats) (i) = tempvec2 [i];
    }
    tempvec1.resize (0);
    tempvec2.resize (0);
} // end function resizeLatLonData



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        RESIZEMATRICES                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void resizeMatrices (dmat* st_dists, dmat* st_dists_x, dmat* st_dists_y, 
        dmat* st_angles, dmat* ntrips, bvec* has_data, stnData station_data)
{
    int tempi [2], nst = station_data.nstations;

    dmat st_dists2 (nst, nst), st_dists_x2 (nst, nst),
         st_dists_y2 (nst, nst), st_angles2 (nst, nst),
         ntrips2 (nst, nst);

    tempi [0] = 0;
    for (int i=0; i<station_data.nstations_full; i++) {
        if ((*has_data) (i)) {
            tempi [1] = 0;
            for (int j=0; j<station_data.nstations_full; j++) {
                if ((*has_data) (j)) {
                    st_dists2 (tempi [0], tempi [1]) = (*st_dists) (i, j);
                    st_dists_x2 (tempi [0], tempi [1]) = (*st_dists_x) (i, j);
                    st_dists_y2 (tempi [0], tempi [1]) = (*st_dists_y) (i, j);
                    st_angles2 (tempi [0], tempi [1]) = (*st_angles) (i, j);
                    ntrips2 (tempi [0], tempi [1]) = (*ntrips) (i, j);
                    tempi [1]++;
                } // end if has_data
            } // end for j
            tempi [0]++;
        } // end if has_data
    } // end for i

    (*st_dists).resize (nst, nst);
    (*st_dists_x).resize (nst, nst);
    (*st_dists_y).resize (nst, nst);
    (*st_angles).resize (nst, nst);
    (*ntrips).resize (nst, nst);
    for (int i=0; i<nst; i++) {
        for (int j=0; j<nst; j++) {
            (*st_dists) (i, j) = st_dists2 (i, j);
            (*st_dists_x) (i, j) = st_dists2 (i, j);
            (*st_dists_y) (i, j) = st_dists_y2 (i, j);
            (*st_angles) (i, j) = st_angles2 (i, j);
            (*ntrips) (i, j) = ntrips2 (i, j);
        }
    }

    st_dists2.resize (0, 0);
    st_dists_x2.resize (0, 0);
    st_dists_y2.resize (0, 0);
    st_angles2.resize (0, 0);
    ntrips2.resize (0, 0);
} // end function resizeMatrices
