/*
 * Calculations.cc
 */

#include "Calculations.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       ALLOCATECLUSTERS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


ivec allocateClusters (int num_clusters, dmat* st_dists, 
        base_generator_type*  generator)
{
    /*
     * Initially allocates num_clusters random points to cluster seeds. Routine then
     * proceeds by (1) selecting a random point from among the non-allocated points;
     * (2) finding closest allocated point to that point; (3) Finding closest
     * unallocated point to this allocated point (which may or may not be the first
     * point); and (4) allocating this point to the cluster of the point identified
     * in step#2. The computational order is thus a few times N, whereas
     * constructing clusters by expanding from random points on the boundaries of
     * allocated points is comparably way more computationally complex, and not
     * worth doing. 
     *
     * NOTE, however, that although this works as outlined, it does produce clusters
     * of very unequal sizes, very often with just one or two clusters having most
     * points, yet others having very, very few points. This arises because as a
     * cluster grows bigger, it becomes proportionately more likely to contain the
     * closest point to any randomly selected point. Growth rates are thus
     * proportional to sizes. The remedy employed here is to grow clusters perfectly
     * evenly by selecting the cluster of the smallest size, and then finding a
     * random point that is closer to that cluster than to any other. Only if this
     * condition fails does the routine revert to simply random point selection
     * regardless of cluster ID.
     */
    int nstations = (*st_dists).size1 (); 
    int tempi, nearest_in, count [2], this_cluster;
    double dmin, tempd;
    bool check;
    std::vector <int> pt_list;
    ivec cluster_ids; // The vector of cluster IDs that is returned at end.
    cluster_ids.resize (nstations);

    boost::uniform_real <> uni_dist (0, 1);
    boost::variate_generator <base_generator_type&,
        boost::uniform_real <> > runif ((*generator), uni_dist);
    /* boost::normal_distribution <> norm_dist (0.0, 1.0);
       boost::variate_generator <base_generator_type&,
       boost::normal_distribution <> > rnorm ((*generator), norm_dist); */
    // Burn generator in
    for (int i=0; i<20; i++) { tempd = runif ();	}

    for (int i=0; i<nstations; i++) cluster_ids (i) = INT_MIN;
    count [0] = nstations;
    // Set up cluster centres
    for (int i=0; i<num_clusters; i++) {
        tempi = floor (runif () * nstations);
        while (cluster_ids (tempi) > INT_MIN) {
            tempi = floor (runif () * nstations);
        }
        cluster_ids (tempi) = i;
        count [0]--;
    }
    // The allocate the remaining points
    while (count [0] > 0) {
        tempi = INT_MAX;
        this_cluster = INT_MIN;
        for (int i=0; i<num_clusters; i++) {
            count [1] = 0;
            for (int j=0; j<nstations; j++) {
                if (cluster_ids (j) == i) { count [1]++;    }
            }
            if (count [1] < tempi) {
                tempi = count [1];
                this_cluster = i;
            }
        } // end for i
        // this_cluster then identifies the smallest cluster, with the following
        // lines finding all points that are closer to that cluster than to any
        // other.
        pt_list.resize (0);
        for (int i=0; i<nstations; i++) {
            if (cluster_ids (i) == INT_MIN) {
                tempd = DOUBLE_MAX; 
                tempi = INT_MIN;
                // The search for closest point to i that is in a cluster, and get
                // cluster_num.
                for (int j=0; j<nstations; j++) {
                    if (cluster_ids (j) > INT_MIN && (*st_dists) (i, j) < tempd) {
                        tempd = (*st_dists) (i, j);
                        tempi = j;
                    }
                } // end for j - tempi is the closest point to i that is in a cluster.
                if (cluster_ids (tempi) == this_cluster) {
                    pt_list.push_back (i);
                }
            }
        } // end for i
        if (pt_list.size () == 0) { // Then just select a random point
            tempi = floor (runif () * count [0]);
            count [1] = 0;
            for (int i=0; i<nstations; i++) {
                if (cluster_ids (i) == INT_MIN) { count [1]++;  }
                if (count [1] == tempi) {
                    tempi = count [1];
                    break;
                }
            }  // end for i
            // tempi is then the random point not in a cluster, and is directly indexed
            // into cluster_ids. The next lines find the nearest cluster.
            dmin = DOUBLE_MAX;
            for (int i=0; i<nstations; i++) {
                if (cluster_ids (i) > INT_MIN && (*st_dists) (tempi, i) < dmin) {
                    dmin = (*st_dists) (tempi, i);
                    nearest_in = i;
                    this_cluster = cluster_ids (i);
                } // end if
            } // end for i
            // Then find point closest to nearest_in that is not in a cluster. This
            // may or may not be the same as tempi above.
            dmin = DOUBLE_MAX;
            for (int i=0; i<nstations; i++) {
                if (cluster_ids (i) == INT_MIN && (*st_dists) (nearest_in, i) < dmin) {
                    dmin = (*st_dists) (nearest_in, i);
                    tempi = i;
                }
            } // end for i
        }
        else { // Pick random point from pt_list 
            tempi = floor (runif () * pt_list.size ());
            tempi = pt_list [tempi];
        } // end else
        // tempi is then simply added to this_cluster
        cluster_ids (tempi) = this_cluster;
        count [0]--;
    } // end while not_in_clust.size () > 0

    pt_list.resize (0);

    return cluster_ids;
} // end function allocateClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       CALCCLUSTERDISTS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

distStats calcClusterDists (int nc, ivec cluster_ids, dmat* st_dists, dmat* ntrips)
{
    int count = 0, nstations = (*st_dists).size1 ();
    double dist_in, dist_out, tempd;
    distStats clustDists;

    clustDists.meanProp = clustDists.sdProp = clustDists.d_in = clustDists.d_out = 0.0;
    for (int i=0; i<nstations; i++) {
        for (int j=0; j<nstations; j++) {
            if ((*ntrips) (i, j) > 0.0) {
                    tempd = (*ntrips) (i, j) * (*st_dists) (i, j);
                if (cluster_ids (i) == cluster_ids (j)) {
                    clustDists.d_in += tempd;
                } else {
                    clustDists.d_out += tempd;
                }
            } // end if ntrips > 0.0
        } // end for j
    } // end for i

    for (int i=0; i<nc; i++) {
        dist_in = dist_out = 0.0;
        for (int j=0; j<nstations; j++) {
            for (int k=0; k<nstations; k++) {
                if ((*ntrips) (j, k) > 0.0) {
                    tempd = (*ntrips) (j, k) * (*st_dists) (j, k);
                    if (cluster_ids (j) == i && cluster_ids (k) == i) {
                        dist_in += tempd;
                    }
                    else if (cluster_ids (j) == i || cluster_ids (k) == i) {
                        dist_out += tempd;
                    }
                } // end if ntrips > 0
            } // end for k
        } // end for j
        if (dist_in > 0.0 || dist_out > 0.0) {
            count++;
            tempd = dist_in / (dist_in + dist_out);
            clustDists.meanProp += tempd;
            clustDists.sdProp += tempd * tempd;
        }
    }
    clustDists.meanProp = clustDists.meanProp / (double) count;
    clustDists.sdProp = clustDists.sdProp / (double) count - 
        clustDists.meanProp * clustDists.meanProp;
    clustDists.sdProp = clustDists.sdProp * (double) count / ((double) count - 1.0);
    clustDists.sdProp = sqrt (clustDists.sdProp);
    clustDists.d_total = clustDists.d_in + clustDists.d_out;

    return clustDists;
} // end function calcClusterDists



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        VECTORANALYSES                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

double vectorAnalyses (vectorData *vec_data, dmat ntrips, stnData station_data, 
        dvec *lons, dvec *lats, dmat *st_dists)
{
    double tempd;
    RegrResults regr1, regr2;
    std::vector <double> vec_fromA, vec_toA, vec_fromB, vec_toB;

    (*vec_data).ivec.resize (0);
    (*vec_data).jvec.resize (0);
    (*vec_data).lonveci.resize (0);
    (*vec_data).latveci.resize (0);
    (*vec_data).lonvecj.resize (0);
    (*vec_data).latvecj.resize (0);
    (*vec_data).distvec.resize (0);
    (*vec_data).r2from.resize (0);
    (*vec_data).r2to.resize (0);

    for (int i=0; i<station_data.nstations; i++) {
        (*vec_data).r2from_mat (i, i) = (*vec_data).r2to_mat (i, i) = DOUBLE_MIN;
    }
    for (int i=0; i<(station_data.nstations - 1); i++) {
        for (int j=(i + 1); j<station_data.nstations; j++) {
            (*vec_data).r2from_mat (i, j) = (*vec_data).r2from_mat (j, i) = 
                (*vec_data).r2to_mat (i, j) = (*vec_data).r2to_mat (j, i) = DOUBLE_MIN;
        }
    }

    int count = 0;
    for (int i=0; i<(station_data.nstations - 1); i++) {
        for (int j=(i + 1); j<station_data.nstations; j++) {
            vec_fromA.resize (0);
            vec_fromB.resize (0);
            vec_toA.resize (0);
            vec_toB.resize (0);
            for (int k=0; k<station_data.nstations; k++) {
                if (k != i && k != j) {
                    if (ntrips (i, k) > 0.0 && ntrips (j, k) > 0.0) {
                        vec_fromA.push_back (ntrips (i, k));
                        vec_fromB.push_back (ntrips (j, k));
                    }
                    if (ntrips (k, i) > 0.0 && ntrips (k, j) > 0.0) {
                        vec_toA.push_back (ntrips (k, i));
                        vec_toB.push_back (ntrips (k, j));
                    }
                } // end if k != i
            } // end for k
            if (vec_fromA.size () > 2 && vec_toA.size () > 2) {
                regr1 = regression (vec_fromA, vec_fromB);
                regr2 = regression (vec_toA, vec_toB);
                if (!isnan (regr1.r2) && !isnan (regr2.r2)) {
                    count++;
                    (*vec_data).r2from.push_back (regr1.r2);
                    (*vec_data).r2from_mat (i, j) = 
                        (*vec_data).r2from_mat (j, i) = regr1.r2;
                    (*vec_data).r2to.push_back (regr2.r2);
                    (*vec_data).r2to_mat (i, j) = 
                        (*vec_data).r2to_mat (j, i) = regr2.r2;
                    (*vec_data).ivec.push_back (i);
                    (*vec_data).jvec.push_back (j);
                    (*vec_data).lonveci.push_back ((*lons) (i));
                    (*vec_data).latveci.push_back ((*lats) (i));
                    (*vec_data).lonvecj.push_back ((*lons) (j));
                    (*vec_data).latvecj.push_back ((*lats) (j));
                    (*vec_data).distvec.push_back ((*st_dists) (i, j));
                }
            }
        } // end for j
        tempd = (double) i / ((double) station_data.nstations - 1.0);
        progLine (tempd);
    } // end for i
    progLine (1.0);
    std::cout << std::endl << count;
    count = station_data.nstations * (station_data.nstations - 1) / 2;
    std::cout << " / " << count << " pair-wise comparisons calculated." <<
        std::endl;

    vec_fromA.resize (0);
    vec_toA.resize (0);
    for (int i=0; i<(station_data.nstations - 1); i++) {
        for (int j=(i + 1); j<station_data.nstations; j++) {
            if (ntrips (i, j) > 0.0 && ntrips (j, i) > 0.0) {
                vec_fromA.push_back (ntrips (i, j));
                vec_toA.push_back (ntrips (j, i));
            }
        } // end for j
    } // end for i
    regr1 = regression (vec_fromA, vec_toA);

    vec_fromA.resize (0);
    vec_fromB.resize (0);
    vec_toA.resize (0);
    vec_toB.resize (0);

    return regr1.r2;
} // end function vectorAnalyses

