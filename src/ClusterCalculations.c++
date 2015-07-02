/***************************************************************************
 *  Project:    BikeClusters
 *  File:       ClusterCalculations.c++
 *  Language:   C++
 *
 *  BikeClusters is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  BikeClusters is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 *  more details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  NeutralClusters.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright   Mark Padgham July 2015
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Uses observed statistics from a public bicycle hire system
 *                  to clusters hire stations into a predefined number of random
 *                  clusters, in order to measure the total distance ridden
 *                  between the clusters based on observed rides.  The total
 *                  actual distance ridden can then be used to convert these
 *                  values to a proportion of all trips confined within the
 *                  randomly generated clusters.  The distribution of resultant
 *                  values provides the statistical values with which the actual
 *                  observed values may be compared.
 *
 *  Project Structure:  
 *      Routines are divided between the two main programs:
 *      1. mainNeutral:     Randomly allocates points to clusters to estimate
 *                          `neutral' values of distances ridden within versus
 *                          between clusters 
 *      2. mainActual:      Reads in file of actual cluster memberships
 *                          (generated in R) to determine actual values of
 *                          distances ridden within versus between clusters.
 *      These in turn use the classes ClusterData (to store data as produced from
 *      routines in bike-correlations) and Clusters (to store results of
 *      subsequent calculations). ClusterData includes:
 *      (i)     GetDirName ()
 *      (ii)    GetStations ();
 *      (iii)   MakeStationIndex ();
 *      (iv)    readDMat ();
 *      (v)     readNTrips ();
 *      ClusterCalculations includes:
 *      (i)     allocateClusters ();
 *      (ii)    readClusters ();
 *      (iii)   calcClusterDists ();
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11
 ***************************************************************************/


#include "ClusterCalculations.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       ALLOCATECLUSTERS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


int Clusters::allocateClusters (base_generator_type*  generator)
{
    /*
     * Initially allocates numClusters random points to cluster seeds. Routine
     * then proceeds by (1) selecting a random point from among the
     * non-allocated points; (2) finding closest allocated point to that point;
     * (3) Finding closest unallocated point to this allocated point (which may
     * or may not be the first point); and (4) allocating this point to the
     * cluster of the point identified in step#2. The computational order is
     * thus a few times N, whereas constructing clusters by expanding from
     * random points on the boundaries of allocated points is comparably way
     * more computationally complex, and not worth doing. 
     *
     * NOTE, however, that although this works as outlined, it does produce
     * clusters of very unequal sizes, very often with just one or two clusters
     * having most points, yet others having very, very few points. This arises
     * because as a cluster grows bigger, it becomes proportionately more likely
     * to contain the closest point to any randomly selected point. Growth rates
     * are thus proportional to sizes. The remedy employed here is to grow
     * clusters perfectly evenly by selecting the cluster of the smallest size,
     * and then finding a random point that is closer to that cluster than to
     * any other. Only if this condition fails does the routine revert to simply
     * random point selection regardless of cluster ID.
     */
    int tempi, nearest_in, count [2], this_cluster;
    double dmin, tempd;
    bool check;
    std::vector <int> pt_list;
    clusterIDs.resize (_numStations);

    boost::uniform_real <> uni_dist (0, 1);
    boost::variate_generator <base_generator_type&,
        boost::uniform_real <> > runif ((*generator), uni_dist);
    /* boost::normal_distribution <> norm_dist (0.0, 1.0);
       boost::variate_generator <base_generator_type&,
       boost::normal_distribution <> > rnorm ((*generator), norm_dist); */
    // Burn generator in
    for (int i=0; i<20; i++) { tempd = runif ();	}

    for (int i=0; i<_numStations; i++) clusterIDs (i) = INT_MIN;
    count [0] = _numStations;
    // Set up cluster centres
    for (int i=0; i<numClusters; i++) {
        tempi = floor (runif () * _numStations);
        while (clusterIDs (tempi) > INT_MIN) {
            tempi = floor (runif () * _numStations);
        }
        clusterIDs (tempi) = i;
        count [0]--;
    }
    // The allocate the remaining points
    while (count [0] > 0) {
        tempi = INT_MAX;
        this_cluster = INT_MIN;
        for (int i=0; i<numClusters; i++) {
            count [1] = 0;
            for (int j=0; j<_numStations; j++) {
                if (clusterIDs (j) == i) { count [1]++;    }
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
        for (int i=0; i<_numStations; i++) {
            if (clusterIDs (i) == INT_MIN) {
                tempd = DOUBLE_MAX; 
                tempi = INT_MIN;
                // The search for closest point to i that is in a cluster, and get
                // cluster_num.
                for (int j=0; j<_numStations; j++) {
                    if (clusterIDs (j) > INT_MIN && dists (i, j) < tempd) {
                        tempd = dists (i, j);
                        tempi = j;
                    }
                } // end for j - tempi is the closest point to i that is in a cluster.
                if (clusterIDs (tempi) == this_cluster) {
                    pt_list.push_back (i);
                }
            }
        } // end for i
        if (pt_list.size () == 0) { // Then just select a random point
            tempi = floor (runif () * count [0]);
            count [1] = 0;
            for (int i=0; i<_numStations; i++) {
                if (clusterIDs (i) == INT_MIN) { count [1]++;  }
                if (count [1] == tempi) {
                    tempi = count [1];
                    break;
                }
            }  // end for i
            // tempi is then the random point not in a cluster, and is directly indexed
            // into clusterIDs. The next lines find the nearest cluster.
            dmin = DOUBLE_MAX;
            for (int i=0; i<_numStations; i++) {
                if (clusterIDs (i) > INT_MIN && dists (tempi, i) < dmin) {
                    dmin = dists (tempi, i);
                    nearest_in = i;
                    this_cluster = clusterIDs (i);
                } // end if
            } // end for i
            // Then find point closest to nearest_in that is not in a cluster. This
            // may or may not be the same as tempi above.
            dmin = DOUBLE_MAX;
            for (int i=0; i<_numStations; i++) {
                if (clusterIDs (i) == INT_MIN && dists (nearest_in, i) < dmin) {
                    dmin = dists (nearest_in, i);
                    tempi = i;
                }
            } // end for i
        }
        else { // Pick random point from pt_list 
            tempi = floor (runif () * pt_list.size ());
            tempi = pt_list [tempi];
        } // end else
        // tempi is then simply added to this_cluster
        clusterIDs (tempi) = this_cluster;
        count [0]--;
    } // end while not_in_clust.size () > 0

    pt_list.resize (0);

    return 0;
} // end function allocateClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          READCLUSTERS                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Clusters::readClusters (bool dir_to)
{
    const std::string dir = "./results/";

    int ipos, count;
    ivec cluster_ids;
    std::string fname, linetxt;
    std::ifstream in_file;

    if (dir_to)
        fname = "results/" + _city + "-clust-to-members-ward.txt";
    else
        fname = "results/" + _city + "-clust-from-members-ward.txt";

    in_file.open (fname.c_str (), std::ifstream::in);
    assert (!in_file.fail ());
    
    clusterIDs.resize (_numStations);
    count = 0;
    while (getline (in_file, linetxt, '\n')) count++;
    in_file.clear ();
    in_file.seekg (0);
    count = 0;
    while (getline (in_file, linetxt, '\n')) {
        for (int i=2; i<numClusters; i++) {
            ipos = linetxt.find (',', 0);
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        }
        ipos = linetxt.find (',', 0);
        clusterIDs (count) = atoi (linetxt.substr (0, ipos).c_str ());
        count++;
    } // end while getline
    in_file.close ();

    return 0;
} // end function allocateClusters

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       CALCCLUSTERDISTS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

distStats Clusters::calcClusterDists ()
{
    int count = 0;
    double dist_in, dist_out, tempd;
    distStats clustDists;

    // Note that distance matrices do contain NA values set as DOUBLE_MIN
    clustDists.meanProp = clustDists.sdProp = 0.0;
    clustDists.d_in = clustDists.d_out = 0.0;
    for (int i=0; i<_numStations; i++) {
        for (int j=0; j<_numStations; j++) {
            if (ntrips (i, j) > 0.0 && dists (i, j) >= 0.0) {
                tempd = ntrips (i, j) * dists (i, j);
                if (clusterIDs (i) == clusterIDs (j))
                    clustDists.d_in += tempd;
                else
                    clustDists.d_out += tempd;
            } // end if ntrips > 0.0
        } // end for j
    } // end for i

    for (int i=0; i<numClusters; i++) {
        dist_in = dist_out = 0.0;
        for (int j=0; j<_numStations; j++) {
            for (int k=0; k<_numStations; k++) {
                if (ntrips (j, k) > 0.0 && dists (j, k) >= 0.0) {
                    tempd = ntrips (j, k) * dists (j, k);
                    if (clusterIDs (j) == i && clusterIDs (k) == i) {
                        dist_in += tempd;
                    }
                    else if (clusterIDs (j) == i || clusterIDs (k) == i) {
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
