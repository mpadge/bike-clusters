/***************************************************************************
 *  Project:    BikeClusters
 *  File:       ClusterCalculations.h
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
 *  BikeClusters.  If not, see <http://www.gnu.org/licenses/>.
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


#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "Utils.h"
#include "ClusterData.h"



class Clusters : public ClusterData
{
    private:
        const int _MAX_CLUST_SIZE = 100, _NUM_REPEATS = 100;
        /*
         * MAX_CLUST_SIZE must be the same as the value used in the R routine
         * "get.clusters" that is used to generate initial cluster memberships.
         * NUM_REPEATS is used for NeutralClusters, to determine how many sets
         * of cluster memberships are used to generate statistics.
         */
    protected:
        const std::string _method;
    public:
        int numClusters;
        std::string clustMethod;

        Clusters (std::string cityStr, std::string methodStr) 
            : _method (methodStr), ClusterData (cityStr)
        {
            clustMethod = "complete";
        }
        ~Clusters ()
        {
            clusterIDs.resize (0);
        }

        int returnMaxClustSize () 
        {
            return _MAX_CLUST_SIZE; 
        }
        int returnNumRepeats () 
        {
            return _NUM_REPEATS;
        }
        std::string returnMethod () 
        {
            return _method; 
        }

        int allocateClusters (base_generator_type * generator);
        int readClusters (bool dir_to, int trialNum = -1);
        distStats calcClusterDists ();
}; // end class Clusters

#endif
