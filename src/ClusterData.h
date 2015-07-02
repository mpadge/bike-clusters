/***************************************************************************
 *  Project:    BikeClusters
 *  File:       ClusterData.h
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


#ifndef CLUSTERDATA_H
#define CLUSTERDATA_H

#include "Utils.h"

#include <dirent.h>
#include <stdlib.h> // for EXIT_FAILURE
#include <string.h>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <iomanip> // for setfill
#include <sys/ioctl.h> // for console width: Linux only!

class ClusterData
{
    protected:
        std::string _dirName;
        const std::string _city;
        int _numStations, _maxStation;
        std::vector <int> _StationIndex;
    public:
        std::string fileName;
        dmat ntrips, r2, dists;
        ivec clusterIDs;

        struct OneStation 
        {
            std::string name; 
            int ID;
            float lon, lat;
        };
        std::vector <OneStation> StationList;

        ClusterData (std::string str)
            : _city (str)
        {
            _dirName = GetDirName ();
            _maxStation = GetStations ();
            _numStations = StationList.size ();
            InitialiseArrays ();
            if (_city.substr (0, 6) != "oyster")
                MakeStationIndex ();
            readDMat ();
            readNTrips ();
        }
        ~ClusterData ()
        {
            StationList.resize (0);
            ntrips.resize (0, 0);
            r2.resize (0, 0);
            dists.resize (0, 0);
            clusterIDs.resize (0);
        }

        std::string returnDirName () { return _dirName; }
        std::string returnCity () { return _city;   }
        int returnNumStations () { return _numStations; }
        int returnMaxStation () { return _maxStation;   }
        
        std::string GetDirName ();
        int GetStations ();
        void MakeStationIndex ();
        int readDMat ();
        int readNTrips ();

        void InitialiseArrays ()
        {
            ntrips.resize (_numStations, _numStations);
            r2.resize (_numStations, _numStations);
            dists.resize (_numStations, _numStations);
            for (int i=0; i<_numStations; i++)
            {
                for (int j=0; j<_numStations; j++)
                {
                    ntrips (i, j) = DOUBLE_MIN;
                    r2 (i, j) = DOUBLE_MIN;
                    dists (i, j) = DOUBLE_MIN;
                }
                dists (i, i) = 0.0;
            }
        }
}; // end class ClusterData

#endif
