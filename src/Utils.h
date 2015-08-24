/***************************************************************************
 *  Project:    BikeClusters
 *  File:       Utils.h
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

#include <stdlib.h> // has abs function
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include <string>
#include <iomanip> // for setfill
#include <sys/ioctl.h> // for console width: Linux only!
#include <ctype.h>
#include <fstream>
#include <assert.h>

#include <boost/config.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/program_options.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#ifndef UTILS_H
#define UTILS_H

#define PI 3.1415926535897932384626433832795

typedef boost::numeric::ublas::vector <int> ivec;
typedef boost::numeric::ublas::matrix <int> imat;
typedef boost::numeric::ublas::vector <double> dvec;
typedef boost::numeric::ublas::matrix <double> dmat;
typedef boost::numeric::ublas::vector <bool> bvec;
typedef boost::numeric::ublas::matrix <bool> bmat;
typedef boost::numeric::ublas::zero_matrix <double> zmat_d;
typedef boost::numeric::ublas::zero_matrix <int> zmat_i;

const double DOUBLE_MAX = std::numeric_limits<double>::max (),
    DOUBLE_MIN = -DOUBLE_MAX,
    FLOAT_MAX = std::numeric_limits <float>::max ();


// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;

struct distStats {
        double meanProp, sdProp, d_in, d_out, d_total;
};

void progLine (double progress, int nc);

#endif
