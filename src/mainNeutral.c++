/***************************************************************************
 *  Project:    BikeClusters
 *  File:       mainNeutral.cc
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

#include "mainNeutral.h"

int main(int argc, char *argv[]) 
{

    int tempi;
    double tempd, sum_mn, sum_sd;
    std::string city = "nyc", cityCaps;
    std::ofstream out_file;
    base_generator_type generator(42u);
    time_t seed;

    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));

    distStats clustDists;

    std::cout << std::endl << "_____________________________________________" << 
        "____________________________________________" << std::endl;
    std::cout << "|\t\t\t\t\t\t\t\t\t\t\t|" << std::endl;
    std::cout << "|\t./ClustersNeutral <city> =  " <<
        "<london/nyc/boston/chicago/washingtondc>\t\t|" << std::endl;
    std::cout << "|\t\t\t\t\t\t\t\t\t\t\t|" << std::endl;
    std::cout << "_____________________________________________" << 
        "____________________________________________" << std::endl <<
        std::endl;

    while (*++argv != NULL)
    {
        city = *argv;
        std::transform (city.begin(), city.end(), city.begin(), ::tolower);
        if (city.substr (0, 2) == "lo")
            city = "london";
        else if (city.substr (0, 2) == "bo")
            city = "boston";
        else if (city.substr (0, 2) == "ch")
            city = "chicago";
        else if (city.substr (0, 2) == "wa" || city.substr (0, 2) == "dc")
            city = "washingtondc";
        else
            city = "nyc";
    }
    cityCaps = city;
    std::transform (cityCaps.begin(), cityCaps.end(), 
            cityCaps.begin(), ::toupper); 
    // method is junk in this call:
    std::string method = "complete";
    Clusters clusters (city, method);
    std::cout << cityCaps << ": Number of stations = " << 
        clusters.returnNumStations () << std::endl;

    std::string fname = city + "-results-neutral.txt";
    out_file.open (fname.c_str(), std::ios::out);
    std::cout << "writing to file:" << fname.c_str () << std::endl;
    std::cout << "Cluster size " << std::endl;
    out_file << "nc,\tdmn,\tdsd" << std::endl;
    for (int nc=2; nc <= clusters.returnMaxClustSize (); nc++) 
    {
        clusters.numClusters = nc;
        sum_mn = sum_sd = 0.0;
        for (int i=0; i<clusters.returnNumRepeats (); i++) 
        {
            tempi = clusters.allocateClusters (&generator);
            clustDists = clusters.calcClusterDists ();
            sum_mn += clustDists.d_in;
            sum_sd += clustDists.d_in * clustDists.d_in;
            tempd = ((double) i + 1.0) / (double) clusters.returnNumRepeats ();
            progLine (tempd, clusters.numClusters);
        }
        sum_mn = sum_mn / (double) clusters.returnNumRepeats ();
        // Note that population variance is calculated, rather than sample
        // variance. The latter requires dividing sum_sd by (n - 1) rather than
        // n, and multiplying _mn * _mn by n / (n - 1).
        sum_sd = sum_sd / (double) clusters.returnNumRepeats () - 
            sum_mn * sum_mn;
        sum_sd = sqrt (sum_sd );
        std::cout << "\r[" << clusters.numClusters <<
            "] | Intra-cluster distance = " << sum_mn << " +/- " <<
            sum_sd << " " << std::endl;
        out_file << nc << ", " << sum_mn << ", " << sum_sd << std::endl;
    } // end for nc 
    out_file.close ();
}
