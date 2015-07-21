/***************************************************************************
 *  Project:    BikeClusters
 *  File:       mainActual.cc
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
#include "mainActual.h"

int main(int argc, char *argv[]) 
{

    bool dir_to;
    int count, dir_to_i;
    double tempd, sum_mn, sum_sd;
    std::string fname, city = "nyc", method="complete";
    std::ofstream out_file;
    base_generator_type generator(42u);
    time_t seed;
    distStats clustDists;

    // These hold averages for k-means of the 10 repeats
    double meanProp, sdProp, d_in, d_out, d_total;

    std::cout << std::endl << "_____________________________________________" << 
        "____________________________________________" << std::endl;
    std::cout << "|\t\t\t\t\t\t\t\t\t\t\t|" << std::endl;
    std::cout << "|\t./ClustersActual with three parameters:\t\t\t\t\t\t|" << std::endl;
    std::cout << "|\t\t1. <city> =  " <<
        "<london/nyc/boston/chicago/washingtondc>\t\t\t|" << std::endl;
    std::cout << "|\t\t2. <method> = <ward/complete/k-means/skater>\t\t\t\t|" <<
        std::endl;
    std::cout << "|\t\t3. <direction> = <0 for TO, otherwise FROM>\t\t\t\t|" <<
        std::endl;
    std::cout << "|\t\t\t\t\t\t\t\t\t\t\t|" << std::endl;
    std::cout << "_____________________________________________" << 
        "____________________________________________" << std::endl <<
        std::endl;

    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));

    count = 0;
    city = "nyc";

    while (*++argv != NULL)
    {
        if (count < 1)
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
        } else if (count < 2) {
            method = *argv;
            std::transform (method.begin(), method.end(), 
                    method.begin(), ::tolower);
            if (method.substr (0, 2) == "co")
                method = "complete";
            else if (method.substr (0, 2) == "wa")
                method = "ward";
            else if (method.substr (0, 2) == "sk")
                method = "skater";
            else if (method.substr (0, 1) == "k")
                method = "k-means";

        } else {
            dir_to_i = atoi (*argv);
            if (dir_to_i == 0) 
                dir_to = true;
            else 
                dir_to = false;
        }
        count++;
    }
    Clusters clusters (city, method);
    std::cout << "Calculating actual distances ridden for:" << std::endl;
    std::cout << "\tCity = " << city << "; clustering method = " <<
        method << "; direction = ";
    if (dir_to)
        std::cout << "TO";
    else
        std::cout << "FROM";
    std::cout << "; number of stations = " << 
        clusters.returnNumStations () << std::endl;

    fname = city + "-results-actual-";
    if (dir_to)
        fname += "to-" + method + ".txt";
    else
        fname += "from-" + method + ".txt";
    std::cout << "writing to file:" << fname.c_str () << std::endl;
    out_file.open (fname.c_str(), std::ios::out);
    out_file << "nc,\tprop.mn,\tprop.sd,\td.in,\td.out,\tprop.total" << std::endl;
    for (int nc=2; nc <= clusters.returnMaxClustSize (); nc++) {
        clusters.numClusters = nc;
        sum_mn = sum_sd = 0.0;
        if (method != "k-means")
        {
            count = clusters.readClusters (dir_to);
            clustDists = clusters.calcClusterDists ();
            out_file << nc << ",\t" << clustDists.meanProp << ",\t" << 
               clustDists.sdProp << ",\t" << clustDists.d_in << ",\t" <<
               clustDists.d_out << ",\t" << clustDists.d_in / clustDists.d_total <<
               std::endl;
        } else {
            meanProp = sdProp = d_in = d_out = d_total = 0.0;
            for (int i=0; i<10; i++)
            {
                count = clusters.readClusters (dir_to, i);
                clustDists = clusters.calcClusterDists ();
                meanProp += clustDists.meanProp;
                sdProp += clustDists.sdProp;
                d_out += clustDists.d_out;
                d_in += clustDists.d_in;
                d_total += clustDists.d_total;
            }
            out_file << nc << ",\t" << meanProp / 10.0 << ",\t" << 
               sdProp / 10.0 << ",\t" << d_in / 10.0 << ",\t" <<
               d_out / 10.0 << ",\t" << d_in / d_total <<
               std::endl;
        }

        tempd = ((double) nc - 1.0) / 
            ((double) clusters.returnMaxClustSize () - 1.0);
        progLine (tempd, clusters.numClusters);
    }
    std::cout << std::endl << "Total distance travelled = " << 
        clustDists.d_total << " km" << std::endl;
    out_file.close ();
}
