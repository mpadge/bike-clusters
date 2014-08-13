/*
 *  Project:    NeutralClusters
 *  File:       NeutralClusters.cc
 *  Language:   C++
 *
 *  NeutralClusters is free software: you can redistribute it and/or modify it under
 *  the terms of the GNU General Public License as published by the Free Software
 *  Foundation, either version 3 of the License, or (at your option) any later
 *  version.
 *
 *  NeutralClusters is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  NeutralClusters.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright Mark Padgham 2014
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Uses observed statistics from a public bicycle hire system to
 *                  clusters hire stations into a predefined number of random
 *                  clusters, in order to measure the total distance ridden between
 *                  the clusters based on observed rides.  The total actual distance
 *                  ridden can then be used to convert these values to a proportion
 *                  of all trips confined within the randomly generated clusters.
 *                  The distribution of resultant values provides the statistical
 *                  values with which the actual observed values may be compared.
 *
 *  Project Structure:  Routines are divided between the four programs:
 *      1. InOut:           Routines for input of raw data and output of results.
 *                          Contains:
 *                          (i)     read_latlons
 *                          (ii)    checkDataFile
 *                          (iii)   checkResultsFile
 *                          (iv)    ReadData
 *                          (v)     WriteData
 *      2. DataProcessing:  Routines for processing of input data, primarily storage
 *                          and re-shaping of matrices. Contains:
 *                          (i)     getAllDists
 *                          (ii)    getGridData
 *                          (iii)   resizeLatLonData
 *                          (iv)    resizeMatrices
 *      3. Calculations:    Contains the two primary functions to (i) allocate the
 *                          random clusters, and (ii) calculate distances ridden
 *                          within those clusters. Contains:
 *                          (i)     allocateClusters
 *                          (ii)    calcClusterDists
 *      4. Utils:           Miscellaneous utility functions, including:
 *                          (i)     calc_angle
 *                          (ii)    getdists 
 *                          (iii)   convert_distance 
 *                          (iv)    timeout 
 *
 *  Limitations:
 *
 *  Dependencies:   libboost
 *
 *  Compiler Options:   none
 */

#include "NeutralClusters.h"


int main (int argc, char* argv[])
{
    bool check;
    char ch;
    double tempd, sum_mn, sum_sd;
    distStats clustDists;
    stnData station_data;
    std::ofstream out_file;
    clock_t timer[2];
    time_t seed;
    base_generator_type generator(42u);

    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));
    timer[0] = clock();

    std::cout.setf (std::ios::fixed, std::ios::floatfield);
    std::cout.precision (2);

    // First read the station_latlons and make a lookup table for station numbers.
    station_data.nstations_full = getNumStations ();
    dvec lon (0), lat (0);
    if (station_data.nstations_full > 0)
    {
        lon.resize (station_data.nstations_full);
        lat.resize (station_data.nstations_full);
        for (int i=0; i<station_data.nstations_full; i++) {
            lon(i) = NAN;
            lat(i) = NAN;	}
        readLatLons (&lon, &lat);
        check = true;
    } else check = false;
    if (check) check = checkDataFile (); // DataFile should exist
    if (check) check = !checkResultsFile ();
    // At this point, check=true only if station_latlons exists, raw_data.txt exists,
    // and either aaaresults_neutral.txt doesn't exist, or has been approved to be
    // written over.

    if (check) {
        // Convert all pair-wise lat-lon distances to actual distances in km.
        dmat st_dists (station_data.nstations_full, station_data.nstations_full),
             st_dists_x (station_data.nstations_full, station_data.nstations_full),
             st_dists_y (station_data.nstations_full, station_data.nstations_full),
             st_angles (station_data.nstations_full, station_data.nstations_full);
        getAllDists (&lon, &lat, &st_dists, &st_dists_x, &st_dists_y, &st_angles);
        getGridData (&lon, &lat, &st_dists, &station_data);

        // Then read the actual data and analyse it as it's read.
        std::cout<<"Loading data ... "; std::cout.flush();
        dmat ntrips (station_data.nstations_full, station_data.nstations_full); 
        // (from, to)
        bvec has_data (station_data.nstations_full);
        station_data.nstations = readData (&ntrips, &has_data, &lon);
        // Resize data to only those stations with has_data
        resizeLatLonData (&lon, &lat, &has_data, station_data);
        resizeMatrices (&st_dists, &st_dists_x, &st_dists_y, &st_angles, &ntrips, 
                &has_data, station_data);
        std::cout << " done. nstations = " << station_data.nstations << " / " << 
            station_data.nstations_full << " have data" << std::endl; 

        ivec cluster_ids;

        std::string fname = "./results/results_neutral.txt";
        out_file.open (fname.c_str(), std::ios::out);
        std::cout << "writing to file:" << fname.c_str () << std::endl;
        out_file << "nc,\tdmn,\tdsd" << std::endl;
        for (int nc=2; nc <= MAX_CLUST_SIZE; nc++) {
            sum_mn = sum_sd = 0.0;
            std::cout << "Cluster size = " << nc << " ... " << std::endl;
            for (int i=0; i<NUM_REPEATS; i++) {
                cluster_ids = allocateClusters (nc, &st_dists, &generator);
                clustDists = calcClusterDists (nc, cluster_ids, &st_dists, &ntrips);
                sum_mn += clustDists.d_in;
                sum_sd += clustDists.d_in * clustDists.d_in;
                tempd = ((double) i + 1.0) / (double) NUM_REPEATS;
                progLine (tempd);
            }
            sum_mn = sum_mn / (double) NUM_REPEATS;
            // Note that population variance is calculated, rather than sample
            // variance. The latter requires dividing sum_sd by (n - 1) rather than
            // n, and multiplying _mn * _mn by n / (n - 1).
            sum_sd = sum_sd / (double) NUM_REPEATS - sum_mn * sum_mn ;
            sum_sd = sqrt (sum_sd );
            std::cout << ": Intra-cluster distance = " << sum_mn << " +/- " <<
                sum_sd << std::endl;
            out_file << nc << ", " << sum_mn << ", " << sum_sd << std::endl;
        } // end for nc 
        out_file.close ();
        cluster_ids.resize (0);

        //writeData (&station_data, &ntrips, &lon, &lat, &st_dists);
    } // end if check lat-lon file exists

    timer[1] = clock() - timer[0];
    std::cout<<"Total Calculation Time = ";
    timeout(timer[1] / ((double)CLOCKS_PER_SEC));
    std::cout<<std::endl;

    return 0;
}; // end main
