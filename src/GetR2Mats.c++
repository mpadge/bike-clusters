/*
 * GetR2Mats.c++
 *
 * Correlates numbers of trips between each pair of stations and all other stations.
 */

#include "GetR2Mats.h"


int main(int argc, char *argv[])
{
    const std::string dir = "/data/analyses/bikes_and_cities/londonbikes/data/";
    int count, ipos, tempi [2], ttime[2];
    double tempd;
    std::string fname;
    DistStruct dists;
    RegrResults regr;
    stnData station_data;
    vectorData vec_data;
    std::ifstream in_file;
    std::ofstream out_file;
    std::string linetxt;
    clock_t timer[2];

    //int whichfile = atoi(argv[2]);

    timer[0] = clock();
    //cout.setf(ios::fixed,ios::floatfield);   // floatfield set to fixed
    //cout.precision(2);

    // First read the station_latlons and make a lookup table for station numbers.
    station_data.nstations_full = getNumStations ();
    boost::numeric::ublas::vector<double> lon (station_data.nstations_full);
    boost::numeric::ublas::vector<double> lat (station_data.nstations_full);
    for (int i=0; i<station_data.nstations_full; i++) {
        lon(i) = NAN;
        lat(i) = NAN;	}
    readLatLons (&lon, &lat);

    // Then convert all pair-wise lat-lon distances to actual distances in km.
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
        station_data.nstations_full << " have data... Performing Analyses:" << 
        std::endl;
    std::cout.flush ();

    // Vector analyses
    vec_data.r2from_mat.resize (station_data.nstations, station_data.nstations);
    vec_data.r2to_mat.resize (station_data.nstations, station_data.nstations);
    tempd = vectorAnalyses (&vec_data, ntrips, station_data, &lon, &lat, &st_dists);

    writeData (&station_data, &vec_data, &ntrips, &lon, &lat, &st_dists);

    timer[1] = clock() - timer[0];
    std::cout<<"Total Calculation Time = ";
    timeout(timer[1] / ((double)CLOCKS_PER_SEC));
    std::cout<<std::endl;

    return 0;
}; // end main
