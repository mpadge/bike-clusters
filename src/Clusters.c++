/*
 * Clusters.cc
 *
 * Reads in a text file with the hierarchical clustering scheme from the actual
 * data, and uses that to calculate the inter-cluster distances ridden as a function
 * of the number of clusters.
 */

#include "Clusters.h"

int main(int argc, char *argv[])
{
    int count, ipos, tempi [2], ttime[2];
    double tempd;
    bool dir_to;
    std::string fname;
    DistStruct dists;
    distStats clustDists;
    RegrResults regr;
    stnData station_data;
    std::ifstream in_file;
    std::ofstream out_file;
    std::string linetxt;
    clock_t timer[2];
    time_t seed;
    base_generator_type generator(42u);

    int dir_to_i = atoi (argv [2]);
    // default of no input sets dir_to = true below, otherwise any non-zero
    // value sets dir_to = false, and analyses movement from
    if (dir_to_i == 0) { dir_to = true; }
    else { dir_to = false;  }
    int clust_method = atoi (argv [3]);
    // same principle: no input sets hclust = true, otherwise false triggers
    // analysis of kmeans clustering results, rather that hierarchical
    // clustering. NOTE: New modified version uses only clust_method thus:
    // clust_method = 0 -> hierarchical clustering with ward distance
    // clust_method = 1 -> hierarchical clustering with complete distance
    // clust_method = 2 -> k-means clustering
    // clust_method = 3 -> skater clustering

    //seed = atoi (argv [3]);
    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));
    timer[0] = clock();

    std::cout.setf (std::ios::fixed, std::ios::floatfield);
    std::cout.precision (2);

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
        station_data.nstations_full << " have data" << std::endl; 

    ivec cluster_ids;
    if (clust_method == 0) {
        if (dir_to) {
            out_file.open ("./results/results_actual_ward_to.txt", std::ios::out);
            std::cout << "Analysing hclust ward for DIR = TO" << std::endl;
        }
        else {
            out_file.open ("./results/results_actual_ward_from.txt", std::ios::out);
            std::cout << "Analysing hclust ward for DIR = FROM" << std::endl;
        }
    } else if (clust_method == 1) {
        if (dir_to) {
            out_file.open ("./results/results_actual_complete_to.txt", std::ios::out);
            std::cout << "Analysing hclust complete for DIR = TO" << std::endl;
        }
        else {
            out_file.open ("./results/results_actual_complete_from.txt", std::ios::out);
            std::cout << "Analysing hclust complete for DIR = FROM" << std::endl;
        }
    } else if (clust_method == 2) { 
        if (dir_to) {
            out_file.open ("./results/results_actual_kmeans_to.txt", std::ios::out);
            std::cout << "Analysing hclust k-means for DIR = TO" << std::endl;
        }
        else {
            out_file.open ("./results/results_actual_kmeans_from.txt", std::ios::out);
            std::cout << "Analysing hclust k-means for DIR = FROM" << std::endl;
        }
    } else { // clust_method == 3
        if (dir_to) {
            out_file.open ("./results/results_actual_skater_to.txt", std::ios::out);
            std::cout << "Analysing skater clusters for DIR = TO" << std::endl;
        }
        else {
            out_file.open ("./results/results_actual_skater_from.txt", std::ios::out);
            std::cout << "Analysing skater clusters for DIR = FROM" << std::endl;
        }
    }

    out_file << "nc,\tprop.mn,\tprop.sd,\td.in,\td.out,\tprop.total" << std::endl;
    cluster_ids.resize (station_data.nstations);
    //std::cout << "Total distance ridden between clusters ... " << std::endl;
    for (int nc=2; nc <= MAX_CLUST_SIZE; nc++) {
        std::cout.flush ();
        cluster_ids = readClusters (nc, station_data.nstations, dir_to, clust_method);
        clustDists = calcClusterDists (nc, cluster_ids, &st_dists, &ntrips);
        out_file << nc << ",\t" << clustDists.meanProp << ",\t" << 
           clustDists.sdProp << ",\t" << clustDists.d_in << ",\t" <<
           clustDists.d_out << ",\t" << clustDists.d_in / clustDists.d_total <<
           std::endl;
        tempd = ((double) nc - 2.0) / ((double) MAX_CLUST_SIZE - 2.0);
        progLine (tempd);
    }

    out_file.close ();
    cluster_ids.resize (0);
    std::cout << std::endl << "Total distance travelled = " << clustDists.d_total << 
        " km" << std::endl;

    //writeData (&station_data, &ntrips, &lon, &lat, &st_dists);
    cluster_ids.resize (0);

    timer[1] = clock() - timer[0];
    std::cout<<"Total Calculation Time = ";
    timeout(timer[1] / ((double)CLOCKS_PER_SEC));
    std::cout<<std::endl;

    return 0;
}; // end main
