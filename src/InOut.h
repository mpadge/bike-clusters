/*
 * InOut.h
 *
 */


#include <fstream>
#include <dirent.h>

int getNumStations ();
void readLatLons (dvec* lons, dvec* lats);
bool checkDataFile ();
bool checkResultsFile ();
int readData (dmat* ntrips, bvec* has_data, dvec* lons);
ivec readClusters (int num_clusters, int nstations, bool dir_to, int method);
void writeData (stnData *station_data, vectorData *vec_data, dmat *ntrips, 
        dvec *lons, dvec *lats, dmat *st_dists);
