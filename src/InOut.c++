/*
 * InOut.c++
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include "Utils.h"

#endif

#ifndef STRUCTS_H
#define STRUCTS_H

#include "Structures.h"

#endif

#include "InOut.h"



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         GETNUMSTATIONS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int getNumStations ()
{
    const std::string dir = "./data/";
    int tempi, ipos, count = 0, nstations = 0;
    std::string fname;
    std::ifstream in_file;
    std::string linetxt;

    fname = dir + "station_latlons.txt";
    in_file.open (fname.c_str (), std::ifstream::in);
    if (in_file.fail ()) {
        std::cout << "***ERROR: Failed to read `station_latlons.txt'***" << 
            std::endl;
    } else {
        getline (in_file, linetxt, '\n');
        count = 0;
        while (getline (in_file, linetxt, '\n')) { count++;	}
        in_file.clear ();
        in_file.seekg (0); // Both lines needed to rewind file.
        for (int i=0; i<=count; i++) {
            getline (in_file, linetxt,'\n');
            ipos = linetxt.find(',',0);
            tempi = atoi (linetxt.substr (0, ipos).c_str());
            if (tempi > nstations) nstations = tempi;
        }
        in_file.close();
    }

    return nstations;
} // end function read_latlons



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           READLATLONS                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void readLatLons (dvec* lons, dvec* lats)
{
    // Presumes file existence confirmed with getNumStations
    const std::string dir = "./data/";
    int count, ipos, tempi, nstations = (*lons).size ();
    std::string fname;
    std::ifstream in_file;
    std::string linetxt;

    for (int i=0; i<nstations; i++) {
        (*lons) (i) = NAN;
        (*lats) (i) = NAN;	}

    fname = dir + "station_latlons.txt";
    in_file.open (fname.c_str (), std::ifstream::in);
    getline (in_file, linetxt, '\n');
    count = 0;
    while (getline (in_file, linetxt, '\n')) { count++;	}
    in_file.clear ();
    in_file.seekg (0); // Both lines needed to rewind file.
    for (int i=0; i<=count; i++) {
        getline (in_file, linetxt,'\n');
        ipos = linetxt.find(',',0);
        tempi = atoi (linetxt.substr (0, ipos).c_str());
        linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        ipos = linetxt.find (',', 0);
        (*lons) (tempi - 1) = atof (linetxt.substr (0, ipos).c_str());
        linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        (*lats) (tempi - 1) = atof (linetxt.c_str());
    }
    in_file.close();
} // end function read_latlons


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          CHECKDATAFILE                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

bool checkDataFile ()
{
    bool check;
    std::ifstream in_file;
    std::string dir = "./data/";

    std::string fname = dir + "raw_data.txt";
    in_file.open (fname.c_str (), std::ifstream::in);
    if (in_file.fail ()) {
        check = false;
        std::cout << "ERROR: File ./data/raw_data.txt does not exist!" << std::endl;
    } else {
        check = true;
    }
    return check;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         CHECKRESULTSFILE                           **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

bool checkResultsFile ()
{
    // Returns false unless file exists and should not be overwritten.
    bool check;
    char ch;
    std::ifstream in_file;
    std::string dir = "./results/";

    std::string fname = dir + "aaaresults_neutral.txt";
    in_file.open (fname.c_str (), std::ifstream::in);
    if (in_file.fail ()) {
        check = false;
    } else {
        check = true;
        std::cout << "WARNING: " << fname << " exists. Overwrite (y/n)? ";
        std::cin >> ch;
        if (tolower (ch) == 'y') check = false;
    }
    return check;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             READDATA                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int readData (dmat* ntrips, bvec* has_data, dvec* lons)
{
    /*
     * Reads the tab-delimited .txt format data provided by the London hire bike
     * scheme. These have 11 columns of:
     * (Journey_Id   Bike_Id Start_Date  Start_Time  End_Date    End_Time    
     * Start_Station  Start_Station_Id    End_Station End_Station_Id
     * Duration_in_Seconds).
     *
     * This routine presumes checkDataFile (), so there is no fail check for .open.
     */
    std::string dir = "./data/";
    int count, ipos, tempi [2], nstations_full = (*lons).size (), nstations;
    std::string fname;
    std::ifstream in_file;
    std::string linetxt;

    fname = dir + "raw_data.txt";
    in_file.open (fname.c_str (), std::ifstream::in);

    getline (in_file, linetxt, '\n');
    count = 0;
    while (getline (in_file, linetxt, '\n')) { count++;	}
    in_file.clear ();
    in_file.seekg (0); // Both lines needed to rewind file.

    (*ntrips) = zmat_d (nstations_full, nstations_full);
    for (int k=0; k<=count; k++) {
        getline (in_file, linetxt,'\n');
        for (int j=0; j<7; j++) {
            ipos = linetxt.find('\t',0);
            linetxt = linetxt.substr (ipos + 1, 
                    linetxt.length () - ipos - 1);
        }
        ipos = linetxt.find ('\t', 0);
        tempi [0] = atoi (linetxt.substr (0, ipos).c_str()); // Start Station ID
        for (int j=0; j<2; j++) {
            linetxt = linetxt.substr (ipos + 1, 
                    linetxt.length () - ipos - 1);
            ipos = linetxt.find ('\t', 0);
        }
        tempi [1] = atoi (linetxt.substr (0, ipos).c_str()); // End Station ID
        if (tempi [0] > 0 && tempi [0] <= nstations_full && tempi [1] > 0 && 
                tempi [1] <= nstations_full &&
                !isnan ((*lons) (tempi [0] - 1)) && 
                !isnan ((*lons) (tempi [1] - 1))) {
            (*ntrips) (tempi [0] - 1, tempi [1] - 1) += 1.0;
        } // end if 
    } // end for k
    in_file.close();

    /* All instances in which ntrips (i, j) == ntrips (j, i) == 0.0 then have no
     * data, so are set to dnix. If only one of these is 0.0, the other one
     * simply has no recorded trips, so is left as 0.0.  The trips also also
     * normalised to the 97 days of data, so it's in trips per day.  Also make a
     * binary index into stations that have non-zero data, and determine
     * nstations -- the number that have data -- which is used from here on
     * instead of NSTATIONS_FULL itself. All data are accordingly resized. */
    for (int i=0; i<nstations_full; i++) {
        (*ntrips) (i, i) = NAN;
    }
    for (int i=0; i<(nstations_full - 1); i++) {
        for (int j=(i + 1); j<nstations_full; j++) {
            if ((*ntrips) (i, j) > 0.0 || (*ntrips) (j, i) > 0.0) {
                (*ntrips) (i, j) = (*ntrips) (i, j) / 97.0;
                (*ntrips) (j, i) = (*ntrips) (j, i) / 97.0;
            }
            else {
                (*ntrips) (i, j) = (*ntrips) (j, i) = NAN;
            }
        }
    }
    nstations = 0;
    for (int i=0; i<nstations_full; i++) {
        count = 0;
        for (int j=0; j<nstations_full; j++) {
            if ((*ntrips) (i, j) >= 0.0) count++;
        }
        if (count > 0) {
            nstations++;
            (*has_data) (i) = true;
        }
        else {
            (*has_data) (i) = false;
        }
    }

    return nstations;
} // end function readData


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          READCLUSTERS                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

ivec readClusters (int num_clusters, int nstations, bool dir_to, int method)
{
    const std::string dir = "../results/";

    int ipos, count;
    ivec cluster_ids;
    std::string fname, linetxt;
    std::ifstream in_file;

    if (method == 0) {
        if (dir_to) { fname = dir + "clust_from_members_ward.txt";  }
        else { fname = dir + "clust_to_members_ward.txt";   }
    } else if (method == 1) {
        if (dir_to) { fname = dir + "clust_from_members_complete.txt";  }
        else { fname = dir + "clust_to_members_complete.txt";   }
    } else if (method == 2) {
        if (dir_to) { fname = dir + "clust_from_members_k-means.txt";  }
        else { fname = dir + "clust_to_members_k-means.txt";   }
    } else if (method == 3) {
        if (dir_to) { fname = dir + "skater_from_members.txt";  }
        else { fname = dir + "skater_to_members.txt";   }
    }
    in_file.open (fname.c_str (), std::ifstream::in);
    if (in_file.fail ()) {
        std::cout << "ERROR: Failed to read " << fname.c_str () << std::endl;
    }
    else {
        cluster_ids.resize (nstations);
        count = 0;
        while (getline (in_file, linetxt, '\n')) count++;
        in_file.clear ();
        in_file.seekg (0);
        count = 0;
        while (getline (in_file, linetxt, '\n')) {
            for (int i=2; i<num_clusters; i++) {
                ipos = linetxt.find (',', 0);
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            }
            ipos = linetxt.find (',', 0);
            cluster_ids (count) = atoi (linetxt.substr (0, ipos).c_str ());
            count++;
        } // end while getline
        in_file.close ();
    } // end else read did not fail.

    return cluster_ids;
} // end function allocateClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            WRITEDATA                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void writeData (stnData *station_data, vectorData *vec_data, dmat *ntrips, 
        dvec *lons, dvec *lats, dmat *st_dists)
{
    // NOTE that this presumes the directory ./results exists, and will crash if
    // not. Directory checks are system dependent, so this way at least remains
    // system independent.
    std::ofstream out_file;

    out_file.open ("./results/ntrips.txt", std::ofstream::out);
    for (int i=0; i<(*station_data).nstations; i++) {
        for (int j=0; j<(*station_data).nstations; j++) {
            out_file << (*ntrips) (i, j);
            if (j < ((*station_data).nstations - 1)) { out_file << ",\t";	}
            else { out_file << std::endl;	}
        }
    }
    out_file.close ();

    out_file.open ("./results/results_list.txt", std::ofstream::out);
    out_file << "i,\tj,\tloni,\tlati,\tlonj,\tlatj,\td,\tr2from,\tr2to" << std::endl;
    for (int i=0; i<(*vec_data).ivec.size (); i++) {
        out_file << (*vec_data).ivec [i] << ",\t" << (*vec_data).jvec [i] << ",\t" << 
            (*vec_data).lonveci [i] << ",\t" << (*vec_data).latveci [i] << ",\t" <<
            (*vec_data).lonvecj [i] << ",\t" << (*vec_data).latvecj [i] << ",\t" <<
            (*vec_data).distvec [i] << ",\t" << (*vec_data).r2from [i] << ",\t" << 
            (*vec_data).r2to [i] << std::endl;
    } // end for i
    out_file.close ();

    out_file.open ("./results/results_st_latlons.txt", std::ofstream::out);
    out_file << "lon,\tlat" << std::endl;
    for (int i=0; i<(*station_data).nstations; i++) {
        out_file << i << ", " << (*lons) (i) << ", " << (*lats) (i) << std::endl;
    }
    out_file.close ();

    out_file.open ("./results/results_r2from.txt", std::ofstream::out);
    for (int i=0; i<(*station_data).nstations; i++) {
        for (int j=0; j<(*station_data).nstations; j++) {
            out_file << (*vec_data).r2from_mat (i, j);
            if (j < ((*station_data).nstations - 1)) { out_file << ",\t";	}
            else { out_file << std::endl;	}
        }
    }
    out_file.close ();

    out_file.open ("./results/results_r2to.txt", std::ofstream::out);
    for (int i=0; i<(*station_data).nstations; i++) {
        for (int j=0; j<(*station_data).nstations; j++) {
            out_file << (*vec_data).r2to_mat (i, j);
            if (j < ((*station_data).nstations - 1)) { out_file << ",\t";	}
            else { out_file << std::endl;	}
        }
    }
    out_file.close ();

    // This is not used, but allows direct comparison of
    // correlations with distances
    out_file.open ("./results/results_dists.txt", std::ofstream::out);
    for (int i=0; i<(*station_data).nstations; i++) {
        for (int j=0; j<(*station_data).nstations; j++) {
            out_file << (*st_dists) (i, j);
            if (j < ((*station_data).nstations - 1)) { out_file << ",\t";	}
            else { out_file << std::endl;	}
        }
    }
    out_file.close ();

} // end writeVectorData
