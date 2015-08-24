/***************************************************************************
 *  Project:    BikeClusters
 *  File:       ClusterData.c++
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


#include "ClusterData.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           GETDIRNAME                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

std::string ClusterData::GetDirName ()
{
    std::ifstream in_file;
    std::string dirtxt;
    std::string configfile = "clusters.cfg"; 
    // Contains name of bike-correlations directory (on single line only)
    struct dirent *ent;

    in_file.open (configfile.c_str (), std::ifstream::in);
    assert (!in_file.fail ());

    getline (in_file, dirtxt, '\n');
    //std::transform (dirtxt.begin(), dirtxt.end(), dirtxt.begin(), ::tolower);

    in_file.close ();
    return dirtxt;
} // end ClusterData::GetDirName



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          GETSTATIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int ClusterData::GetStations ()
{
    int ipos, tempi, count;
    bool tube;
    OneStation oneStation;
    std::string fname;
    std::ifstream in_file;
    std::string linetxt;

    StationList.resize (0);
    count = 0;
    oneStation.name = "";

    if (_city == "london" || _city == "nyc" || _city == "washingtondc")
    {
        fname = _dirName + "data/station_latlons_" + _city + ".txt";
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header
        while (getline (in_file, linetxt,'\n'))
        {
            if (linetxt.length () > 1) 
            {
                ipos = linetxt.find(',',0);
                tempi = atoi (linetxt.substr (0, ipos).c_str());
                if (tempi > count) 
                    count = tempi;
                oneStation.ID = tempi;
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                ipos = linetxt.find (',', 0);
                oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                ipos = linetxt.find (',', 0);
                oneStation.lon = atof (linetxt.substr (0, ipos).c_str());
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                oneStation.name = linetxt;
                StationList.push_back (oneStation);
            }
        }
        in_file.close();
    } 
    else if (_city == "boston")
    {
        // TODO: Store names
        fname = _dirName + "data/hubway_stations.csv";
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header
        while (getline (in_file, linetxt,'\n'))
        {
            ipos = linetxt.find(',',0);
            tempi = atoi (linetxt.substr (0, ipos).c_str());
            if (tempi > count) 
                count = tempi;
            oneStation.ID = tempi;
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            for (int i=0; i<3; i++)
            {
                ipos = linetxt.find (',', 0);
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            }
            ipos = linetxt.find (',', 0);
            oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.lon = atof (linetxt.substr (0, ipos).c_str());
            StationList.push_back (oneStation);
        }
        in_file.close();
    }
    else if (_city == "chicago")
    {
        fname = _dirName + "data/Divvy_Stations_2014-Q3Q4.csv";
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header

        while (getline (in_file, linetxt,'\n'))
        {
            ipos = linetxt.find(',',0);
            tempi = atoi (linetxt.substr (0, ipos).c_str());
            if (tempi > count) 
                count = tempi;
            oneStation.ID = tempi;
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.name = linetxt.substr (0, ipos);
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.lon = atof (linetxt.substr (0, ipos).c_str());
            StationList.push_back (oneStation);
        }
        in_file.close();
    }

    return count;
} // end ClusterData::GetStations


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          MAKESTATIONINDEX                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void ClusterData::MakeStationIndex ()
{
    // First station is #1 and last is _maxStation, so _StationIndex has 
    // len (_maxStns + 1), with _StationIndex [sti.ID=1] = 0 and
    // _StationIndex [sti.ID=_maxStation] = _numStations.
    OneStation sti;

    _StationIndex.resize (_maxStation + 1);
    for (std::vector <int>::iterator pos=_StationIndex.begin();
            pos != _StationIndex.end(); pos++)
        *pos = INT_MIN;
    for (int i=0; i<StationList.size (); i++) 
    {
        sti = StationList [i];
        _StationIndex [sti.ID] = i;
    }
} // end ClusterData::MakeStationIndex


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              READDMAT                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int ClusterData::readDMat ()
{
    std::string distFile;
    if (_city == "nyc")
        distFile = _dirName + "data/station_dists_nyc.txt";
    else 
        distFile = _dirName + "results/stationDistsMat-" + _city + ".csv";

    int count, ipos, tempi [2];
    double d;
    std::ifstream in_file;
    std::string linetxt;
    in_file.open (distFile.c_str (), std::ifstream::in);
    assert (!in_file.fail ());
    in_file.clear ();
    in_file.seekg (0); 
    count = 0;

    if (_city == "nyc")
    {
        while (getline (in_file, linetxt, '\n')) {
            ipos = linetxt.find(',',0);
            tempi [0] = atoi (linetxt.substr (0, ipos).c_str()); // Start Station ID
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find(',',0);
            tempi [1] = atoi (linetxt.substr (0, ipos).c_str()); // End Station ID
            tempi [0] = _StationIndex [tempi[0]];
            tempi [1] = _StationIndex [tempi[1]];
            if (tempi[0] < 0 | tempi[0] > returnNumStations ())
                std::cout << "ERROR: stn[0]#" << tempi[0] << std::endl;
            if (tempi[1] < 0 | tempi[1] > returnNumStations ())
                std::cout << "ERROR: stn[1]#" << tempi[1] << std::endl;
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            d = atof (linetxt.c_str());
            dists (tempi [0], tempi [1]) = dists (tempi [1], tempi [0]) = d;
            count++;
        }
    }
    else 
    {
        /*
         * All other distances, including between train stations, are read
         * direct from matrix
         */
        while (getline (in_file, linetxt, '\n'))
            count++;
        if (count != _numStations)
        {
            std::cout << "ERROR: " << distFile << " does not have " <<
                _numStations << " columns and rows!" << std::endl;
            return -1;
        }
        else
        {
            in_file.clear ();
            in_file.seekg (0); 
            for (int i=0; i<_numStations; i++)
            {
                getline (in_file, linetxt, '\n');
                for (int j=0; j<(_numStations - 1); j++)
                {
                    ipos = linetxt.find (',',0);
                    dists (i, j) = dists (j, i) = 
                        atof (linetxt.substr (0, ipos).c_str());
                    linetxt = linetxt.substr (ipos + 1, 
                            linetxt.length () - ipos - 1);
                } // end for j over columns
                dists (i, _numStations - 1) = dists (_numStations - 1, i) =
                    atof (linetxt.c_str ());
            } // end for i over rows
        } // end else count == _numStations
        // london at this point has some rows of dists with no positive entries
    }

    in_file.close ();

    for (int i=0; i<_numStations; i++)
    {
        d = 0.0;
        for (int j=0; j<_numStations; j++)
            if (dists (i, j) > 0.0)
                d += dists (i, j);
        if (d > 0.0)
            has_data (i) = true;
    }

    return 0;
} // end ClusterData::ReadDMat


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            READNTRIPS                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int ClusterData::readNTrips ()
{
    int ipos, ix, jx;
    std::string fname;
    std::ifstream in_file;
    std::string linetxt;

    if (_city == "london")
        fname = _dirName + "results/NumTrips_london.csv";
    else if (_city == "nyc")
        fname = _dirName + "results/NumTrips_nyc_00.csv";
    else if (_city == "boston")
        fname = _dirName + "results/NumTrips_boston_00.csv";
    else if (_city == "chicago")
        fname = _dirName + "results/NumTrips_chicago_00.csv";
    else if (_city == "washingtondc")
        fname = _dirName + "results/NumTrips_washingtondc.csv";

    in_file.open (fname.c_str (), std::ifstream::in);
    assert (!in_file.fail ());
    in_file.clear ();
    in_file.seekg (0); 

    ix = 0;
    while (getline (in_file, linetxt, '\n')) { 
        jx = 0;
        while ((ipos = linetxt.find (",")) > 0)
        {
            ntrips (ix, jx++) = atof (linetxt.substr (0, ipos).c_str ());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        }
        assert (jx == (_numStations - 1));
        ntrips (ix++, jx) = atof (linetxt.c_str ());
    } 
    in_file.close();
    // just like dists, london at this point has some rows of ntrips with no
    // positive entries (although only 2, while there are 6 with no dists).

    double n;
    for (int i=0; i<_numStations; i++)
    {
        n = 0.0;
        for (int j=0; j<_numStations; j++)
            if (ntrips (i, j) > 0.0)
                n += dists (i, j);
        // This time only reset has_data for stations with dist data, but no
        // ntrips (I don't think there are any of these).
        if (n == 0.0 && has_data (i) == true)
            has_data (i) = false;
    }

    return 0;
} // end function ClusterData::readNTrips

