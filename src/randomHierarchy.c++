/***************************************************************************
 *  Project:    BikeClusters
 *  File:       randomHierarchy.c++
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
 *  randomHierarchy.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright   Mark Padgham July 2015
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    This routine generates random clusters to quantify the
 *                  probabilities of the largest continuous components being of
 *                  equal or larger size than observed clusters. The observed
 *                  values are entered directly into this code. 
 *
 *                  Results in a list of T-statistics for each chosen number of
 *                  total clusters, and contiguous groups from selected numbers
 *                  within those clusters. (That is, for each number of
 *                  clusters, a given number are subsequently partitioned, with
 *                  the observed data revealing that some number less than this
 *                  latter number are contiguous.)
 *
 *  Limitations:
 *
 *  Dependencies:       libboost CGAL gmp
 *
 *  Compiler Options:   -std=c++11 -lCGAL -lgmp
 ***************************************************************************/


#include "randomHierarchy.h"

/**********************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main (int argc, char *argv []) 
{

    bool dir_to;
    int tempi = 0, dir_to_i;
    double tempd, mn, sd, tt;
    std::string fname, city = "nyc", method="complete";
    std::ofstream out_file;
    time_t seed;
    base_generator_type generator(42u);

    while (*++argv != NULL)
    {
        if (tempi < 1)
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
        } else if (tempi < 2) {
            method = *argv;
            std::transform (method.begin(), method.end(), method.begin(), 
                    ::tolower);
            if (method.substr (0, 2) == "co")
                method = "complete";
            else
                method = "k-means";
        } else {
            dir_to_i = atoi (*argv);
            if (dir_to_i == 0) 
                dir_to = true;
            else 
                dir_to = false;
        }
        tempi++;
    }

    //seed = atoi (argv [3]);
    generator.seed (static_cast <unsigned int> (seed));

    ClusterData clusterData (city, dir_to, method);

    // Allocate the random points with CGAL in square with values in +/- first
    // argument, and size passed implicitly as npts.
    CGAL::Random_points_in_square_2<Point> g(1.0); // random points generator
    Points_with_id points;
    Delaunay tr;
    for (int i=0; i<clusterData.npts; i++) 
    {
        points.push_back (std::make_pair (*g, i));
        g++;
    }
    // points can then be accessed by .first, which has [0:1] for (x,y), and
    // .second, which is the ID number
    clusterData.getNeighbours (&points);
    clusterData.getdists (&points);

    ivec counts (clusterData.numClusts);
    for (int i=0; i<clusterData.numClusts; i++) 
    {
        mn = sd = 0.0;
        counts (i) = 0;
        for (int j=0; j<clusterData.num_repeats; j++) 
        {
            int mincount = 0;
            while (mincount < 5) 
            {
                clusterData.allocateClusters (clusterData.clusterNumbers [i], 
                        &generator);
                clusterData.getTable (&clusterData.cluster_ids);
                mincount = 0;
                for (int k=0; k<clusterData.table.size (); k++) 
                    if (clusterData.table (k) > mincount) 
                        mincount = clusterData.table (k);
            } // end while mincout < 5

            tempi = clusterData.getContiguousClusters (&points,
                    clusterData.clusterNumbers (i),
                    clusterData.numTot (i), &generator);
            mn += (double) tempi;
            sd += (double) tempi * (double) tempi;

            if (tempi >= clusterData.numContig (i)) 
                counts (i)++;
        } // end for j
        mn = mn / (double) clusterData.num_repeats;
        sd = sd / ((double) clusterData.num_repeats - 1.0) - mn * mn;
        tt = (mn - (double) clusterData.numContig (i)) * 
            sqrt ((double) clusterData.num_repeats) / sqrt (sd);

        // Calculate cumulative probability:
        tempd = 1.0;
        for (int j=0; j<=i; j++)
            tempd = tempd * (double) counts (j) / 
                (double) clusterData.num_repeats;

        std::cout << "[" << i << ", nc=" << clusterData.clusterNumbers (i) << 
            "]: contig = " << mn << " +/- " << sd << " / observed = " <<
            clusterData.numContig (i) << "; (N, p, p_cum) = (" << counts (i) <<
            ", " << (double) counts (i) / (double) clusterData.num_repeats << 
            ", " << tempd << "); T = " << tt << std::endl;
    } // end for i

    /* 
     * Note that for N=1000 trials, the critical T-value is around 1.65. All
     * values are above this, so the overall p-values are really, really low.
     * The T-values can be converted to probabalities with pt, which gives the
     * area to the left of T. For one-tailed tests, 1-pt will thus give the area
     * to the right.  > 1 - pt (T, 1000)
     */

    return 0;
}



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         SUBFUNCTIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       READNUMCLUSTERS                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void ClusterData::readNumClusters ()
{
    int ipos, coli = 0; 
    std::string linetxt;
    std::ifstream in_file;

    in_file.open (_fname.c_str(), std::ifstream::in);
    assert (!in_file.fail ());
    if (!_dir_to)
        coli = 6;
    if (_method != "complete")
        coli += 3;
    getline (in_file, linetxt, '\n'); // header

    numClusts = 0;
    std::vector <int> clusterNumbersTemp, numContigTemp, numTotTemp;
    while (getline (in_file, linetxt, '\n'))
    {
        for (int i=0; i<coli; i++)
        {
            ipos = linetxt.find (',', 0);
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        }
        ipos = linetxt.find (',', 0);
        clusterNumbersTemp.push_back (atoi (linetxt.substr (0, ipos).c_str()));
        if (clusterNumbersTemp.back () > 0)
            numClusts++;
        linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        ipos = linetxt.find (',', 0);
        numContigTemp.push_back (atoi (linetxt.substr (0, ipos).c_str()));
        linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        if (_dir_to)
        {
            ipos = linetxt.find (',', 0);
            numTotTemp.push_back (atoi (linetxt.substr (0, ipos).c_str()));
        } else {
            numTotTemp.push_back (atoi (linetxt.c_str()));
        }
    }
    in_file.close ();

    clusterNumbers.resize (numClusts);
    numContig.resize (numClusts);
    numTot.resize (numClusts);
    for (int i=0; i<numClusts; i++)
    {
        clusterNumbers (i) = clusterNumbersTemp [i];
        numContig (i) = numContigTemp [i];
        numTot (i) = numTotTemp [i];
    }

    clusterNumbersTemp.resize (0);
    numContigTemp.resize (0);
    numTotTemp.resize (0);
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       ALLOCATECLUSTERS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void ClusterData::allocateClusters (int num_clusters, 
        base_generator_type * generator)
{
    /*
     * Initially allocates num_clusters random points to cluster seeds. Routine
     * then proceeds by (1) selecting a random point from among the
     * non-allocated points; (2) finding closest allocated point to that point;
     * (3) Finding closest unallocated point to this allocated point (which may
     * or may not be the first point); and (4) allocating this point to the
     * cluster of the point identified in step#2. The computational order is
     * thus a few times N, whereas constructing clusters by expanding from
     * random points on the boundaries of allocated points is comparably way
     * more computationally complex, and not worth doing. 
     *
     * NOTE, however, that although this works as outlined, it does produce
     * clusters of very unequal sizes, very often with just one or two clusters
     * having most points, yet others having very, very few points. This arises
     * because as a cluster grows bigger, it becomes proportionately more likely
     * to contain the closest point to any randomly selected point. Growth rates
     * are thus proportional to sizes. The remedy employed here is to grow
     * clusters perfectly evenly by selecting the cluster of the smallest size,
     * and then finding a random point that is closer to that cluster than to
     * any other. Only if this condition fails does the routine revert to simply
     * random point selection regardless of cluster ID.
     */
    int tempi, nearest_in, count [2], this_cluster;
    double dmin, tempd;
    bool check;
    std::vector <int> pt_list;
    cluster_ids.resize (npts);

    boost::uniform_real <> uni_dist (0, 1);
    boost::variate_generator <base_generator_type&,
        boost::uniform_real <> > runif ((*generator), uni_dist);
    /* boost::normal_distribution <> norm_dist (0.0, 1.0);
       boost::variate_generator <base_generator_type&,
       boost::normal_distribution <> > rnorm ((*generator), norm_dist); */
    // Burn generator in
    for (int i=0; i<20; i++) 
        tempd = runif ();

    for (int i=0; i<npts; i++) 
        cluster_ids (i) = INT_MIN;
    count [0] = npts;
    // Set up cluster centres
    for (int i=0; i<num_clusters; i++) 
    {
        tempi = floor (runif () * npts);
        while (cluster_ids (tempi) > INT_MIN)
            tempi = floor (runif () * npts);
        cluster_ids (tempi) = i;
        count [0]--;
    }

    // The allocate the remaining points
    while (count [0] > 0) 
    {
        tempi = 9999;
        this_cluster = INT_MIN;
        for (int i=0; i<num_clusters; i++) 
        {
            count [1] = 0;
            for (int j=0; j<npts; j++)
                if (cluster_ids (j) == i) 
                    count [1]++;
            
            if (count [1] < tempi) 
            {
                tempi = count [1];
                this_cluster = i;
            }
        } // end for i
        // this_cluster then identifies the smallest cluster, with the following
        // lines finding all points that are closer to that cluster than to any
        // other.
        pt_list.resize (0);
        for (int i=0; i<npts; i++)
            if (cluster_ids (i) == INT_MIN) 
            {
                tempd = 999999.0;
                tempi = INT_MIN;
                // The search for closest point to i that is in a cluster, and
                // get cluster_num.
                for (int j=0; j<npts; j++)
                    if (cluster_ids (j) > INT_MIN && distmat (i, j) < tempd) {
                        tempd = distmat (i, j);
                        tempi = j;
                    }
                // end for j - tempi is the closest point to i that is in a
                // cluster.
                if (cluster_ids (tempi) == this_cluster)
                    pt_list.push_back (i);
            }

        if (pt_list.size () == 0) // Then just select a random point
        {
            tempi = floor (runif () * count [0]);
            count [1] = 0;
            for (int i=0; i<npts; i++) 
            {
                if (cluster_ids (i) == INT_MIN) 
                    count [1]++;
                if (count [1] == tempi) 
                {
                    tempi = count [1];
                    break;
                }
            }  // end for i

            // tempi is then the random point not in a cluster, and is directly
            // indexed into cluster_ids. The next lines find the nearest
            // cluster.
            dmin = 999999.9;
            for (int i=0; i<npts; i++)
                if (cluster_ids (i) > INT_MIN && distmat (tempi, i) < dmin) 
                {
                    dmin = distmat (tempi, i);
                    nearest_in = i;
                    this_cluster = cluster_ids (i);
                } // end if

            // Then find point closest to nearest_in that is not in a cluster.
            // This may or may not be the same as tempi above.
            dmin = 999999.9;
            for (int i=0; i<npts; i++)
                if (cluster_ids (i) == INT_MIN && 
                        distmat (nearest_in, i) < dmin) 
                {
                    dmin = distmat (nearest_in, i);
                    tempi = i;
                }
        }
        else // Pick random point from pt_list 
        {
            tempi = floor (runif () * pt_list.size ());
            tempi = pt_list [tempi];
        } // end else

        // tempi is then simply added to this_cluster
        cluster_ids (tempi) = this_cluster;
        count [0]--;
    } // end while not_in_clust.size () > 0

    pt_list.resize (0);
}; // end function allocateClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        GETNEIGHBOURS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void ClusterData::getNeighbours (Points_with_id *points)
{
    int tempi, tempj;

    Delaunay triangulation;
    triangulation.insert ((*points).begin (), (*points).end ());
    nbs.resize (npts, npts);
    for (int i=0; i<npts; i++) 
        nbs (i, i) = false;
    for (int i=0; i<(npts - 1); i++)
        for (int j=(i+1); j<npts; j++) 
            nbs (i, j) = false;

    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
            fit != triangulation.finite_faces_end(); ++fit) 
    {
        Delaunay::Face_handle face = fit;
        for (int i=0; i<2; i++) 
        {
            tempi = face->vertex(i)->info();
            for (int j=(i+1); j<3; j++) 
            {
                tempj = face -> vertex (j) -> info ();
                nbs (tempi, tempj) = true;
            } // end for j
        } // end for i
    }
} // end function getNeighbours


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                     GETCONTIGUOUSCLUSTERS                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int ClusterData::getContiguousClusters (Points_with_id *pts, 
        int nclusters, int nnew, base_generator_type *generator)
{
    int tempi, maxClusterSize = 0, clustSize;
    bool check;
    std::vector <int> clustIdsTemp;

    boost::uniform_real <> uni_dist (0, 1);
    boost::variate_generator <base_generator_type&,
        boost::uniform_real <> > runif ((*generator), uni_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        tempi = round (100.0 * runif ());
    // Then seleect random clusters
    clustIdsTemp.resize (0);
    for (int i=0; i<nnew; i++) 
    {
        tempi = floor (runif () * (double) nclusters);
        clustIdsTemp.push_back (tempi);
    }

    for (int i=0; i<(nnew - 1); i++) 
    {
        clustSize = 1; // includes i
        for (int j=(i+1); j<nnew; j++) 
        {
            check = false;
            for (int k=0; k<npts; k++)
                if (cluster_ids (k) == clustIdsTemp [i])
                    for (int m=0; m<npts; m++)
                        if (nbs (k, m) && cluster_ids (m) == clustIdsTemp [j]) 
                        {
                            if (!check) 
                            {
                                check = true;
                                clustSize++;
                            }
                            break;
                        }
        }
        if (clustSize > maxClusterSize) 
            maxClusterSize = clustSize;
    }
    clustIdsTemp.resize (0);

    return maxClusterSize;
} // end function getContiguousClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              TABLE                                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void ClusterData::getTable (ivec *cluster_ids)
{
    int maxc = 0;

    for (int i=0; i<npts; i++)
        if ((*cluster_ids) (i) > maxc)
            maxc = (*cluster_ids) (i);

    table.resize (maxc + 1);
    for (int i=0; i<maxc; i++) 
        table (i) = 0;

    for (int i=0; i<npts; i++)
        table ((*cluster_ids) (i))++;
}; // end function table


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            GETDISTS                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void ClusterData::getdists (Points_with_id *pts)
{
    double x [2], y [2];
    distmat.resize (npts, npts);

    // points can then be accessed by .first, which has [0:1] for (x,y), and
    // .second, which is the ID number
    for (int i=0; i<npts; i++) 
        distmat (i, i) = 0.0;
    for (int i=0; i<(npts - 1); i++) 
    {
        x [0] = (*pts) [i].first [0];
        y [0] = (*pts) [i].first [1];
        for (int j=(i+1); j<npts; j++) 
        {
            x [1] = (*pts) [j].first [0] - x [0];
            y [1] = (*pts) [j].first [1] - y [0];

            distmat (i, j) = sqrt (x [1] * x [1] + y [1] * y [1]);
            distmat (j, i) = distmat (i, j);
        } // end for j
    } // end for i
}
