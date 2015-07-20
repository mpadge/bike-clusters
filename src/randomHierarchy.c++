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

/*
 * Results from 10000 repeats are:
 * -----------------------------------------------------------------------
 Calculating 10000 random clusters for direction = FROM:
 [0, nc=3]: contig = 1.9988 +/- 0.0016 / observed = 2; (N, p, p_cum) = (9988, 0.9988, 0.9988); T = -3.0017
 [1, nc=8]: contig = 2.4736 +/- 0.3077 / observed = 3; (N, p, p_cum) = (5025, 0.5025, 0.5019); T = -94.8899
 [2, nc=11]: contig = 2.8367 +/- 0.5579 / observed = 4; (N, p, p_cum) = (1947, 0.1947, 0.0977); T = -155.7457
 [3, nc=15]: contig = 3.0089 +/- 0.7558 / observed = 5; (N, p, p_cum) = (518, 0.0518, 0.0051); T = -229.0283
 [4, nc=20]: contig = 1.7268 +/- 0.4073 / observed = 3; (N, p, p_cum) = (1042, 0.1042, 0.0005); T = -199.4982
 [5, nc=24]: contig = 1.2303 +/- 0.1774 / observed = 1; (N, p, p_cum) = (10000, 1.0000, 0.0005); T = 54.6738
 
 Calculating 10000 random clusters for direction = TO:
 [0, nc=4]: contig = 1.9469 +/- 0.0507 / observed = 2; (N, p, p_cum) = (9469, 0.9469, 0.9469); T = -23.5908
 [1, nc=8]: contig = 2.4476 +/- 0.3181 / observed = 3; (N, p, p_cum) = (4827, 0.4827, 0.4571); T = -97.9449
 [2, nc=11]: contig = 2.1843 +/- 0.3834 / observed = 2; (N, p, p_cum) = (8837, 0.8837, 0.4039); T = 29.7626
 [3, nc=15]: contig = 3.0018 +/- 0.7860 / observed = 4; (N, p, p_cum) = (2694, 0.2694, 0.1088); T = -112.5934
 [4, nc=20]: contig = 1.7296 +/- 0.4094 / observed = 2; (N, p, p_cum) = (6237, 0.6237, 0.0679); T = -42.2591
 [5, nc=23]: contig = 2.0413 +/- 0.5391 / observed = 2; (N, p, p_cum) = (7819, 0.7819, 0.0531); T = 5.6251
 *
 */

#include "randomHierarchy.h"

/**********************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main (int argc, char *argv []) {
    int tempi;
    double tempd, mn, sd, tt;
    std::string fname;
    std::ofstream out_file;
    clock_t timer[2];
    time_t seed;
    base_generator_type generator(42u);

    // Default direction is from, otherwise for any diri != 0, direction = to.
    int diri = atoi (argv [2]);

    //seed = atoi (argv [3]);
    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));
    timer [0] = clock();

    std::cout.setf (std::ios::fixed, std::ios::floatfield);
    std::cout.precision (4);

    // Allocate the random points with CGAL:
    const int numPoints = npts;
    // Random_points in square makes square with values in +/- first argument, and
    // size passed implicitly as numPoints.
    CGAL::Random_points_in_square_2<Point> g(1.0); // random points generator
    Points_with_id points;
    Delaunay tr;
    for (int i=0; i<numPoints; i++) 
    {
        points.push_back (std::make_pair (*g, i));
        g++;
    }
    // points can then be accessed by .first, which has [0:1] for (x,y), and
    // .second, which is the ID number
    bmat nbs = getNeighbours (&points);
    dmat distmat = getdists (&points);
    ivec cluster_ids;
    cluster_ids.resize (npts);

    // I've not yet figured out how to use boost::assign to fill a boost vector - it
    // only seems to work with std::vector!
    ivec numClusters (6), numContig (6), numTot (6);
    if (diri == 0) 
    {
        std::cout << "Calculating " << num_repeats << 
            " random clusters for direction = FROM:" << std::endl;
        numClusters (0) = 3;
        numClusters (1) = 8;
        numClusters (2) = 11;
        numClusters (3) = 15;
        numClusters (4) = 20;
        numClusters (5) = 24;
        numContig (0) = 2;
        numContig (1) = 3;
        numContig (2) = 4;
        numContig (3) = 5;
        numContig (4) = 3;
        numContig (5) = 1;
        numTot (0) = 2;
        numTot (1) = 3;
        numTot (2) = 4;
        numTot (3) = 5;
        numTot (4) = 3;
        numTot (5) = 2;
    } else {
        std::cout << "Calculating " << num_repeats << 
            " random clusters for direction = TO:" << std::endl;
        numClusters (0) = 4;
        numClusters (1) = 8;
        numClusters (2) = 11;
        numClusters (3) = 15;
        numClusters (4) = 20;
        numClusters (5) = 23;
        numContig (0) = 2;
        numContig (1) = 3;
        numContig (2) = 2;
        numContig (3) = 4;
        numContig (4) = 2;
        numContig (5) = 2;
        numTot (0) = 2;
        numTot (1) = 3;
        numTot (2) = 3;
        numTot (3) = 5;
        numTot (4) = 3;
        numTot (5) = 4;
    }

    ivec counts (numClusters.size ());
    for (int i=0; i<numClusters.size (); i++) 
    {
        mn = sd = 0.0;
        counts (i) = 0;
        for (int j=0; j<num_repeats; j++) 
        {
            int mincount = 0;
            while (mincount < 5) 
            {
                cluster_ids = allocateClusters (numClusters [i], &distmat, &generator);
                ivec counts = table (&cluster_ids);
                mincount = 0;
                for (int k=0; k<counts.size (); k++) 
                    if (counts (k) > mincount) mincount = counts (k);
            } // end while mincout < 5

            tempi = getContiguousClusters (&cluster_ids, &points, &nbs, 
                    numClusters (i), numTot (i), &generator);
            mn += (double) tempi;
            sd += (double) tempi * (double) tempi;

            if (tempi >= numContig (i)) 
                counts (i)++;
        } // end for j
        mn = mn / (double) num_repeats;
        sd = sd / ((double) num_repeats - 1.0) - mn * mn;
        tt = (mn - (double) numContig (i)) * sqrt ((double) num_repeats) / sqrt (sd);

        // Calculate cumulative probability:
        tempd = 1.0;
        for (int j=0; j<=i; j++)
            tempd = tempd * (double) counts (j) / (double) num_repeats;

        std::cout << "[" << i << ", nc=" << numClusters (i) << "]: contig = " <<
            mn << " +/- " << sd << " / observed = " << numContig (i) << 
            "; (N, p, p_cum) = (" << counts (i) << ", " <<
            (double) counts (i) / (double) num_repeats << ", " << tempd << 
            "); T = " << tt << std::endl;
    } // end for i

    /* 
     * Note that for N=1000 trials, the critical T-value is around 1.65. All
     * values are above this, so the overall p-values are really, really low.
     * The T-values can be converted to probabalities with pt, which gives the
     * area to the left of T. For one-tailed tests, 1-pt will thus give the area
     * to the right.  > 1 - pt (T, 1000)
     */

    nbs.resize (0, 0);
    distmat.resize (0, 0);

    timer[1] = clock() - timer[0];
    std::cout<<"Total Calculation Time = ";
    timeout(timer[1] / ((double)CLOCKS_PER_SEC));
    std::cout<<std::endl;

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
 **                        ALLOCATEPOINTS                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


dmat allocatePoints (base_generator_type * generator)
{
    double tempd;
    dmat xy (npts, 2);

    boost::uniform_real <> uni_dist (0, 1);
    boost::variate_generator <base_generator_type&,
        boost::uniform_real <> > runif ((*generator), uni_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        tempd = runif ();

    for (int i=0; i<npts; i++) 
    {
        xy (i, 0) = runif ();
        xy (i, 1) = runif ();
    } // end for i

    return xy;
} // end function allocatePoints


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       ALLOCATECLUSTERS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

ivec allocateClusters (int num_clusters, dmat *distmat, 
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
    ivec cluster_ids; // The vector of cluster IDs that is returned at end.
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
                    if (cluster_ids (j) > INT_MIN && (*distmat) (i, j) < tempd) {
                        tempd = (*distmat) (i, j);
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
                if (cluster_ids (i) > INT_MIN && (*distmat) (tempi, i) < dmin) 
                {
                    dmin = (*distmat) (tempi, i);
                    nearest_in = i;
                    this_cluster = cluster_ids (i);
                } // end if

            // Then find point closest to nearest_in that is not in a cluster.
            // This may or may not be the same as tempi above.
            dmin = 999999.9;
            for (int i=0; i<npts; i++)
                if (cluster_ids (i) == INT_MIN && 
                        (*distmat) (nearest_in, i) < dmin) 
                {
                    dmin = (*distmat) (nearest_in, i);
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

    return cluster_ids;
}; // end function allocateClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        GETNEIGHBOURS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

bmat getNeighbours (Points_with_id *points)
{
    int tempi, tempj;

    Delaunay triangulation;
    triangulation.insert ((*points).begin (), (*points).end ());
    bmat nbs (npts, npts);
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

    return nbs;
} // end function getNeighbours


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                     GETCONTIGUOUSCLUSTERS                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int getContiguousClusters (ivec *clIds, Points_with_id *pts, bmat *nbs, 
        int nclusters, int nnew, base_generator_type *generator)
{
    int tempi, maxClusterSize = 0, clustSize;
    bool check;
    std::vector <int> clustIds;

    boost::uniform_real <> uni_dist (0, 1);
    boost::variate_generator <base_generator_type&,
        boost::uniform_real <> > runif ((*generator), uni_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        tempi = round (100.0 * runif ());
    // Then seleect random clusters
    clustIds.resize (0);
    for (int i=0; i<nnew; i++) 
    {
        tempi = floor (runif () * (double) nclusters);
        clustIds.push_back (tempi);
    }

    for (int i=0; i<(nnew - 1); i++) 
    {
        clustSize = 1; // includes i
        for (int j=(i+1); j<nnew; j++) 
        {
            check = false;
            for (int k=0; k<npts; k++)
                if ((*clIds) (k) == clustIds [i])
                    for (int m=0; m<npts; m++)
                        if ((*nbs) (k, m) && (*clIds) (m) == clustIds [j]) 
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
    clustIds.resize (0);

    return maxClusterSize;
} // end function getContiguousClusters


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              TABLE                                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

ivec table (ivec *cluster_ids)
{
    int maxc = 0;
    ivec counts;

    for (int i=0; i<npts; i++)
        if ((*cluster_ids) (i) > maxc)
            maxc = (*cluster_ids) (i);

    counts.resize (maxc + 1);
    for (int i=0; i<maxc; i++) 
        counts (i) = 0;

    for (int i=0; i<npts; i++)
        counts ((*cluster_ids) (i))++;

    return counts;
}; // end function table


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            GETDISTS                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

dmat getdists (Points_with_id *pts)
{
    double x [2], y [2];
    dmat distmat (npts, npts);

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

    return distmat;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             TIMEOUT                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void timeout(double tseconds)
{
    int hh = floor(tseconds / 3600.0);

    if (hh == 0) 
        std::cout<<"00:";
    else if (hh < 10) 
        std::cout<<"0"<<hh<<":";
    else 
        std::cout<<hh<<":";

    double trem = tseconds - (double) hh * 3600.0;
    int mm = floor(trem / 60.0);

    if (mm == 0) 
        std::cout<<"00:";
    else if (mm < 10) 
        std::cout<<"0"<<mm<<":";
    else 
        std::cout<<mm<<":";

    double ss = trem - (double) mm * 60.0;

    if (ss == 0.0) 
        std::cout<<"00:";
    else if (ss < 10) 
        std::cout<<"0"<<ss;
    else 
        std::cout<<ss;
} // end function timeout
