/***************************************************************************
 *  Project:    BikeClusters
 *  File:       randomHierarchy.h
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


#ifndef RHIER_H
#define RHIER_H

#include "Utils.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Kernel::Point_2                                                Point;
typedef std::vector< std::pair <Point, unsigned> > Points_with_id;

const int npts = 100, num_repeats = 1000;

dmat allocatePoints (base_generator_type * generator);
ivec allocateClusters (int num_clusters, dmat *distmat, 
        base_generator_type *generator);
bmat getNeighbours (Points_with_id *pts);
int getContiguousClusters (ivec *clIds, Points_with_id *pts, bmat *nbs, 
        int nclusters, int nnew, base_generator_type *generator);
ivec table (ivec *cluster_ids);
dmat getdists (Points_with_id *pts);
void timeout (double tseconds);


#endif