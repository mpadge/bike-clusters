/*
 * Utils.h
 */

#include <stdlib.h> // has abs function
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <iomanip> // for setfill
#include <sys/ioctl.h> // for console width: Linux only!

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

#define PI 3.1415926535897932384626433832795

typedef boost::numeric::ublas::vector <int> ivec;
typedef boost::numeric::ublas::matrix <int> imat;
typedef boost::numeric::ublas::vector <double> dvec;
typedef boost::numeric::ublas::matrix <double> dmat;
typedef boost::numeric::ublas::vector <bool> bvec;
typedef boost::numeric::ublas::matrix <bool> bmat;
typedef boost::numeric::ublas::zero_matrix <double> zmat_d;

const double DOUBLE_MAX = std::numeric_limits<double>::max (),
    DOUBLE_MIN = std::numeric_limits <double>::min ();

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;

struct myTime{
    int hh, mm;
    float ss;	};

struct DistStruct{
    double dx, dy, d;	};


struct RegrResults {
    double r2, slope, intercept, SS, tval;      };

double calc_angle (double x, double y);
DistStruct getdists (double xa, double ya, double xb, double yb);
DistStruct convert_distance (double dist, double midx, double midy);
void timeout (double tseconds);
RegrResults regression (std::vector <double> x, std::vector <double> y);
void progLine (double progress);
