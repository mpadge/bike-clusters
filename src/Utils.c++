/*
 * Utils.cc
 */

#ifndef UTILS_H
#define UTILS_H

#include "Utils.h"

#endif

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          CALC_ANGLE                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

double calc_angle (double x, double y)
{
    /* 
     * Calculates an angle between 0 and 360 degrees for the vector (x, y) from the
     * origin. 
     */
    double angle;

    angle = atan (y / x) * 180.0 / PI;
    if (x < 0.0) {
        angle += 180.0;
    }
    else if (y < 0.0) {
        angle += 360.0;
    }
    return angle;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           GETDISTS                                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

DistStruct getdists (double xa, double ya, double xb, double yb)
{
    /*
     * Calculates a distance in km between the two lat-lon points defined by (xa,ya)
     * and (xb,yb).  dists.dx and dists.dy are signed according to a -> b, so
     * they're positive if lat or lon of b is greater than of a; negative otherwise.
     */
    double x, y, d[2];
    DistStruct dists;

    x = (xa - xb) * PI / 180.0;
    y = (ya - yb) * PI / 180.0;
    d[0] = sin(y / 2.0) * sin(y / 2.0) + cos(ya * PI / 180.0) *
        cos(yb * PI / 180.0) * sin(x / 2.0) * sin(x / 2.0);
    d[0] = 2.0 * atan2 (sqrt (d[0]), sqrt (1.0 - d[0]));
    dists.d = d[0] * 6371.0;
    // Longitudinal distance only, for which sin(y) = y = 0:
    d[0] = cos(ya * PI / 180.0) * cos(ya * PI / 180.0) * 
        sin(x / 2.0) * sin(x / 2.0);
    d[0] = 2.0 * atan2 (sqrt (d[0]), sqrt (1.0 - d[0]));
    d[0] = d[0] * 6371.0;
    d[1] = cos(yb * PI / 180.0) * cos(yb * PI / 180.0) * 
        sin(x / 2.0) * sin(x / 2.0);
    d[1] = 2.0 * atan2 (sqrt (d[1]), sqrt (1.0 - d[1]));
    d[1] = d[1] * 6371.0;
    dists.dx = (d[0] + d[1]) / 2.0;
    // Latitudinal distance only, for which sin(x) = x = 0:
    d[0] = sin(y / 2.0) * sin(y / 2.0);
    d[0] = 2.0 * atan2 (sqrt (d[0]), sqrt (1.0 - d[0]));
    dists.dy = d[0] * 6371.0;
    // Then re-scale the 2 of them
    d[0] = sqrt(dists.dx * dists.dx + dists.dy * dists.dy);
    dists.dx = dists.dx * dists.d / d[0];
    dists.dy = dists.dy * dists.d / d[0];
    // Then put signs on them
    if (xa > xb) { dists.dx = -dists.dx;	}
    if (ya > yb) { dists.dy = -dists.dy;	}

    return dists;
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       CONVERT_DISTANCE                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

DistStruct convert_distance (double dist, double midx, double midy)
{
    /* 
     * Converts a km distance, d, to increments of latitude and longitude calculated
     * at the [lon, lat] midpt.  Output uses a dist structure for convenience, but
     * with (dx, dy) = (longitude, latitude) increments. 
     */
    //const double midx = -0.131429, midy = 51.6526, tol=1e-6;
    const double tol=1e-6;
    int count;
    double incr [3], x [2], d [2];
    DistStruct dists;

    incr [0] = 0.00001;
    incr [2] = 1.0;
    x [0] = midx - incr [0];
    x [1] = midx + incr [0];
    dists = getdists (x [0], midy, x [1], midy);
    d [0] = dists.d;
    x [0] = midx - incr [2];
    x [1] = midx + incr [2];
    dists = getdists (x [0], midy, x [1], midy);
    d [1] = dists.d;

    double err = 999999.0;
    count = 0;
    while (err > tol) {
        incr [1] = (incr [0] + incr [2]) / 2.0;
        x [0] = midx - incr [1];
        x [1] = midx + incr [1];
        dists = getdists (x [0], midy, x [1], midy);
        if (dists.d > dist && d [1] > dist) {
            incr [2] = incr [1];
        }
        else {
            incr [0] = incr [1];
        }
        err = fabs (dist - dists.d);
        count++;
        if (count > 1e6) {
            std::cout << "ERROR: Distance of " << d << 
                " unable to be converted to " <<
                " longitude increment!" << std::endl;
            break;
        }
    }
    double lon_span = 2.0 * incr [1]; // That's necessary because dist is re-used below

    incr [0] = 0.00001;
    incr [2] = 1.0;
    x [0] = midy - incr [0];
    x [1] = midy + incr [0];
    dists = getdists (midx, x [0], midx, x [1]);
    d [0] = dists.d;
    x [0] = midy - incr [2];
    x [1] = midy + incr [2];
    dists = getdists (midx, x [0], midx, x [1]);
    d [1] = dists.d;

    err = 999999.0;
    count = 0;
    while (err > tol) {
        incr [1] = (incr [0] + incr [2]) / 2.0;
        x [0] = midy - incr [1];
        x [1] = midy + incr [1];
        dists = getdists (midx, x [0], midx, x [1]);
        if (dists.d > dist && d [1] > dist) {
            incr [2] = incr [1];
        }
        else {
            incr [0] = incr [1];
        }
        err = fabs (dist - dists.d);
        count++;
        if (count > 1e6) {
            std::cout << "ERROR: Distance of " << d << 
                " unable to be converted to " <<
                " latitude increment!" << std::endl;
            break;
        }
    }
    dists.dy = 2.0 * incr [1];
    dists.dx = lon_span;

    return dists;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          REGRESSION                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

RegrResults regression(std::vector <double> x, std::vector <double> y)
{
    double sx, sx2, sy, sy2, sxy, t1, t2, xmn, ymn;
    RegrResults regr_results;

    sx = 0.0; sx2 = 0.0;
    sy = 0.0; sy2 = 0.0; sxy = 0.0;
    int count = 0, n = x.size ();
    for (int i=0; i<n; i++) {
        if (!isnan (x [i]) && !isnan (y [i])) {
            count++;
            sx += x [i];
            sx2 += x [i] * x [i];
            sy += y [i];
            sy2 += y [i] * y [i];
            sxy += x [i] * y [i];       }       }
    xmn = sx / (double) count;
    ymn = sy / (double) count;
    if (count > 0) {
        t1 = (sxy - sx * sy / (double) count);
        t2 = (sx2 - sx * sx / (double) count) * (sy2 - sy * sy / (double) count);
        regr_results.r2 = t1 / sqrt(t2); // the R-value
        regr_results.slope = t1 / (sx2 - sx * sx / (double) count); // Slope
        regr_results.intercept = sy / (double) count -
            regr_results.slope * sx / (double) count; // Intercept

        regr_results.r2 = regr_results.r2 * regr_results.r2;
        if (regr_results.slope < 0.0) { regr_results.r2 = -regr_results.r2;     }

        // Then calculate SS and tval
        sy2 = 0.0; sx2 = 0.0; count = 0;
        regr_results.SS = 0.0;
        for (int i=0; i<n; i++) {
            if (!isnan (x [i]) && !isnan (y [i])) {
                    count++;
                    t1 = regr_results.slope * x [i] + regr_results.intercept;
                    regr_results.SS += (y [i] - t1) * (y [i] - t1);
                    sx2 += (x [i] - xmn) * (x [i] - xmn);
                    sy2 += (y [i] - ymn) * (y [i] - ymn);
            }
        } // end for i
        if (count > 0) { // tval calculation
            regr_results.SS = regr_results.SS / (double) count;
            regr_results.tval = sqrt(sy2 / ((double) count - 2.0)) / sqrt(sx2);
            regr_results.tval = regr_results.slope / regr_results.tval;
        } else {
            regr_results.SS = regr_results.tval = NAN;
        }
    } // end if count > 0
    else {
        regr_results.r2 = regr_results.slope = regr_results.intercept =
            regr_results.SS = regr_results.tval = NAN;
    }
    return regr_results;
}
// end function regression



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           TIMEOUT                                  **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void timeout(double tseconds)
{
    int hh = floor(tseconds / 3600.0);
    if (hh == 0) { std::cout<<"00:";	}
    else if (hh < 10) { std::cout<<"0"<<hh<<":";	}
    else { std::cout<<hh<<":";	}
    double trem = tseconds - (double) hh * 3600.0;
    int mm = floor(trem / 60.0);
    if (mm == 0) { std::cout<<"00:";	}
    else if (mm < 10) { std::cout<<"0"<<mm<<":";	}
    else { std::cout<<mm<<":";	}
    double ss = trem - (double) mm * 60.0;
    if (ss == 0.0) { std::cout<<"00:";	}
    else if (ss < 10) { std::cout<<"0"<<ss;	}
    else { std::cout<<ss;	}
} // end function convtime


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           PROGLINE                                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void progLine (double progress)
{
    struct winsize w;
    ioctl (0, TIOCGWINSZ, &w); 
    int ncols = w.ws_col; // Number of columns in console

    int barlen = ncols - 10;
    int proglen = floor (barlen * progress);
    int gaplen = barlen - proglen;

    std::cout << "|";
    for (int i=0; i<proglen; i++) std::cout << "-";
    for (int i=0; i<gaplen; i++) std::cout << " ";
    std::cout << "| " << (int) floor (progress * 100.0) << "%\r";
    std::cout.flush ();
}
