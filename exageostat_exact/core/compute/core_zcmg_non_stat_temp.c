/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dcmg.c
 *
 * Generate covariance matrix of a set of locations in 2D using first order non-statinary Matern kernel.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2019-05-22
 *
 **/
#include "../include/exageostatcore.h"

// This function converts decimal degrees to radians
static double deg2rad(double deg)
{
    return (deg * PI / 180);
}
//  This function converts radians to decimal degrees
static double rad2deg(double rad)
{
    return (rad * 180 / PI);
}
/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
static double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d)
{
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r) / 2);
    v = sin((lon2r - lon1r) / 2);
    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

static double calculateDistance(double x1, double y1, double x2, double y2, int distance_metric)
{

    if (distance_metric == 1)
        return distanceEarth(x1, y1, x2, y2);
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}

static double calculateMahalanobisDistance(double x1, double y1, double x2, double y2, double a11, double a12, double a21, double a22)
{

    double diffx = x1 - x2;
    double diffy = y1 - y2;

    double el1 = a11 * diffx + a21 * diffy;
    double el2 = a12 * diffx + a22 * diffy;

    double ans = el1 * diffx + el2 * diffy;

    return ans;
}

/***************************************************************************/ /**
 *
 *  core_dcmg - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Matern Kernel).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.    
 *******************************************************************************
 *
 * @param[out] A
 *           The m-by-n matrix on which to compute the covariance matrix.
 *
 * @param[in] m
 *          The number of rows in the tile A. 
 *
 * @param[in] n
 *          The number of cols in the tile A. 
 *
 * @param[in] m0
 *          Global row index of the tile A.
 *
 * @param[in] n0
 *          Global col index of the tile A.
 *
 * @param[in] l1
 *          Location struct of the first input.
 *
 * @param[in] l2
 *          Location struct of the second input.
 *
 * @param[in] localtheta
 *          Parameter vector that is used to generate the output covariance matrix.
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_dcmg_non_stat_temp(double *A, int m, int n, int m0, int n0, location *l1, location *l2, double *localtheta, int distance_metric)
{

    double l1x, l1y, l2x, l2y;

    double expr = 0.0;
    double con, sigma_square, beta, nu;

    nu = (localtheta[8] + localtheta[9]) / 2;

    // Standard Deviation at s_i and s_j
    double std_dev_i = localtheta[0];
    double std_dev_j = localtheta[1];

    // Determinant of sigma_i and sigma_j;
    double det_sigma_i = localtheta[3] * localtheta[4];
    double det_sigma_j = localtheta[6] * localtheta[7];

    // Calculating elements of sigma_i
    double a_i = localtheta[3] * cos(localtheta[2]) * cos(localtheta[2]) - localtheta[4] * sin(localtheta[2]) * sin(localtheta[2]);
    double b_i = -localtheta[3] * sin(localtheta[2]) * cos(localtheta[2]) - localtheta[4] * sin(localtheta[2]) * cos(localtheta[2]);
    double c_i = localtheta[3] * sin(localtheta[2]) * cos(localtheta[2]) + localtheta[4] * sin(localtheta[2]) * cos(localtheta[2]);
    double d_i = -localtheta[3] * sin(localtheta[2]) * sin(localtheta[2]) + localtheta[4] * cos(localtheta[2]) * cos(localtheta[2]);

    // Calculating elements of sigma_j
    double a_j = localtheta[6] * cos(localtheta[5]) * cos(localtheta[5]) - localtheta[7] * sin(localtheta[5]) * sin(localtheta[5]);
    double b_j = -localtheta[6] * sin(localtheta[5]) * cos(localtheta[5]) - localtheta[7] * sin(localtheta[5]) * cos(localtheta[5]);
    double c_j = localtheta[6] * sin(localtheta[5]) * cos(localtheta[5]) + localtheta[7] * sin(localtheta[5]) * cos(localtheta[5]);
    double d_j = -localtheta[6] * sin(localtheta[5]) * sin(localtheta[5]) + localtheta[7] * cos(localtheta[5]) * cos(localtheta[5]);

    // Calculating determinant of (sigma_i + sigma_j)/2
    double det_A = (a_i + a_j) * (d_i + d_j) / 4 - (b_i + b_j) * (c_i + c_j) / 4;

    // Calculating elements of inverse of (sigma_i + sigma_j)/2 
    double a_k = (1 / det_A) * ((d_i + d_j) / 2);
    double b_k = -(1 / det_A) * ((c_i + c_j) / 2);
    double c_k = -(1 / det_A) * ((b_i + b_j) / 2);
    double d_k = (1 / det_A) * ((a_i + a_j) / 2);

    // Nugget effect
    double tau = localtheta[10];

    con = std_dev_i * std_dev_j * pow(det_sigma_i, 0.25) * pow(det_sigma_j, 0.25) * pow(det_A, -0.5);

    for (int j = 0; j < n; j++)
    {
        l1x = l1->x[j + n0];
        l1y = l1->y[j + n0];

        for (int i = 0; i < m; i++)
        {
            l2x = l2->x[i + m0];
            l2y = l2->y[i + m0];

            // Calculation of Qij as mentioned in the paper
            double Qij =  calculateMahalanobisDistance(l1x, l1y, l2x, l2y, a_k, b_k, c_k, d_k);
            

            // Need to use Mahalanobis Distance
            expr = 2 * sqrt(nu * Qij);

            if (Qij == 0)
                A[i + j * m] = con * pow(expr, nu) * gsl_sf_bessel_Knu(nu, expr) + tau; // Need to add the first term of the indicator function
            else
                A[i + j * m] = con * pow(expr, nu) * gsl_sf_bessel_Knu(nu, expr);
        }
    }
}

void core_sdcmg_nono_stat(double *A, int m, int n, int m0, int n0, location *l1, location *l2, location *lm, double *localtheta, int distance_metric)
{

    int i, j;
    double l1x, l1y, l2x, l2y, lmx, lmy;
    double theta_0i, theta_0j, theta_1i, theta_1j, theta_2i, theta_2j;
    double diffx, diffy;
    double expr = 0.0;
    double con, sigma_square, beta, nu;

    for (j = 0; j < n; j++)
    {
        l1x = l1->x[j + n0];
        l1y = l1->y[j + n0];
        diffx = abs(l1x - lmx);
        diffy = abs(l1y - lmy);
        theta_0i = localtheta[0] + (localtheta[1] * diffx) + (localtheta[2] * diffy);
        theta_1i = localtheta[3] + (localtheta[4] * diffx) + (localtheta[5] * diffy);
        theta_2i = localtheta[6] + (localtheta[7] * diffx) + (localtheta[8] * diffy);

        for (i = 0; i < m; i++)
        {
            l2x = l2->x[i + m0];
            l2y = l2->y[i + m0];
            diffx = abs(l2x - lmx);
            diffy = abs(l2y - lmy);
            theta_0j = localtheta[0] + (localtheta[1] * diffx) + (localtheta[2] * diffy);
            theta_1j = localtheta[3] + (localtheta[4] * diffx) + (localtheta[5] * diffy);
            theta_2j = localtheta[6] + (localtheta[7] * diffx) + (localtheta[8] * diffy);
            sigma_square = (theta_0i + theta_0j) / 2;
            beta = (theta_1i + theta_1j) / 2;
            nu = (theta_2i + theta_2j) / 2;
            con = pow(2, (nu - 1)) * tgamma(nu);
            con = 1.0 / con;
            con = sigma_square * con;
            //MLE calculation
            expr = calculateDistance(l1x, l1y, l2x, l2y, distance_metric) / beta;
            if (expr == 0)
                A[i + j * m] = sigma_square; /* + 1e-4*/
            else
                A[i + j * m] = con * pow(expr, nu) * gsl_sf_bessel_Knu(localtheta[2], expr); // Matern Function
        }
    }
}

void core_scmg_nono_stat(float *A, int m, int n, int m0, int n0, location *l1, location *l2, location *lm, double *localtheta, int distance_metric)
{

    int i, j;
    double l1x, l1y, l2x, l2y, lmx, lmy;
    double theta_0i, theta_0j, theta_1i, theta_1j, theta_2i, theta_2j;
    double diffx, diffy;
    double expr = 0.0;
    double con, sigma_square, beta, nu;

    for (j = 0; j < n; j++)
    {
        l1x = l1->x[j + n0];
        l1y = l1->y[j + n0];
        diffx = abs(l1x - lmx);
        diffy = abs(l1y - lmy);
        theta_0i = localtheta[0] + (localtheta[1] * diffx) + (localtheta[2] * diffy);
        theta_1i = localtheta[3] + (localtheta[4] * diffx) + (localtheta[5] * diffy);
        theta_2i = localtheta[6] + (localtheta[7] * diffx) + (localtheta[8] * diffy);

        for (i = 0; i < m; i++)
        {
            l2x = l2->x[i + m0];
            l2y = l2->y[i + m0];
            diffx = abs(l2x - lmx);
            diffy = abs(l2y - lmy);
            theta_0j = localtheta[0] + (localtheta[1] * diffx) + (localtheta[2] * diffy);
            theta_1j = localtheta[3] + (localtheta[4] * diffx) + (localtheta[5] * diffy);
            theta_2j = localtheta[6] + (localtheta[7] * diffx) + (localtheta[8] * diffy);
            sigma_square = (theta_0i + theta_0j) / 2;
            beta = (theta_1i + theta_1j) / 2;
            nu = (theta_2i + theta_2j) / 2;
            con = pow(2, (nu - 1)) * tgamma(nu);
            con = 1.0 / con;
            con = sigma_square * con;
            //MLE calculation
            expr = calculateDistance(l1x, l1y, l2x, l2y, distance_metric) / beta;
            if (expr == 0)
                A[i + j * m] = (float)sigma_square; /* + 1e-4*/
            else
                A[i + j * m] = (float)(con * pow(expr, nu) * gsl_sf_bessel_Knu(localtheta[2], expr)); // Matern Function
        }
    }
}
