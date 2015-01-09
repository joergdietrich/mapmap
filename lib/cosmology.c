/* $Id: cosmology.c,v 1.1 2006/09/06 08:05:38 cvs Exp $
 * $Date: 2006/09/06 08:05:38 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

#include <math.h>
#include <stdio.h>

#include <gsl/gsl_integration.h>

#include "cosmology.h"

/*
  This module computes various cosmological quantities like comoving,
  angular diameter, luminosity distance, lookback time etc.. Distance
  definitions are from Hogg 1999, astro-ph/9905116.
*/

/* Initialize the structure fixing the cosmological model */
void cosmo_init(cos_param *p, double m, double l, double h) {
    p->omega_m = m;
    p->omega_l = l;
    p->omega_k = 1. - m - l;
    p->h = h ;
    p->dh = 3000./h;
}


/* Hogg's function E(z); Hogg's eq. (14) */
double Ez(double z, void *p) {
    double e;
    cos_param *pars = (cos_param *)p;
    double m = pars->omega_m;
    double l = pars->omega_l;
    double k = pars->omega_k;
    
    e = sqrt(m*pow((1+z), 3) + k*pow((1+z), 2) + l);
    return e;
}


/* 1./E(z) */
double ooEz(double z, void *p) {
    return 1./Ez(z, p);
}


/* 1./(E(z) * (1+z)) */
double ooEzpz(double z, void *p) {
    return 1./(Ez(z, p)*(1+z));
}


/* Computes the line of sight comoving distance from redshifts z1 to z2. */
double dcom_los(double z1, double z2, cos_param p) {
    double dclos, err;
    size_t neval;
    gsl_function F;
    
    F.function = &ooEz;
    F.params = &p;

    if (z1>=z2) {
	fprintf(stderr, "dcom_los: z2 must be greater than z1. (z1, z2) = (%lf, %lf)\n",
		z1, z2);
	return -1;
    }
    gsl_integration_qng(&F, z1, z2, 1e-6, 1e-6, &dclos, &err, &neval);
    return p.dh * dclos;
}

/* Computes the line of transversal comoving distance (proper motion distacnce) 
   from redshifts z1 to z2. */
double dcom_tra(double z1, double z2, cos_param p) {
    double dcl, dct;

    dcl = dcom_los(z1, z2, p);
    if (p.omega_k == 0.0)
	dct = dcl;
    else if (p.omega_k > 0)
	dct = p.dh / sqrt(p.omega_k) * sinh(sqrt(p.omega_k) * dcl/p.dh);
    else
	dct = p.dh / sqrt(fabs(p.omega_k)) * sin(sqrt(fabs(p.omega_k)) * dcl/p.dh);
    return dct;
}


/* Computed the angular diameter distance between points at redshift z1 and z2 */
double dang(double z1, double z2, cos_param p) {
    double dct;

    dct = dcom_tra(z1, z2, p);
    return dct/(1.+z2);
}


/* Computes luminosity distance between points at redshifts z1 and z2.

      WARNING!                                          WARNING!  
                   This function is untested for z1>0!
      WARNING!                                          WARNING! 

*/

double dlum(double z1, double z2, cos_param p) {
    double dct;

    dct = dcom_tra(z1, z2, p);
    return (1.+z2)/(1.+z1) * dct;
}


/* Computes the comoving volume element d V_c in a solid angle d Omaga at redshift z. */
double covol(double z, cos_param p) {
    double da;

    da = dang(0, z, p);
    return p.dh * (1.+z)*(1.+z) * da*da/Ez(z, &p);
}


/* This function returns the lookback time in units of the Hubble time. */
double tlook(double z, cos_param p) {
    double tl, err;
    size_t neval;
    gsl_function F;

    F.function = &ooEzpz;
    F.params = &p;

    gsl_integration_qng(&F, 0, z, 1e-6, 1e-6, &tl, &err, &neval);
    return tl;
}
