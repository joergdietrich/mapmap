/* $Id: cosmology.h,v 1.1 2006/09/06 08:05:50 cvs Exp $
 * $Date: 2006/09/06 08:05:50 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

#ifndef _COSMOLOGY_H 
#define _COSMOLOGY_H 1

typedef struct {
    float omega_m;
    float omega_l;
    float omega_k;
    float h;
    float dh;
} cos_param;


void cosmo_init(cos_param *p, double m, double l, double h);
double Ez(double z, void *p);
double ooEz(double z, void *p);
double ooEzpz(double z, void *p);
double dcom_los(double z1, double z2, cos_param p);
double dcom_tra(double z1, double z2, cos_param p);
double dang(double z1, double z2, cos_param p);
double dlum(double z1, double z2, cos_param p);
double covol(double z, cos_param p);
double tlook(double z, cos_param p);


#endif
