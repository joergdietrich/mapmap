#ifndef _FILTER_H
#define _FILTER_H 1

static double sqrarg __attribute__ ((unused));
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


#define G (6.67259e-11)      /* Gravitational constant m^3 kg^-1 s^-2 */
#define MSUN (2e30)          /* Solar mass in kg */
#define MPC (3.1e22)         /* Megaparsec in m */
#define C (299792458)        /* Speed of light in m/s */


#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <callocmatrix.h>
#include <cosmology.h>
#include <getaline.h>

#include "define.h"
#include "prefs.h"


typedef struct {
    float c;  /* concentration parameter */
    float m200; /* Mass inside r_200 */
} nfw_struct;


typedef struct {
    double xc;
    double r;
    double c;
} weight_param;


/* skip linear white space */
#define SKIPLWS(p) while (*(p) && isspace((unsigned char) *(p))) { (p)++; }


extern cos_param cospar;  /* Cosmological parameters */

double get_rhoc(double, cos_param );
double nfw_shear(nfw_struct, double, double, double, cos_param);
double get_Sigmac(double, double, cos_param);
double nfw_g(double);
double nfw_kappa(nfw_struct, double, double, double, cos_param);

double weight(double , void *);
double schneider_weight(double);
double schirmer_weight(double, double);
double file_weight(double);
void linint(double *, double *, unsigned long, double, double *, int *);



#endif
