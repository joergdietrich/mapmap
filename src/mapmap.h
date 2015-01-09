/* $Id: mapmap.h,v 1.9 2007/08/22 11:37:43 jdietric Exp $
 * $Date: 2007/08/22 11:37:43 $
 * $Author: jdietric $
 * $Revision: 1.9 $
 */

#ifndef _MAPMAP_H_
#define _MAPMAP_H_ 1

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fitsio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../config.h"

#include <callocmatrix.h>
#include <getaline.h>
#include <cosmology.h>
#ifdef USING_PLPA
#include <plpa.h>
#endif

#include "define.h"
#include "types.h"
#include "prefs.h"
#include "misc.h"
#include "filter.h"
#include "message.h"



#define ZMAX (5.0)
#define ZSTEPS (500)

double zindex[ZSTEPS];
double ds[ZSTEPS];
double zweightarr[ZSTEPS];


/* A global variable */
cos_param cospar;  /* Cosmological parameters */


/* functions from mapmap.c */

long get_nrows(fitsfile *);
void readfitscol(fitsfile *, char *, unsigned long, double *);
void make_out_names(char *, char *, char *, char *);
void sky2xy(double **, double **, long, mapstruct *, wcsstruct *);
void calcmap(mapstruct, mapstruct, double *, double *, double *, double *,
	     double *, int *, long, int);
double z_weight(int);
void fill_zweight_array(double);
void fill_ds_array(void);
void cross(double *, double *, unsigned long);
void randomize(double *, double *, unsigned long);
double get_z_ang(double, double);
int writefitsmap(mapstruct, char *, char *);


#endif
