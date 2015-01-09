#ifndef _PREFS_H_
#define _PREFS_H_ 1

#include <limits.h>

#include <fitsio.h>

typedef struct 
{
    double h;        /* scaled Hubble constant H_0 = 100 km/s/Mpc */
    double Omega_M;  /* Matter density */
    double Omega_L;  /* Cosmological constant */
    double spix;     /* Pixel size */
    double zcl;      /* Fiducial cluster redshift */
    char   filter[PATH_MAX]; /* filter function type */ 
    double fradius;  /* filter radius in pixels */
    double fradius_orig; /* Keep the configured value in this variable */
    char   ffile[PATH_MAX]; /* File with filter function to read 
			     * if filter== FILE */
    double ovrsht;   /* Overshooting of the filter scale */
    double xc;       /* Mischa's xc */
    char   tblname[FLEN_KEYWORD];      /* Fits table name of objects table */
    char   xname[FLEN_KEYWORD];        /* Column name for right ascencion */
    char   yname[FLEN_KEYWORD];        /* Column name for declination */
    char   e1name[FLEN_KEYWORD];       /* Column name for e1 component */
    char   e2name[FLEN_KEYWORD];       /* Column name for e2 component */
    int    doweight;                   /* Use weighting of ellipticities */
    char   weightname[FLEN_KEYWORD];   /* Column name for magnitude */ 
    int    dotomo;                     /* Perform lens tomography */
    char   zname[FLEN_KEYWORD];        /* Column name for redshift */
    char   outname[PATH_MAX];          /* Output name */
    char   mapsuf[PATH_MAX];           /* Suffix for Map maps */
    int    smap;                       /* Create significance maps */
    char   smapsuf[PATH_MAX];          /* Suffix for smaps */
    int    invert;                     /* invert images */
    char   invsuf[PATH_MAX];           /* suffix for inverted images */
    int    cross;                      /* rotate ellipticities by 45 deg */
    char   xsuf[PATH_MAX];             /* Suffix for cross Map maps */
    int    mock;                       /* Create mock field */
    char   mocksuf[PATH_MAX];          /* Suffix for mock maps */
    double zmin;                       /* Minimum redshift for tomography */
    double zmax;                       /* Maximum redshift for tomography */
    int    nplanes;                    /* Number of lens planes */ 
    int    verbose;         /* level of verbosity: 0 quiet
			     *                     1 normal output
			     *                     2 debug */
    int    nthreads;        /* Number of threads if compiled with OpenMP */
    char   incat[PATH_MAX]; /* Path to the input catalog */
} prefstruct;

prefstruct prefs;

/* Function in prefs.c */
void dumpprefs(void);
void readprefs(char *, char **, char **, int);
int findkeys(char *, char keyw[][16]);

#endif
