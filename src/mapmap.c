#include "mapmap.h"

#include <unistd.h>

int main(int argc, char *argv[]) {
    int a, narg, opt;
    int i;
    int status = 0;
    long nrows;
    unsigned long in;
    char outname1[PATH_MAX];
    char outname2[PATH_MAX];
    char outname3[PATH_MAX];
    char outname4[PATH_MAX];
    char *history;
    double *xg1, *xg2, *e1, *e2, *sigmae, *z;
    int *zbin = NULL;
    mapstruct map;
    mapstruct smap;
    wcsstruct newwcs;
    static char prefsname[PATH_MAX];
    char fits_err[FITS_ERR_LEN];
    char **argkey, **argval;
    fitsfile *fptr;

    msg_init(MSG_INFO, NULL);
    if(argc < 2) {
        msg_message(MSG_INFO, 0, "%s-%s", BANNER, VERSION);
	msg_message(MSG_INFO, 0, "Written and copyright by %s", COPYRIGHT);
	msg_message(MSG_FATAL, 0, "%s", "SYNTAX:");
	msg_message(MSG_FATAL, 0, "%s", SYNTAX);
	exit(EXIT_SUCCESS);
    }

    argkey = critcalloc(argc, sizeof(char *), "Could not allocate argkey: ");
    argval = critcalloc(argc, sizeof(char *), "Could not allocate argval: ");

    /* default config file */
    strcpy(prefsname, "mapmap_default.conf");

    /* This is all preferences and option reading taking from SWarp */
    narg = 0;
    for(a=1; a<argc; a++) {
        if (*(argv[a]) == '-') {
            opt = (int)argv[a][1];
            if (strlen(argv[a]) < 3 || opt == '-') {
                if (opt == '-')
                    opt = (int)tolower((int)argv[a][2]);
                switch(opt) {
                case 'c':
                    if (a < (argc-1))
                        strcpy(prefsname, argv[++a]);
                    break;
                case 'd':
                    dumpprefs();
                    exit(EXIT_SUCCESS);
                    break;
                case 'h':
                default:
                    msg_message(MSG_FATAL, 0, "SYNTAX:\n%s", SYNTAX);
                }
            } else {
                argkey[narg] = &argv[a][1];
                argval[narg++] = argv[++a];
            }
        } else {
            snprintf(prefs.incat, PATH_MAX, "%s",  argv[a]);
        }
    }
    
    readprefs(prefsname, argkey, argval, narg);

    msg_init(prefs.verbose, NULL);
    msg_message(MSG_INFO, 0, "%s-%s", BANNER, VERSION);
    msg_message(MSG_INFO, 0, "by %s", COPYRIGHT);
    msg_message(MSG_DEBUG, 0, "Verbosity level: %d", prefs.verbose);

    history = malloc(sizeof(char));
    *history = '\0';
    for(i=0; i<argc; i++) {
	history = realloc(history, 
			  (strlen(history)+strlen(argv[i])+2)*sizeof(char));
	strcat(history, argv[i]);
	strcat(history, " ");
    }

#ifdef _OPENMP
    if (prefs.nthreads)
	omp_set_num_threads(prefs.nthreads);
    prefs.nthreads = omp_get_max_threads();
    msg_message(MSG_NOTICE, 0, "Running OpenMP version with %d threads.", 
		prefs.nthreads);
#else
    msg_message(MSG_NOTICE, 0, "This build of mapmap is single-threaded.");
    prefs.nthreads = 1;
#endif
#ifdef USING_PLPA
    msg_message(MSG_NOTICE, 0, "This version is PLPA enabled.");
#endif
    cosmo_init(&cospar, prefs.Omega_M, prefs.Omega_L, prefs.h);
    
    if (prefs.cross && prefs.mock)
	msg_message(MSG_WARNING, 0, "CROSS and MOCK option are set. This "\
		    "is probably not what you want.");

    if (prefs.cross && prefs.invert)
	msg_message(MSG_WARNING, 0, "Do you expect rotation by -45 deg to be" \
		    " different from rotation by \n  45 deg?");
    
    if (prefs.dotomo && prefs.smap) {
	msg_message(MSG_FATAL, 0, "Tomography implies that no significance "\
		    "map is produced.");
	exit(EXIT_FAILURE);
    }

    fits_open_file(&fptr, prefs.incat, READONLY, &status);
    if (status) {
	snprintf(fits_err, FITS_ERR_LEN, 
		 "fits_open_file: error opening file: %s", prefs.incat);
	fits_error(fits_err, status);
    }
    nrows = get_nrows(fptr);
    msg_message(MSG_DEBUG, 0, "%s contains %ld objects", prefs.incat, nrows);

    xg1 = critcalloc(nrows, sizeof(double), "Could not allocate xg1: ");
    readfitscol(fptr, prefs.xname, nrows, xg1);

    xg2 = critcalloc(nrows, sizeof(double), "Could not allocate xg2: ");
    readfitscol(fptr, prefs.yname, nrows, xg2);

    e1 = critcalloc(nrows, sizeof(double), "Could not allocate e1: ");
    readfitscol(fptr, prefs.e1name, nrows, e1);

    e2 = critcalloc(nrows, sizeof(double), "Could not allocate e2: ");
    readfitscol(fptr, prefs.e2name, nrows, e2);

    if (prefs.doweight) {
	sigmae = critcalloc(nrows, sizeof(double), 
			    "Could not allocate sigmae: ");
	readfitscol(fptr, prefs.weightname, nrows, sigmae);
    } else
	/* shut compiler warnings up */
	sigmae = NULL;
    
    if (prefs.dotomo) {
	z = critcalloc(nrows, sizeof(double),
		       "Could not allocate z: ");
	readfitscol(fptr, prefs.zname, nrows, z);
	zbin = critcalloc(nrows, sizeof(int), "Could not allocate zbin: ");
	/* instead of using the exact redshifts, we use redshift bins
	   of size ZSTEP (defined at compile time). Here we create an
	   array in which we record for every galaxy which bin number
	   it is in. */
	for (i=0; i<nrows; i++) 
	    zbin[i] = (z[i]+1./(2.*ZSTEPS))/(ZMAX/ZSTEPS);
    }
    else
	z = 0;

    fits_close_file(fptr, &status);
    fptr = NULL;
    msg_message(MSG_INFO, 0, "Read %ld objects", nrows);

    sky2xy(&xg1, &xg2, nrows, &map, &newwcs);
    map.wcs = &newwcs;
    smap.naxes[0] = map.naxes[0];
    smap.naxes[1] = map.naxes[1];
    smap.naxes[2] = map.naxes[2];
    smap.wcs = &newwcs;

    make_out_names(outname1, outname2, outname3, outname4);

    if (prefs.cross)
	cross(e1, e2, nrows);

    if (prefs.mock) 
	randomize(e1, e2, nrows);

    if (!prefs.dotomo)
	prefs.nplanes = 1;

    map.arr = critcalloc(map.naxes[0]*map.naxes[1]*prefs.nplanes, 
			 sizeof(double), "Could not allocate map array: ");
    smap.arr = critcalloc(smap.naxes[0]*smap.naxes[1]*prefs.nplanes, 
			  sizeof(double), "Could not allocate smap array: ");
    
    if (prefs.dotomo)
	fill_ds_array();

    for(i=0; i<prefs.nplanes; i++) {
	if (prefs.dotomo)
	    msg_message(MSG_NOTICE, 0, "Lens plane %d/%d", i+1, 
			prefs.nplanes);
	calcmap(map, smap, xg1, xg2, e1, e2, sigmae, zbin, nrows, i);
    }
    
    writefitsmap(map, outname1, history);
    if (prefs.smap)
	writefitsmap(smap, outname2, history);

    if (prefs.invert) {
	msg_message(MSG_DEBUG, 0, "Inverting galaxy orientation");
	for (in=0; in < (unsigned)(map.naxes[0]*map.naxes[1]); in++) {
	    map.arr[in] = -map.arr[in];
	    if (prefs.smap)
		smap.arr[in] = -smap.arr[in];
	}
	msg_message(MSG_DEBUG, 0, "done.");
	writefitsmap(map, outname3, history);
	if (prefs.smap)
	    writefitsmap(smap, outname4, history);
    }

    msg_finish();
    return EXIT_SUCCESS;
}


void make_out_names(char *outname1, char *outname2, char *outname3, 
		    char *outname4 ) {
    int outlen1 =0;
    int outlen2 = 0;
    int outlen3 =0;
    int outlen4 = 0;

    strncpy(outname1, prefs.outname, PATH_MAX);
    outlen1 = strlen(prefs.outname);
    strncat(outname1, prefs.mapsuf, PATH_MAX-outlen1); 
    outlen1 += strlen(prefs.outname);
    strncpy(outname2, prefs.outname, PATH_MAX);
    outlen2 = strlen(prefs.outname);
    strncat(outname2, prefs.smapsuf, PATH_MAX-outlen1); 
    outlen2 += strlen(prefs.smapsuf);

    
    if (prefs.cross) {
	strncat(outname1, prefs.xsuf, PATH_MAX-outlen1);
	outlen1 += strlen(prefs.xsuf);
	strncat(outname2, prefs.xsuf, PATH_MAX-outlen2);
	outlen2 += strlen(prefs.xsuf);
    }

    if (prefs.mock) {
	strncat(outname1, prefs.mocksuf, PATH_MAX-outlen1);
	outlen1 += strlen(prefs.mocksuf);
	strncat(outname2, prefs.mocksuf, PATH_MAX-outlen2);
	outlen2 += strlen(prefs.mocksuf);
    }

    if (prefs.invert) {
	outlen3 = outlen1;
	outlen4 = outlen2;
	strncpy(outname3, outname1, PATH_MAX);
	strncpy(outname4, outname2, PATH_MAX);
	strncat(outname3, prefs.invsuf, PATH_MAX-outlen3);
	outlen3 += strlen(prefs.invsuf);
	strncat(outname4, prefs.invsuf, PATH_MAX-outlen3);
	outlen4 += strlen(prefs.invsuf);
	strncat(outname3, ".fits", PATH_MAX-outlen3);
	strncat(outname4, ".fits", PATH_MAX-outlen4);
    }

    strncat(outname1, ".fits", PATH_MAX-outlen1);
    strncat(outname2, ".fits", PATH_MAX-outlen2);
    return;
}


void sky2xy(double **ra, double **dec, long nrows, mapstruct *map, 
	    wcsstruct *wcs) {

    double minra, maxra, cra;
    double mindec, maxdec, cdec;
    double delta, cdelt, sizera, sizedec;
    long i;
    
    min_max_center(*ra, nrows, &minra, &maxra, &cra);
    min_max_center(*dec, nrows, &mindec, &maxdec, &cdec);
    sizera = (maxra - minra) * cos(M_PI/180. * cdec);
    sizedec = (maxdec-mindec);
    delta = get_z_ang(prefs.spix, prefs.zcl);
    cdelt = 180./M_PI * delta;
    msg_message(MSG_DEBUG, 0, "Min Ra: %lf\tMax Ra: %lf\tCen Ra: %lf", 
		minra, maxra, cra);
    msg_message(MSG_DEBUG, 0, "Min Dec: %lf\tMax Dec: %lf\tCen Dec: %lf", 
		mindec, maxdec, cdec);
    msg_message(MSG_DEBUG, 0, "CDELT: %lf", cdelt);
	
    map->naxes[0] = (long)(ceil(sizera/cdelt));
    map->naxes[1] = (long)(ceil(sizedec/cdelt));
    if (prefs.dotomo)
	map->naxes[2] = prefs.nplanes;
    else
	map->naxes[2] = 1;
    if (prefs.verbose) {
	if (prefs.dotomo)
	    msg_message(MSG_INFO, 0, "Output dimensions: %ld x %ld x %ld pix^3", 
			map->naxes[0], map->naxes[1], map->naxes[2]);
	else
	    msg_message(MSG_INFO, 0, "Output dimensions: %ld x %ld pix^2", 
			map->naxes[0], map->naxes[1]);
    }
    for (i=0; i<nrows; i++) {
	(*ra)[i] = (cra - (*ra)[i]) * cos(cdec*M_PI/180.)/cdelt 
	    + sizera/(2.*cdelt);
	(*dec)[i] = ((*dec)[i] - cdec)/cdelt + sizedec/(2.*cdelt);
    }
    wcs->cdelt1 = -cdelt;
    wcs->cdelt2 = cdelt;
    wcs->crval1 = cra;
    wcs->crval2 = cdec;
    wcs->crpix1 = map->naxes[0]/2.;
    wcs->crpix2 = map->naxes[1]/2.;
    wcs->ctype1 = strdup("RA---TAN");
    wcs->ctype2 = strdup("DEC--TAN");
    wcs->equinox = 2000;
    if (prefs.dotomo) {
	wcs->cdelt3 = (prefs.zmax-prefs.zmin)/(prefs.nplanes-1);
	wcs->crval3 = prefs.zmin;
	wcs->crpix3 = 1.0;
	wcs->ctype3 = strdup("LINEAR");
    }
    
    /* Filter is in kpc */
    prefs.fradius_orig = prefs.fradius;
    prefs.fradius = 180./M_PI * get_z_ang(prefs.fradius, prefs.zcl);
    /* Now it is in degrees */
    msg_message(MSG_INFO, 0, "Filter scale corresponds to %2.2lf arcmin.",
		60*prefs.fradius);
    prefs.fradius /= cdelt;
    /* and now in pixels */
    msg_message(MSG_DEBUG, 0, "This is %2.2lf pixels.", prefs.fradius);
}


/* Get the number of rows in a table */
long get_nrows(fitsfile *fptr) {
    long nrows;
    int status = 0;
    char fits_err[FITS_ERR_LEN];

    fits_movnam_hdu(fptr, BINARY_TBL, prefs.tblname, 0, &status);
    if (status == BAD_HDU_NUM) { /* try different table type */
	status = 0;
	fits_movnam_hdu(fptr, ASCII_TBL, prefs.tblname, 0, &status);
    }
    if (status) {
	snprintf(fits_err, 128, "fits_movnam_hdu: no such table: %s\n", 
		 prefs.tblname);
	fits_error(fits_err, status);
    }
    fits_get_num_rows(fptr, &nrows, &status);
    if (status)
	fits_error("fits_get_num_rows: could not determine number of entries in table\n", status);

    return nrows;
}


void calcmap(mapstruct map, mapstruct smap, 
	     double *xg1, double *xg2, /* x1, x2 coordinates of galaxies (pixels) */
	     double *e1, double *e2, /* ellipticity components */
	     double *sigmae,         /* Sigma ellipticity for weighting */
	     int *zbin,              /* redshift */
	     long ng,                /* number of galaxies */
	     int lp  /* Lens plane we are currently calculating*/
	     ) {

    double rsq;
    double dx, dy, dxsq, dysq, x, xsq;
    double et, wgt, ewgtsq, e1sq, e2sq;
    double normfac, err;
    double zd, zinc; /* redshift of deflector (lens plane), Delta z between 
		        planes */
    double zwgt;
    double **myxg1, **myxg2, **mye1, **mye2, **mysigmae;
    int **myzbin;
    double **mymap, **mysmap;
    long ixmin, ixmax, iymin, iymax;
    long ix, iy;
    long i, in, myi;
    int iam;
    gsl_function F;
    gsl_integration_workspace *ws;
    weight_param p;

    if (strcmp(prefs.filter, "FILE"))
	p.xc = prefs.xc;
    else { 
	p.r = prefs.fradius;
	p.c = map.wcs->cdelt2;
    }

    F.function = &weight;
    F.params = &p;
    ws = gsl_integration_workspace_alloc(8192);
    gsl_integration_qag(&F, 0, prefs.ovrsht, 1e-6, 1e-6, 8192, 
			GSL_INTEG_GAUSS41, ws, &normfac, &err);
    gsl_integration_workspace_free(ws);

    if (prefs.dotomo) {
	zinc = (prefs.zmax-prefs.zmin)/(prefs.nplanes-1);
	zd = prefs.zmin + lp * zinc;
	msg_message(MSG_INFO, 0, "redshift of lens plane: %.2f", zd);
	fill_zweight_array(zd);
    }
    rsq = SQR(prefs.fradius);

    mymap = critcalloc(prefs.nthreads, sizeof(double *), "mymap");
    mysmap = critcalloc(prefs.nthreads, sizeof(double *), "mysmap");
    myxg1 = critcalloc(prefs.nthreads, sizeof(double *), "mysmap");
    myxg2 = critcalloc(prefs.nthreads, sizeof(double *), "mysmap");
    mye1 = critcalloc(prefs.nthreads, sizeof(double *), "mysmap");
    mye2 = critcalloc(prefs.nthreads, sizeof(double *), "mysmap");
    if (prefs.doweight)
	mysigmae = critcalloc(prefs.nthreads, sizeof(double *), "mysigmae");
    if (prefs.dotomo) 
	myzbin = critcalloc(prefs.nthreads, sizeof(int *), "myzbin");


#pragma omp parallel for schedule(static,1024) default(none) shared(ng, prefs, rsq, xg1, xg2, e1, e2, sigmae, zbin, map, smap, stdout, lp, normfac, p, mymap, mysmap, myxg1, myxg2, mye1, mye2, mysigmae, myzbin, gsl_interp_linear, zd) private(ix, iy, ixmax, ixmin, iymax, iymin, in, dx, dxsq, xsq, dy, dysq, x, e1sq, e2sq, wgt, et, ewgtsq, marg1, marg2, sqrarg, zwgt, iam, myi)
    for(i=0; i<ng; i++) {
	ixmax = (int)ceil(MIN(xg1[i]+prefs.ovrsht*prefs.fradius, map.naxes[0]));
	ixmin = (int)floor(MAX(xg1[i]-prefs.ovrsht*prefs.fradius, 0));
	iymax = (int)ceil(MIN(xg2[i]+prefs.ovrsht*prefs.fradius, map.naxes[1]));
	iymin = (int)floor(MAX(xg2[i]-prefs.ovrsht*prefs.fradius, 0));
#ifdef _OPENMP
	iam = omp_get_thread_num();
#else
	iam = 0;
#endif

	    
	/* Here first bind each thread to one CPU and then we allocate
	   all the per thread memory and touch it first (calloc!).
	 */

	if (!mymap[iam]) {
	    myi = 0;
#ifdef USING_PLPA
#pragma omp critical
	    {
		mplpa_api_type_t p;
		mplpa_cpu_set_t msk;
		int res;
		if (mplpa_api_probe(&p) != 0 ||  p != MPLPA_PROBE_OK) {
		    msg_message(MSG_ERROR, 0, "PLPA failed in thread %d.", iam);
		}
		PLPA_CPU_ZERO(&msk);
		PLPA_CPU_SET(iam, &msk);
		res = mplpa_sched_setaffinity((pid_t)0, (size_t)32, &msk);
		switch(res) {
		case ENOSYS:
		    msg_message(MSG_ERROR, 0, 
				"Thread %d failed to set processor affinity (not supported on this system).", iam);
		    break;
		case EINVAL:
		    msg_message(MSG_ERROR, 0, 
				"Thread %d failed to set processor affinity (other reasons).", iam);
		    break;
		default:
		    break;
		}
	    }
#endif

	    mymap[iam] = critcalloc(map.naxes[0]*map.naxes[1], sizeof(double), 
				    "Could not allocate mymap array: ");
	    mysmap[iam] = critcalloc(map.naxes[0]*map.naxes[1], 
				     sizeof(double), 
				     "Could not allocate mysmap array: ");
	    myxg1[iam] = critcalloc(ng, sizeof(double), "myxg1: ");
	    memcpy(myxg1[iam], xg1, sizeof(double)*ng);
	    myxg2[iam] = critcalloc(ng, sizeof(double), "myxg2: ");
	    memcpy(myxg2[iam], xg2, sizeof(double)*ng);
	    mye1[iam] = critcalloc(ng, sizeof(double), "mye1: ");
	    memcpy(mye1[iam], e1, sizeof(double)*ng);
	    mye2[iam] = critcalloc(ng, sizeof(double), "mye2: ");
	    memcpy(mye2[iam], e2, sizeof(double)*ng);
	    if (prefs.doweight) {
		mysigmae[iam] = critcalloc(ng, sizeof(double), "mysigmae: ");
		memcpy(mysigmae[iam], sigmae, sizeof(double)*ng);
	    }
	    if (prefs.dotomo) {
		myzbin[iam] = critcalloc(ng, sizeof(int), "myzbin: ");
		memcpy(myzbin[iam], zbin, sizeof(int)*ng);
	    }
	}
	
	if (prefs.dotomo) {
	    zwgt = z_weight(myzbin[iam][i]);
	}

	for(ix=ixmin; ix<ixmax; ix++) {
	    dx = xg1[i] - (double)ix;
	    dxsq = SQR(dx); 
	    for(iy=iymin; iy<iymax; iy++) {
		dy = xg2[i] - (double)iy;
		dysq = SQR(dy); 
		xsq = dxsq + dysq;
		if (xsq <= prefs.ovrsht*rsq) {
		    x = sqrt(xsq/rsq);
		    et = -(e1[i] * (dxsq-dysq) + 2. * e2[i]*dx*dy)/xsq;
		    in = ix + iy*map.naxes[0];
		    if (strcmp(prefs.filter, "FILE"))
			wgt = weight(x, &p)/normfac;
		    else 
			wgt = weight(x, &p)/normfac;

		    e1sq = SQR(e1[i]); 
		    e2sq = SQR(e2[i]); 
		    if (prefs.doweight) 
			ewgtsq = SQR(sigmae[i]); 
		    else 
			ewgtsq = 1;
		    if (prefs.dotomo) {
			mymap[iam][in] += et * wgt * zwgt;
			mysmap[iam][in] += SQR(zwgt);
		    }
		    else {
			mymap[iam][in] += et * wgt;
			mysmap[iam][in] += (e1sq+e2sq) * ewgtsq * SQR(wgt);
		    }
		}
	    }
	}
	myi++;
#ifndef _OPENMP
	if (!(i%100)) 
	    msg_message(MSG_INFO, 1, "Processed %ld out of %ld galaxies", 
			i, ng);
#else 
	if (!((myi)%100) && omp_get_thread_num() == 0)
	    msg_message(MSG_DEBUG, 1, "Thread 0 processed %ld/%ld objects", 
			myi, ng/prefs.nthreads);
#endif
    }
    for (iam=0; iam < prefs.nthreads; iam++) {
	for (ix=0; ix < map.naxes[0]; ix++) {
	    for(iy=0; iy < map.naxes[1]; iy++) {
		in = ix + iy*map.naxes[0];
		map.arr[in+lp*map.naxes[0]*map.naxes[1]] += mymap[iam][in];
		smap.arr[in+lp*smap.naxes[0]*smap.naxes[1]] += mysmap[iam][in];
	    }
	}
	free(mymap[iam]);
	free(mysmap[iam]);
	free(myxg1[iam]);
	free(myxg2[iam]);
	free(mye1[iam]);
	free(mye2[iam]);
	if (prefs.doweight)
	    free(mysigmae[iam]);
	if (prefs.dotomo) {
	    free(myzbin[iam]);
	}
    }

    free(mymap);
    free(mysmap);
    free(myxg1);
    free(myxg2);
    free(mye1);
    free(mye2);
    if (prefs.doweight)
	free(mysigmae);
    if (prefs.dotomo) {
	free(myzbin);
    }

    msg_message(MSG_INFO, 2, "Done");
    if (prefs.smap) {
	for(in=0; in<(map.naxes[0]*map.naxes[1]); in++) {
	    if (smap.arr[in])
		smap.arr[in] =  M_SQRT2l * map.arr[in] / sqrt(smap.arr[in]);
	}
    } else if (prefs.dotomo) {
	for(i=0; i<(map.naxes[0]*map.naxes[1]); i++) {
	    in = i + lp*map.naxes[0]*map.naxes[1];
	    map.arr[in] = SQR(map.arr[in])/smap.arr[in];
	}
    }
}




/* Compute redshift weight for tomography. */
double z_weight(int bin) {

    if (bin >= ZSTEPS) {
	/* This galaxy is beyond the maximum redshift */
	return 0;
    }
    return zweightarr[bin];
}


void fill_zweight_array(double zd) {
    int i;
    double z;

    msg_message(MSG_INFO, 0, "Filling new weight array z=%.3lf ", zd);
    for(i=0; i<ZSTEPS; i++) {
	z = ZMAX/(ZSTEPS-1.0) * i;
	if (z>zd)
	    zweightarr[i] = dang(0, zd, cospar) * dang(zd, z, cospar) / ds[i];
	else
	    zweightarr[i] = 0;
    }    
    msg_message(MSG_INFO, 0, "Done");
}


void fill_ds_array(void) {
    int i;
    double z;
    
    ds[0] = 0;
    zindex[0] = 0;
    msg_message(MSG_INFO, 0, "Initial filling of ds array");
    for(i=1; i<ZSTEPS; i++) {
	z = ZMAX/(ZSTEPS-1.0) * i;
	zindex[i] = z;
	ds[i] = dang(0, z, cospar);
    }
    msg_message(MSG_INFO, 0, "done");
}



void randomize(double *e1, double *e2, unsigned long ng) {
    double p, c, s;
    double tmp;
    unsigned long i;
    
    srand48(21101975);

    for (i=0; i<ng; i++) {
	p = drand48() * 2. * M_PI;
	c = cos(p);
	s = sin(p);
	tmp = c*e1[i] + s*e2[i];
	e2[i] = -s*e1[i]+c*e2[i];
	e1[i] = tmp;
    }
    return;
}


/* Rotate ellipticities by 45 deg. */
void cross(double *e1, double *e2, unsigned long ng) {
    double tmp;
    unsigned long i;
    
    for (i=0; i<ng; i++) {
	tmp = e2[i];
	e2[i] = e1[i];
	e1[i] = -tmp;
    }
}


/* Return angular separation in radians for object projected
   separation proj_dist at redshift z */
double get_z_ang(double proj_dist, double z) {
    double ang_dist, angle;
    
    ang_dist = dang(0, z, cospar);
    angle = proj_dist/(1000*ang_dist);
    return angle;
}


int writefitsmap(mapstruct map, char *outname, char *history) {
    int status = 0;
    fitsfile *fptr;

    msg_message(MSG_NOTICE, 0, "Writing fits file %s", outname);
    remove(outname);
    fits_create_file(&fptr, outname, &status);
    if (prefs.dotomo)
	fits_create_img(fptr, FLOAT_IMG, 3, map.naxes, &status);
    else
	fits_create_img(fptr, FLOAT_IMG, 2, map.naxes, &status);
    if (map.wcs) {
	fits_write_key_str(fptr, "CTYPE1", map.wcs->ctype1, 
			   "Projection type for RA", &status);
	fits_write_key_str(fptr, "CTYPE2", map.wcs->ctype2, 
			   "Projection type for Dec", &status);

	fits_write_key_flt(fptr, "CRPIX1", map.wcs->crpix1, -5, 
			   "Reference pixel x", &status);
	fits_write_key_flt(fptr, "CRPIX2", map.wcs->crpix2, -5,
			   "Reference pixel y", &status);

	fits_write_key_flt(fptr, "CRVAL1", map.wcs->crval1, -9,
			   "RA at reference pixel", &status);
	fits_write_key_flt(fptr, "CRVAL2", map.wcs->crval2, -9,
			   "Dec at reference pixel", &status);


	fits_write_key_flt(fptr, "CDELT1", map.wcs->cdelt1, 9, 
			   "Pixel scale in x", &status);
	fits_write_key_flt(fptr, "CDELT2", map.wcs->cdelt2, 9, 
			   "Pixel scale in y", &status);

	fits_write_key_flt(fptr, "EQUINOX", map.wcs->equinox, -5,
			   "Equinox", &status);
	fits_write_key_str(fptr, "FILTER", prefs.filter,
			   "filter function", &status);
	fits_write_key_flt(fptr, "FRADIUS", prefs.fradius_orig, -4,
			   "filter radius in kpc", &status);
	fits_write_key_flt(fptr, "PIXSCALE", prefs.spix/1000., -10,
			   "pixel scale in Mpc at redshift", &status);
	fits_write_key_flt(fptr, "REDSHIFT", prefs.zcl, -5,
			   "redshift of fiducial model", &status);
	if (prefs.dotomo) {
	    fits_write_key_str(fptr, "CTYPE3", map.wcs->ctype3, 
			       "Projection type for redshift", &status);
	    fits_write_key_flt(fptr, "CRPIX3", map.wcs->crpix3, -5,
			       "Reference pixel z", &status);
	    fits_write_key_flt(fptr, "CRVAL3", map.wcs->crval3, -9,
			       "Redshift at reference pixel", &status);
	    fits_write_key_flt(fptr, "CDELT3", map.wcs->cdelt3, 9, 
			       "Pixel scale in z", &status);
	}
    }
    fits_write_history(fptr, history, &status);

    if (status) 
	fits_error("fits_read_col: Error reading column from table.", 
		   status);

    fits_write_img(fptr, TDOUBLE, 1, map.naxes[0]*map.naxes[1]*map.naxes[2], 
		   map.arr, &status);
    fits_close_file(fptr, &status);
    msg_message(MSG_INFO, 0, "Done.");

    return status;
}


void readfitscol(fitsfile *fptr, char *colnam, unsigned long nrows, 
		 double *arr) {
    int status = 0;
    int colnum;
    char fits_err[FITS_ERR_LEN];
    
    fits_get_colnum(fptr, CASESEN, colnam, &colnum, &status);
    if (status) {
	snprintf(fits_err, 128, "fits_get_colnum: No such column: %s", 
		 colnam);
	fits_error(fits_err, status);
    }
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, arr, 0, &status);
    if (status) 
	fits_error("fits_read_col: Error reading column from table.", 
		   status);
}


