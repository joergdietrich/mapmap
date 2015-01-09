#include "filter.h"

/* Get critical density in M_sun /Mpc^3 at redshift z */
double get_rhoc(double z, /* redshift */ 
		cos_param cospar) {
    double rhoc;
    double az, hz;

    az = 1./(1+z);
    hz = sqrt(SQR(cospar.h*100) * (cospar.omega_m/pow(az, 3) 
				   + (1-(cospar.omega_m+cospar.omega_l)/
				      pow(az, 2) + cospar.omega_l))); 
    rhoc = 3 * SQR(hz)/(8*M_PI*G);
    rhoc *= 1e6 * MPC / MSUN; 
    return rhoc;
}

double nfw_shear(nfw_struct nfw,
		 double z_nfw,      /* redshift of NFW halo*/		 
		 double zb,         /* redshift of background source */
		 double r,          /* position in Mpc */
		 cos_param cospar) {

    double val, fac, rs, r200, x;
    double rhoc, deltac, Sigmac;

    Sigmac = get_Sigmac(z_nfw, zb, cospar);
    deltac = 200./3. * pow(nfw.c, 3)/(log(1+nfw.c)-nfw.c/(1+nfw.c));
    rhoc = get_rhoc(z_nfw, cospar);
    r200 = pow(3*nfw.m200/(800*M_PI*rhoc), 1./3.);
    rs = r200/nfw.c;
    x = r/rs; 

    fac = rs*deltac*rhoc/Sigmac;
    val = nfw_g(x);
    return val*fac;
}


/* Return critical surface mass density for a lens at redshift zd and
   background source at redshift zs.
 */
double get_Sigmac(double zd, double zs, cos_param cospar) {
    double dd, ds, dds;
    double Sigmac;

    dd = dang(0, zd, cospar);
    ds = dang(0, zs, cospar);
    dds = dang(zd, zs, cospar);

    Sigmac = SQR(C)/(4.*M_PI*G)*ds/(dd*dds);
    Sigmac *= MPC/MSUN;
}
    


/* the g function of Wright and Brainerd (2003) */
double nfw_g(double x) {
    double val = 0;   /* Shut compiler up */
    double xsq;

    xsq = SQR(x);
    if (fabs(x-1) < 1.0e-9) 
	val = 10./3. + 4*log(0.5);
    else if (x < 1 ) {
	val = 8.*atanh(sqrt((1-x)/(1+x)))/(xsq*sqrt(1-xsq));
	val += 4./xsq * log(x/2.);
	val -= 2./(xsq-1);
	val += 4.*atanh(sqrt((1-x)/(1+x)))/((xsq-1)*sqrt(1-xsq));
    } else if (x > 1) {
	val = 8.*atan(sqrt((x-1)/(1+x)))/(xsq*sqrt(xsq-1));
	val += 4./xsq * log(x/2.);
	val -= 2./(xsq-1);
	val += 4.*atan(sqrt((x-1)/(1+x)))/pow(xsq-1, 1.5);
    }
    return val;
}


double nfw_kappa(nfw_struct nfw,
		 double z_nfw,      /* redshift of NFW halo*/		 
		 double zb,         /* redshift of background source */
		 double r,          /* radius in Mpc */
		 cos_param cospar) {

    double val, fac, rs, r200, x;
    double rhoc, deltac, Sigmac;

    Sigmac = get_Sigmac(z_nfw, zb, cospar);
    deltac = 200./3. * pow(nfw.c, 3)/(log(1+nfw.c)-nfw.c/(1+nfw.c));
    rhoc = get_rhoc(z_nfw, cospar);
    r200 = pow(3*nfw.m200/(800*M_PI*rhoc), 1./3.); /*r200 in Mpc */
    rs = r200/nfw.c;
    x = r/rs;

    fac = 2.*rs*deltac*rhoc/(SQR(x)-1.);
    if (fabs(x-1) < 1.0e-9) {
	val = 2.*rs*deltac*rhoc/3.;
    } else if (x<1) {
	val = 1. - 2./sqrt(1-SQR(x)) * atanh(sqrt((1-x)/(1+x)));
	val *= fac;
    } else {
	val = 1. - 2./sqrt(SQR(x)-1) * atan(sqrt((x-1)/(1+x)));
	val *= fac;
    }

    return val/Sigmac;
}


double weight(double x, void *param) {
    weight_param *p = (weight_param *)param;

    if (!strcmp(prefs.filter, "SCHNEIDER")) 
	return schneider_weight(x);
    else if (!strcmp(prefs.filter, "SCHIRMER")) {
	double xc;
	xc = p->xc;
	return schirmer_weight(x, xc);
    }
    else if (!strcmp(prefs.filter, "FILE")) {
	double r, c;

	r = p->r;
	c = p->c;
	x *= r*c*60.0;
	return file_weight(x);
    }
    else
	exit_failure("Unknown filter type\n");
    
    /* shut compiler warning up */
    return 0;
}


double schneider_weight(double x) {
    double q;

    q = SQR(x) - pow(x, 4);
    return q;
}


double schirmer_weight(double x, double xc) {
    double q;
    float a, b, c, d;

    a = 6;
    b = 150;
    c = 47;
    d = 50;

    q = tanh(x/xc)/(x/xc) * 1./(1+exp(a-b*x) + exp(-c+d*x));

    return q;
}

double file_weight(double x) {
    static double *scale;
    static double *value;
    static unsigned int len;
    
    double val;
    int status = 0;
    char err[MAXCHAR];
    
    /* First time, read the file */
    if (!scale && !value) {
	FILE *f;
	char *l, *p, *q;
	float s, v;

	f = fopen(prefs.ffile, "r");
	if (!f) {
	    snprintf(err, MAXCHAR, "Could not open filter file %s", 
		     prefs.ffile);
	    exit_failure(err);
	}
	len = 1;
	scale = malloc(sizeof(double));
	value = malloc(sizeof(double));
	if (!scale || !value)
	    exit_failure("Could not allocate scale/value");

	while((l=getaline(f))) {
	    p = l;
	    /* skip comments */
	    q = strchr(p, '#');
	    if (q) {
		*q = '\0';
	    }
	    if (*p == '\0')
		continue;
	    SKIPLWS(p);
	    q = p;
	    while(q && !isspace(*q))
		*q++;
	    *q++ = '\0';
	    SKIPLWS(q);
	    s = atof(p);
	    v = atof(q);
	    scale = critrealloc(scale, (len+1)*sizeof(double), 
				"realloc(scale) failed");
	    value = critrealloc(value, (len+1)*sizeof(double),
				"realloc(value) failed"); 
	    scale[len] = s;
	    value[len] = v;
	    len++;
	}
	fclose(f);
	len--;
    }
    linint(scale, value, len, x, &val, &status);
    if (status)
	return 0;
    return val;
}


void linint(double *x, double *y, unsigned long len, double xm, double *ym, 
	    int *status) {
    unsigned long i;

    i = 1;
    while (x[i]<xm && i<=len) 
	i++;
    if (i == 1 || i == len) {
	fprintf(stderr, "Target %lf out of range. Cannot interpolate.\n", xm);
	*status = 1;
	return;
    }
    *ym = y[i-1] + (xm-x[i-1])* (y[i]-y[i-1])/(x[i]-x[i-1]);
}

