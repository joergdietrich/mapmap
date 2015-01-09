/* $Id: types.h,v 1.2 2007/08/14 16:02:25 jdietric Exp $
 * $Date: 2007/08/14 16:02:25 $
 * $Author: jdietric $
 * $Revision: 1.2 $
 */

#ifndef _TYPES_H_
#define _TYPES_H_ 1

struct wcs {
    float crval1;
    float crval2;
    float crval3;
    float crpix1;
    float crpix2;
    float crpix3;
    float cdelt1;
    float cdelt2;
    float cdelt3;
    char *ctype1;
    char *ctype2;
    char *ctype3;
    float equinox;
};
typedef struct wcs wcsstruct;


struct img {
    double *arr;        /* array containing data */
    long naxes[3];      /* The length of the axes */
    wcsstruct *wcs;     /* The wcs information */
};
typedef struct img mapstruct;



#endif
