/* $Id: misc.c,v 1.1 2006/09/06 07:54:01 cvs Exp $
 * $Date: 2006/09/06 07:54:01 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

#include <stdio.h>

#include <fitsio.h>

#include "define.h"
#include "message.h"


void exit_failure(const char *err) {
    msg_message(MSG_FATAL, 0, (char *)err);
    perror(err);
    exit(EXIT_FAILURE);
}


void fits_error(const char *err, int status)  {
    msg_message(MSG_FATAL, 0, (char *)err);
    fits_report_error(stderr, status);
    exit(EXIT_FAILURE);
}

void min_max_center(double *arr, long len, double *min, double *max, 
                    double *center) {
    long i;
    
    if (min)
        *min = BIG;
    if (max)
        *max = -BIG;
    if (center)
        *center = 0;

    for(i=0; i<len; i++) {
        if (min && arr[i] < *min) 
            *min = arr[i];
        if (max && arr[i] > *max)
            *max = arr[i];
        if (center)
            *center += arr[i];
    }
    if (center)
        *center /= len;
}
