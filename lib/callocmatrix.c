/* $Id: callocmatrix.c,v 1.1 2006/09/11 14:18:18 cvs Exp $
 * $Date: 2006/09/11 14:18:18 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

#include <stdlib.h>

#include "misc.h"

extern void exit_failure(const char *err);

void *critcalloc(size_t nmemb, size_t size, const char *err) {
    void *ptr;

    ptr = calloc(nmemb, size);
    if (!ptr)
	exit_failure(err);
    return ptr;
}


void *critrealloc(void *ptr, size_t size, const char *err) {
    ptr = realloc(ptr, size);
    if (!ptr)
	exit_failure(err);
    return ptr;
}
	
/* calloc a matrix with m lines and n columns */

double **callocmatrix(unsigned long m, unsigned long n) {
  
  double **matrix;
  unsigned long i;
  
  matrix = (double **)calloc(m, sizeof(double *));
  if (!matrix) 
    return NULL;

  matrix[0] = (double *)calloc(m*n, sizeof(double));
  if (!matrix[0])
    return NULL;

  for(i=1; i<m; i+=1)
    matrix[i] = matrix[i-1] + n;

  return matrix;
}


/* calloc a matrix with m lines and n columns */

unsigned short int **callocusimatrix(unsigned long m, unsigned long n) {
  
  unsigned short int **matrix;
  unsigned long i;
  
  matrix = (unsigned short int **)calloc(m, sizeof(unsigned short int *));
  if (!matrix) 
    return NULL;

  matrix[0] = (unsigned short int *)calloc(m*n, sizeof(unsigned short int));
  if (!matrix[0])
    return NULL;

  for(i=1; i<m; i+=1)
    matrix[i] = matrix[i-1] + n;

  return matrix;
}

void freematrix(double **matrix) {
  
  free((void *)matrix[0]);
}

 
void freeusimatrix(unsigned int **matrix) {
    free((void *)matrix[0]);
}
