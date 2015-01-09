/* $Id: callocmatrix.h,v 1.1 2006/09/11 14:18:40 cvs Exp $
 * $Date: 2006/09/11 14:18:40 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

void *critcalloc(size_t, size_t, const char *);
void *critrealloc(void *, size_t, const char *);
double **callocmatrix(unsigned long, unsigned long);
unsigned short int **callocusimatrix(unsigned long, unsigned long);
void freematrix(double **);
void freeusimatrix(unsigned short int **);
