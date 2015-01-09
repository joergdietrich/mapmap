/* $Id: define.h,v 1.5 2007/08/14 08:13:34 jdietric Exp $
 * $Date: 2007/08/14 08:13:34 $
 * $Author: jdietric $
 * $Revision: 1.5 $
 */

 /*
 				define.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	mapmap
*
*	Author:		J. P. Dietrich, AIfA, University of Bonn
*
*	Contents:	global definitions.
*
*	Last modify:	2006-01-20
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------ what, who, when and where ------------------------*/

#define		BANNER		"mapmap"
#define		COPYRIGHT	"J.P. Dietrich <jdietric@eso.org>"
#define		INSTITUTE	"ESO Garching"


/*----------------------------- Internal constants --------------------------*/

#define	MAXCHAR			256		/* max. number of characters */
#define	BIG			1e+30		/* a huge number */


/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)


/*------------------------------- Other Macros -----------------------------*/

static double sqrarg __attribute__ ((unused));
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double marg1 __attribute__ ((unused)), marg2 __attribute__ ((unused));
#define MIN(a,b) (marg1=(a), marg2=(b), (marg1<marg2) ? marg1 : marg2)
#define MAX(a,b) (marg1=(a), marg2=(b), (marg1>marg2) ? marg1 : marg2)

#define FITS_ERR_LEN (128)

#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif

#ifndef SQRT2l
# define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
#endif


#define SYNTAX \
"mapmap catalog [-c <configuration_file>] [-<keyword> <value>]\n" \
" or, to dump a default configuration file:\n" \
" mapmap -d \n"

