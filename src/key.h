/* $Id: key.h,v 1.1 2006/09/06 07:51:17 cvs Exp $
 * $Date: 2006/09/06 07:51:17 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

/*--------------------------------- constants -------------------------------*/

#define         FIND_STRICT     0
#define         FIND_NOSTRICT   1


/*--------------------------- structure definitions -------------------------*/
/* Preference keyword */

typedef struct
{
    char          name[16];
    enum  {P_BOOL, P_INT, P_FLOAT, P_STRING} type;
    void          *ptr;                 /* Pointer to the keyword value */
    int           imin, imax;           /* Range for int's */
    float         dmin, dmax;           /* Range for floats */
    char          keylist[16][16];      /* List of keywords */
    int           nlistmin;             /* Minimum number of list members */
    int           nlistmax;             /* Maximum number of list members */
    int           *nlistptr;            /* Ptr to store the nb of read params*/
    int           flag;
} pkeystruct;
