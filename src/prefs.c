/* $Id: prefs.c,v 1.1 2006/09/06 08:00:04 cvs Exp $
 * $Date: 2006/09/06 08:00:04 $
 * $Author: cvs $
 * $Revision: 1.1 $
 */

/*
 				prefs.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	match_filt
*
*	Author:		E.BERTIN (IAP), Joerg Dietrich
*
*	Contents:	Functions to handle the configuration file.
*
*	Last modify:	
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	<stdio.h>
#include	<stdlib.h>
#include        <ctype.h>
#include        <string.h>

#include        "define.h"
#include        "misc.h"
#include 	"prefs.h"
#include	"preflist.h"
#include        "message.h"

/********************************* dumpprefs ********************************/
/*
Print the default preference parameters.
*/
void	dumpprefs(void)
{
    char	**dp;
    
    dp = default_prefs;
    while (**dp)
	if (**dp != '*')
	    printf("%s\n",*(dp++));
	else
	    dp++;
    return;
}

/********************************* readprefs ********************************/
/*
Read a configuration file in ``standard'' format (see the SExtractor
documentation)
*/
void    readprefs(char *filename, char **argkey, char **argval, int narg)
{

    FILE *infile;
    char *cp, str[PATH_MAX], *keyword, *value, **dp;
    int i, ival, nkey, warn, argi, flagc, flagd, flage;
    float dval;

    if ((infile = fopen(filename,"r")) == NULL)
	{
	    flage = 1;
	    msg_message(MSG_WARNING, 0, "%s not found, using internal defaults",
			filename);
	}
    else
      flage = 0;
    
    /*Build the keyword-list from pkeystruct-array */
    
    for (i=0; key[i].name[0]; i++)
	strcpy(keylist[i], key[i].name);
    keylist[i][0] = '\0';
    
    /*Scan the configuration file*/

    argi=0;
    flagc = 0;
    flagd = 1;
    dp = default_prefs;
    for (warn=0;;) {
	if (flagd) {
	    if (**dp) {
		if (**dp=='*')
		    strcpy(str, *(dp++)+1);
		else
		    strcpy(str, *(dp++));
	    }
	    else
		flagd = 0;
	}
	if (!flagc && !flagd)
	    if (flage || !fgets(str, MAXCHAR, infile))
		flagc=1;
	
	if (flagc) {
	    if (argi<narg) {
		sprintf(str, "%s %s", argkey[argi], argval[argi]);
		argi++;
	    }
	    else
		break;
	}
	
	keyword = strtok(str, notokstr);
	if (keyword && keyword[0]!=0 && keyword[0]!=(char)'#') {
	    if (warn>=10) {
		msg_message(MSG_FATAL, 0, "No valid keyword found in %s", 
			    filename);
		exit(EXIT_FAILURE);
	    }
	    nkey = findkeys(keyword, keylist);
	    if (nkey!=RETURN_ERROR) {
		value = strtok((char *)NULL, notokstr);
		switch(key[nkey].type)
		    {
		    case P_BOOL:
			if (!value || value[0]==(char)'#') {
			    msg_message(MSG_FATAL,0,"keyword %s has no value!",
					keyword);
			    exit(EXIT_FAILURE);
			}
			if ((cp = strchr("yYnN", (int)value[0])) != NULL)
			    *(int *)(key[nkey].ptr) = (tolower((int)*cp)=='y')?1:0;
			else {
			    msg_message(MSG_FATAL, 0, keyword, 
					"%s value must be Y or N", keyword);
			    exit(EXIT_FAILURE);
			}
			break;
			
		    case P_INT:
			if (!value || value[0]==(char)'#') {
			    msg_message(MSG_FATAL,0,"keyword %s has no value!",
					keyword);
			    exit(EXIT_FAILURE);
			}
			ival = atoi(value);
			if (ival>=key[nkey].imin && ival<=key[nkey].imax)
			    *(int *)(key[nkey].ptr) = ival;
			else {
			    msg_message(MSG_FATAL,0,"keyword %s out of range",
					keyword);
			    exit(EXIT_FAILURE);
			}
			break;
			
		    case P_FLOAT:
			if (!value || value[0]==(char)'#') {
			    msg_message(MSG_FATAL,0,"keyword %s has no value!",
					keyword);
			    exit(EXIT_FAILURE);
			}
			dval = atof(value);
			if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
			    *(double *)(key[nkey].ptr) = dval;
			else {
			    msg_message(MSG_FATAL,0,"keyword %s out of range",
					keyword);
			    msg_message(MSG_FATAL, 0, "nmin: %lf\tmax: %lf", 
					key[nkey].dmin, key[nkey].dmax);
			    exit(EXIT_FAILURE);
			}
			break;
			
		    case P_STRING:
			if (!value || value[0]==(char)'#') {
			    msg_message(MSG_FATAL, 0, "%s string is empty",
					keyword);
			    exit(EXIT_FAILURE);
			}
			strcpy((char *)key[nkey].ptr, value);
			break;
		    default:
			msg_message(MSG_FATAL, 0, 
				    "Internal error: Type Unknown in readprefs()");
			exit(EXIT_FAILURE);
			break;
		    }
		key[nkey].flag = 1;
	    } else {
		msg_message(MSG_WARNING, 0, "keyword %s unknown", keyword);
		warn++;
	    }
	}
    }
    
    for (i=0; key[i].name[0]; i++)
	if (!key[i].flag) {
	    msg_message(MSG_FATAL, 0, "configuration keyword %s missing", 
			key[i].name);
	}
    if (!flage)
	fclose(infile);
    
    return;
}

/********************************* findkeys **********************************/
/*
find an item within a list of keywords.
*/
int findkeys(char *str, char keyw[][16])
{
  int i;

  for (i=0; keyw[i][0]; i++)
      if (!strcasecmp(str, keyw[i]))
	  return i;
  
  return RETURN_ERROR;
}

			    
