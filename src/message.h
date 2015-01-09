#ifndef _MESSAGE_H_
#define _MESSAGE_H_ 1

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "../config.h"

#include "prefs.h"

enum severity {MSG_FATAL, MSG_ERROR, MSG_WARNING, MSG_NOTICE, MSG_INFO,
	       MSG_DEBUG};

typedef struct {
#ifdef HAVE_GETTIMEOFDAY
    double init_time;
#else
    time_t init_time;
#endif
    char *logfile;
    FILE *f;
    int verbosity;
} msg_state;


void msg_init();
void msg_finish();
msg_state *msg_set_get_state(int set, int verbosity, char *logfile);
void msg_message(int severity, int overwrite, char *format, ...);
void msg_format_time(char *buf, float t);


#endif
