#include "message.h"

char *msg_severity[] = {"FATAL  ", 
			"ERROR  ", 
			"WARNING", 
			"NOTICE ", 
			"INFO   ", 
			"DEBUG  "};


/*
  Initialize the msg interface. Set init_time, verbosity log level,
  open logfile if file logging is requested.
*/
void msg_init(int verbosity, char *logfile) {
    msg_set_get_state(1, verbosity, logfile);
    return;
}

void msg_finish() {
    msg_state *state;

    state = msg_set_get_state(0, 0, NULL);
    if (!state)
	return;
    if (state->f) {
	msg_message(MSG_DEBUG, 0, "Closing log file");
	fclose(state->f);
    }
    if (state->logfile) 
	free(state->logfile);

    free(state);
    return;
}

msg_state *msg_set_get_state(int set, int verbosity, char *logfile) {

    static msg_state *state;
#ifdef HAVE_GETTIMEOFDAY
    struct timeval tv;
#endif

    if (!state || set) {
	if (state)
	    msg_finish();
	state = malloc(sizeof(msg_state));
#ifdef HAVE_GETTIMEOFDAY
	gettimeofday(&tv, NULL);
	state->init_time = tv.tv_sec + tv.tv_usec/1000000.0;
#else	
	state->init_time = time(NULL);
#endif
	state->verbosity = verbosity;
	if (!set) {
	    state->verbosity = MSG_WARNING;
	    msg_message(MSG_WARNING, 0, "Trying to message interface without initializing it. Continuing with default settings.");
	}
	if (logfile) {
	    state->logfile = strdup(logfile);
	    state->f = fopen(state->logfile, "w");
	    if (!state->f) {
		msg_message(MSG_ERROR, 0, "Could not open logfile %s", 
			    state->logfile);
		state->logfile = NULL;
	    }
	} else {
	    state->logfile = NULL;
	    state->f = NULL;
	}
    }
    return state;
}
	
void msg_message(int severity, int overwrite, char *format, ...) {
    
    msg_state *state;
    double now;
    char buf[16], msg[1024];
    va_list ap;
    FILE *stream;
#ifdef HAVE_GETTIMEOFDAY
    struct timeval tv;
#endif

    state = msg_set_get_state(0, 0, NULL);
    if (severity <= state->verbosity) {
	va_start(ap, format);
#ifdef HAVE_GETTIMEOFDAY
	gettimeofday(&tv, NULL);
	now = tv.tv_sec + tv.tv_usec/1000000.0;
#else
	now = (float)time(NULL);
#endif
	msg_format_time(buf, now - state->init_time);
	vsnprintf(msg, 1024, format, ap);
	if (severity <= MSG_ERROR)
	    stream = stderr;
	else
	    stream = stdout;
	if (overwrite == 1)
	    fputc('\r', stream);
	else if (overwrite == 2)
	    fputc('\n', stream);
	fprintf(stream, "[%s%09s] %s", msg_severity[severity], buf, msg);
	if (overwrite!=1)
	    fputc('\n', stream);
	else 
	    fflush(stream);
	if (state->f)
	    fprintf(state->f, "[%s%09s] %s\n", msg_severity[severity], buf, 
		    msg);
    }
}
	

/*
  print the formatted time t into the preallocated buffer buf (min
   size 16) 

   char *buf - preallocated array (min size 16) 
   time_t t - * time in seconds to be formatted
*/
void msg_format_time(char *buf, float t) {
    int day, hour, min, sec, cent;

    t = floor(t * 100.0 + 0.5);
    cent = (int)(t - floor(t / 100.0) * 100.0);
    t = floor(t / 100.0);
    sec = (int)(t - floor(t / 60.0) * 60.0);
    t = floor(t / 60.0);
    min = (int)(t - floor(t / 60.0) * 60.0);
    t = floor(t / 60.0);
    hour = (int)(t - floor(t / 24.0) * 24.0);
    t = floor(t) / 24.0;
    day = (int)(floor(t));
    if (day >= 10) 
	snprintf(buf, 16, "%-2dd %02dh", day, hour);
    else if (day > 0) 
	snprintf(buf, 16, "%1dd%02d:%02d", day, hour, min);
    else if (hour >= 10) 
	snprintf(buf, 16, "%-2dh %02dm", hour, min);
    else if (hour > 0) 
	snprintf(buf, 16, "%02d:%02d:%02d", hour, min, sec);
    else if (min >= 10) 
	snprintf(buf, 16, "%-2d:%02d.%02d", min, sec, cent/10);
    else 
	snprintf(buf, 16, "%02d:%02d.%02d", min, sec, cent);
}
