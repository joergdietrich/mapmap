AM_CFLAGS = -I$(top_srcdir)/lib
AM_CFLAGS += -W -Wall
if PLPA_OPT
BUILD_OPT = $(top_srcdir)/plpa-1.1/src/libplpa/libplpa_included.la
endif
LDADD = -lcosmology -lgetaline -lcallocmatrix $(BUILD_OPT)
AM_LDFLAGS = -g
bin_PROGRAMS = mapmap
mapmap_LDFLAGS = -L../lib
mapmap_SOURCES = mapmap.c misc.c prefs.c filter.c message.c define.h key.h filter.h \
	mapmap.h misc.h message.h \
	preflist.h prefs.h types.h

