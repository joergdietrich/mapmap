/* $Id: preflist.h,v 1.4 2007/08/13 16:17:32 jdietric Exp $
 * $Date: 2007/08/13 16:17:32 $
 * $Author: jdietric $
 * $Revision: 1.4 $
 */

#include "key.h"

#include "prefs.h"

pkeystruct key[] = {
    {"HUBBLE", P_FLOAT, &prefs.h, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"OMEGA_MATTER", P_FLOAT, &prefs.Omega_M, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"OMEGA_LAMBDA", P_FLOAT, &prefs.Omega_L, 0, 0, -BIG, BIG, {""}, 0, 1, NULL, 0},

    {"PIXEL_SIZE", P_FLOAT, &prefs.spix, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"REDSHIFT", P_FLOAT, &prefs.zcl, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"FILTER", P_STRING, &prefs.filter, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"FILTER_RADIUS", P_FLOAT, &prefs.fradius, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"FILTER_FILE", P_STRING, &prefs.ffile, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"OVERSHOOT", P_FLOAT, &prefs.ovrsht, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"SCHIRMER_XC", P_FLOAT, &prefs.xc, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"TABLE", P_STRING, &prefs.tblname, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"RA", P_STRING, &prefs.xname, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"DEC", P_STRING, &prefs.yname, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"E1", P_STRING, &prefs.e1name, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"E2", P_STRING, &prefs.e2name, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"WEIGHT", P_BOOL, &prefs.doweight, 0, 0, 0.0, 0, {""}, 0, 0, NULL, 0},
    {"WEIGHT_NAME", P_STRING, &prefs.weightname, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"TOMOGRAPHY", P_BOOL, &prefs.dotomo, 0, 0, 0.0, 0, {""}, 0, 0, NULL, 0},
    {"Z_NAME", P_STRING, &prefs.zname, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"OUTPUT", P_STRING, &prefs.outname, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"MAP_SUFFIX", P_STRING, &prefs.mapsuf, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"SMAP", P_BOOL, &prefs.smap, 0, 0, 0.0, 0, {""}, 0, 0, NULL, 0},
    {"SMAP_SUFFIX", P_STRING, &prefs.smapsuf, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"CROSS", P_BOOL, &prefs.cross, 0, 0, 0.0, 0, {""}, 0, 0, NULL, 0},
    {"CROSS_SUFFIX", P_STRING, &prefs.xsuf, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"INVERT", P_BOOL, &prefs.invert, 0, 0, 0.0, 0, {""}, 0, 0, NULL, 0},
    {"INVERT_SUFFIX", P_STRING, &prefs.invsuf, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"MOCK", P_BOOL, &prefs.mock, 0, 0, 0.0, 0, {""}, 0, 0, NULL, 0},
    {"MOCK_SUFFIX", P_STRING, &prefs.mocksuf, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"ZMIN", P_FLOAT, &prefs.zmin, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"ZMAX", P_FLOAT, &prefs.zmax, 0, 0, 0.0, BIG, {""}, 0, 1, NULL, 0},
    {"N_PLANES", P_INT, &prefs.nplanes, 0, 512, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"VERBOSE", P_INT, &prefs.verbose, 0, 5, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"NTHREADS", P_INT, &prefs.nthreads, 0, 512, 0.0, 0.0, {""}, 0, 1, NULL, 0},
    {"", 0, NULL, 0, 0, 0.0, 0.0, {""}, 0, 1, NULL, 0}
};


char                    keylist[sizeof(key)/sizeof(pkeystruct)][16];
const char              notokstr[] = {" \t=,;\n\r\""};


char *default_prefs[] = {
    "# Default configuration file for mapmap",
    "# JPD 18.1.2006",
    "#",
    "#------------------------------ Cosmology --------------------------------------",
    "HUBBLE        0.7    # scaled Hubble constant H_0 = 100 h km/s/Mpc",
    "OMEGA_MATTER  0.3    # Matter density",
    "OMEGA_LAMBDA  0.7    # Cosmological constant",
    " ",
    "#---------------------------- Cluster Model ------------------------------------",
    "PIXEL_SIZE     50    # Pixel size in kpc",
    "REDSHIFT       0.3   # Redshift of fiducial cluster model",
    "FILTER         SCHIRMER # The type of filter function",
    "FILTER_RADIUS  1500  # Filter radius in kpc",
    "FILTER_FILE    NONE  # File to read if FILTER == FILE",
    "OVERSHOOT      1.0   # The filter function is applied out to OVERSHOOT*FILTER_RADIUS",
    "SCHIRMER_XC    0.15  # xc value in Schirmer filter",
    " ",
    "#---------------------------- Input Catalog ------------------------------------",
    "TABLE         OBJECTS     # Table name containing the galaxies",
    "RA            ALPHA_J2000 # Column name for right ascencion",
    "DEC           DELTA_J2000 # Column name for declination",
    "E1            e1iso_RG    # Column name for e1 component",
    "E2            e2iso_RG    # Column name for e2 component",
    "WEIGHT        N           # Weigh ellipticity values",
    "WEIGHT_NAME   Sigmae      # Column name for weight",
    "TOMOGRAPHY    N           # Perform weak lensing tomography",
    "Z_NAME        z           # Column name for redshift",
    " ",
    "#---------------------------- Program Setup ------------------------------------",
    "OUTPUT        mapmap      # Default output basename",
    "MAP_SUFFIX    _map        # Suffix for Map maps",
    "SMAP          Y           # Make significance map",
    "SMAP_SUFFIX   _smap       # Suffix for significance maps",
    "CROSS         N           # Rotate all ellipticities by 45 deg",
    "CROSS_SUFFIX  _cross      # Suffix for cross-component images",
    "INVERT        N           # Output inverted (negative) images",
    "INVERT_SUFFIX _invert     # Suffix for inverted images",
    "MOCK          N           # Create mock field (random rotation)",    
    "MOCK_SUFFIX   _mock       # Suffix for mock images",
    "ZMIN          0.1         # Minimum redshift for tomography",
    "ZMAX          1.0         # Maximum redshift for tomography",
    "N_PLANES      15          # Number of redshift planes",
    "VERBOSE       3           # Verbosity of output",
    "NTHREADS      0           # Number of OpenMP threads (0: use system default)",
    ""
};

