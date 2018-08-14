/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2012: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is dual-licensed as follows:                                 *
*                                                                           *
* (1) You may use TMesh as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use TMesh as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#ifndef _BASICS_H
#define _BASICS_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "coordinates.h"

namespace T_MESH
{

#ifdef EXTENSIBLE_TMESH
#define TMESH_VIRTUAL virtual
#else
#define TMESH_VIRTUAL
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TMESH_VERSION	"2.0"
#define TMESH_YEAR		2012

class TMesh
{
 public:

 static void (*display_message)(const char *, int);
 static char *app_name;
 static char *app_version;
 static char *app_year;
 static char *app_authors;
 static char *app_url;
 static char *app_maillist;
 static coord maximum_coord_value; // No other coord number can be larger than this

 static const char *filename; // This might be null. If not, it represents the file we are currently working with.

 static bool quiet;

 static void init(void (*)(const char *, int) = NULL);

 static void info(const char *, ...);
 static void warning(const char *, ...);
 static void error(const char *, ...);
 static void begin_progress();
 static void report_progress(const char *, ...);
 static void end_progress();

 //! When called with a nonzero argument 'ts', launches a chronometer with a timeout of 'ts' seconds.
 //! Later calls without arguments check the chronometer, and if it is over 'ts' the program exits.
 static void exitOnTimeout(clock_t ts = 0);

 //! Appends a line to the file "tmesh.log"
 static void addMessageToLogFile(const char *msg);

 //! Formats a message headed with date/time/filename, appends it to "tmesh.log", and exits with error
 static void logToFileAndExit(const char *msg);

 //! When called without arguments prints the elapsed time from the latest reset.
 static void printElapsedTime(bool reset = false);

 static void useRationals(bool u);
 static bool isUsingRationals();
 static void useFiltering(bool u);
 static bool isUsingFiltering();

 //! Returns the status before the switch
 static bool useRationals() { bool t = isUsingRationals(); useRationals(true); return t; }

 static void setFilename(const char *fname) { filename = fname; }

 static double tri_orientation(double *pa, double *pb, double *pc);
 static double tet_orientation(double *pa, double *pb, double *pc, double *pd);

private:

	static double _spl, _eps, _reb, _ccwebA, _ccwebB, _ccwebC, _o3ebA, _o3ebB, _o3ebC;
	static double _iccebA, _iccebB, _iccebC, _ispebA, _ispebB, _ispebC;

	static int _fesze(int elen, double *e, int flen, double *f, double *h);
	static int _seze(int elen, double *e, double b, double *h);
	static double _estm(int elen, double *e);
	static double _adaptive2dorientation(double *pa, double *pb, double *pc, double detsum);
	static double _adaptive3dorientation(double *pa, double *pb, double *pc, double *pd, double permanent);
};

#define DISPMSG_ACTION_SETWIDGET	1
#define DISPMSG_ACTION_PUTNEWLINE	2
#define DISPMSG_ACTION_PUTPROGRESS	3
#define DISPMSG_ACTION_PUTMESSAGE	4
#define DISPMSG_ACTION_ERRORDIALOG	5

typedef unsigned char	UBYTE;
typedef   signed char	 BYTE;
#ifndef _INC_WINDOWS
typedef unsigned short UINT16;
typedef   signed short	INT16;
#endif

#ifdef IS64BITPLATFORM
typedef long int j_voidint;
#else
typedef int	 j_voidint;
#endif 

#define UBYTE_MAX	255

#ifndef UINT16_MAX
#define UINT16_MAX	65535
#endif

#define FABS(a) (((a)<0)?(-(a)):(a))
#define LOG2(a) (log(a)/log(2))
#define PI2	(M_PI/2.0)
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)(((a)>(b))?(a):(b))
#endif


//////// Swaps two pointers. ///////////////////////////////

inline void p_swap(void **a, void **b) {void *t = *a; *a = *b; *b = t;}

/////////////////////////////////////////////////////////////////////////////////////////////

} //namespace T_MESH

#endif //_BASICS_H

