/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
*/
/*
  Copyright (c) 1997 by Xerox Corporation.  All rights reserved.
 
  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 
  Permission is hereby granted to use or copy this program for any
  purpose, provided the above notices are retained on all copies.
  Permission to modify the code and to distribute modified code is
  granted, provided the above notices are retained, and a notice that
  the code was modified is included with the above copyright notice.
*/
/* 
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */
#ifdef _MSC_VER
#define NO_TIMER
#endif

#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#else

#ifndef NO_TIMER
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>
#endif

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double SuperLU_timer_()
{
#ifdef NO_TIMER
    /* no sys/times.h on WIN32 */
    double tmp;
    tmp = 0.0;
#else
    struct tms use;
    double tmp;
    times(&use);
    tmp = use.tms_utime;
    tmp += use.tms_stime;
#endif
    return (double)(tmp) / CLK_TCK;
}

#endif

