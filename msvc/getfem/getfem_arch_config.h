#ifndef GETFEM_GETFEM_ARCH_CONFIG_H
#define GETFEM_GETFEM_ARCH_CONFIG_H 1
 
/* Define to 1 if you have the `qhull' library (-lqhull). */
//#define GETFEM_HAVE_LIBQHULL
#define GETFEM_USES_BLAS
#define GETFEM_VERSION "4.0"
#if !defined(__GNUC__)
#include <math.h>
inline double round(double d) { return floor(d + 0.5);}
#endif

#endif
