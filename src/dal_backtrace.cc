//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_backtrace.cc : For debugging information
//           
// Date    : June 01, 2003.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1995-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

#include "dal_backtrace.h"
#ifdef GETFEM_HAVE_BACKTRACE
# include <execinfo.h>
#endif
#ifdef GETFEM_HAVE_CXXABI_H
# include <cxxabi.h>
#endif
namespace dal {
  /* demangles a c++ mangled function name, such as the ones
     returned by backtrace_symbols or typeid(x).name()
     If you call this function with a non-mangled name (i.e. "main"),
     you will get strange results.
   */
  std::string demangle(const char *
#ifdef GETFEM_HAVE_CXXABI_H
		       s
#endif
		       ) {
#ifdef GETFEM_HAVE_CXXABI_H
    int status;
    /* documented in http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen/namespaceabi.html */
    char *sd = abi::__cxa_demangle(s, NULL, NULL, &status);
    if (sd == NULL || status) {
      if (sd) free(sd);
      return std::string(""); // + " [could not be demangled]";
    } else {
      std::string res(sd); free(sd); return res;
    }
#else
    return std::string("");
#endif
  }

#ifdef GETFEM_HAVE_BACKTRACE
  void dump_glibc_backtrace() {
    static int cnt = 0;
    int i,n;
    void* trace[256];
    char** strings;
    if (cnt++ == 0) {
      n = backtrace(trace, 256);
      strings = backtrace_symbols (trace, n);
      if (strings == NULL) {
	fprintf(stderr, "backtrace unavailable ... no more memory ?\n"); return;
      }
      fprintf(stderr,"Backtrace dump follows:\n");
      for (i = 0; i < n; ++i) {
	char s[256]; strncpy(s,strings[i],256);s[255]=0;
	char *a = strchr(s,'('), *b = 0;
	if (a) b = strchr(a,'+');
	if (!a || !b) {
	  fprintf(stderr,"%2d : %s\n", i, s);
	} else {
	  *a = 0; *b = 0;
	  fprintf(stderr,"%2d : %s(%s+%s  %s\n",
		  i, s, a+1, b+1, demangle(a+1).c_str());
	}
      }
      free (strings);
      cnt--;
    } else { /* on n'est jamais trop prudent */
      fprintf(stderr, "ooops, a recursive bug in dump_glibc_backtrace\n");
    }
  }
#endif
}
