#ifndef DAL_BACKTRACE
#define DAL_BACKTRACE

#include <getfem_arch_config.h>
#include <string>
#include <cstdio>

namespace dal {
  std::string demangle(const char *mangled_name);
  inline void dump_backtrace() {
#ifdef GETFEM_HAVE_BACKTRACE
    void dump_glibc_backtrace();
    dump_glibc_backtrace();
#endif
  }
}

#endif
