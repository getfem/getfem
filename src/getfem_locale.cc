/*===========================================================================

 Copyright (C) 2012-2017 Andriy Andreykiv.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <iostream>

#include <getfem/getfem_locale.h>
#include <getfem/getfem_omp.h>

namespace getfem{

  #ifdef _WIN32

    standard_locale::standard_locale()
      : cinloc(std::cin.getloc()){
      if (!me_is_multithreaded_now()){
          cloc=setlocale(LC_NUMERIC, 0);
          setlocale(LC_NUMERIC,"C");
      }
    }

    standard_locale::~standard_locale() {
      if (!me_is_multithreaded_now())
          setlocale(LC_NUMERIC, cloc.c_str());

    }

  #else

    standard_locale::standard_locale()
      : cloc(setlocale(LC_NUMERIC, 0)), cinloc(cin.getloc()){
      setlocale(LC_NUMERIC,"C"); cin.imbue(std::locale("C"));
    }

    standard_locale::~standard_locale(){
      setlocale(LC_NUMERIC, cloc.c_str()); cin.imbue(cinloc);
    }

  #endif

}  /* end of namespace getfem.                                             */