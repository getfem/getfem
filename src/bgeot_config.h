/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_config.h : basic configuration.                        */
/*     									   */
/* Date : December 20, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */



#ifndef __BGEOT_CONFIG_H
#define __BGEOT_CONFIG_H

#include <getfem_arch_config.h>

#include <dal_except.h>

#ifdef GETFEM_HAVE_QDLIB
#  define NO_INLINE
#  ifdef GETFEM_QDLIB_USE_QUAD
#    include <qd.h>
#  else
#    include <dd.h>
#  endif
#  include <x86.h>
#endif

namespace bgeot
{
  static const size_t ST_NIL = size_t(-1);
  typedef dal::uint8_type dim_type;
  typedef dal::uint16_type short_type;
  typedef size_t size_type;
  typedef double scalar_type;
#ifndef GETFEM_HAVE_QDLIB
  typedef double long_scalar_type;
  typedef double opt_long_scalar_type;
# define LONG_SCALAR_EPS 1E-16
# define LONG_SCAL(xx) long_scalar_type(xx)
#else
#  ifdef GETFEM_QDLIB_USE_QUAD
  typedef qd_real long_scalar_type;
  typedef qd_real opt_long_scalar_type;
# define LONG_SCALAR_EPS 1E-64
#  else
  typedef dd_real long_scalar_type;
  typedef dd_real opt_long_scalar_type;
# define LONG_SCALAR_EPS 1E-32
#  endif
#  define LONG_SCAL(xx) long_scalar_type(#xx) /* string assignment to preserve the precision */
#endif
  using dal::dimension_error;
  using dal::file_not_found_error;
  using dal::internal_error;
  using dal::not_linear_error;
  using dal::to_be_done_error;
  using dal::failure_error;

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_CONFIG_H                                                */
