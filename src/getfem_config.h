/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_config.h : basic configuration.                       */
/*     									   */
/*                                                                         */
/* Date : November 19, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __GETFEM_CONFIG_H
#define __GETFEM_CONFIG_H

#include <dal_std.h>
#include <bgeot_config.h>

/// GEneric Tool for Finite Element Methods.
namespace getfem
{

  static const size_t ST_NIL = size_t(-1);
  typedef size_t           size_type;
  typedef dal::uint8_type  dim_type;
  typedef dal::uint16_type short_type;
  typedef double scalar_type;
  typedef bgeot::base_vector base_vector;
  typedef bgeot::base_matrix base_matrix;
  typedef bgeot::base_tensor base_tensor;
  typedef bgeot::base_poly base_poly;
  typedef bgeot::base_node base_node;

}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_CONFIG_H  */
