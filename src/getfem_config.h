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

#include <bgeot_tensor.h>
#include <bgeot_poly.h>

/// GEneric Tool for Finite Element Methods.
namespace getfem
{
  using bgeot::ST_NIL;
  using bgeot::size_type;
  using bgeot::dim_type;
  using bgeot::short_type;
  using bgeot::short_type;
  using bgeot::scalar_type;
  
  using bgeot::base_vector;
  using bgeot::base_matrix;
  using bgeot::base_tensor;
  using bgeot::base_poly;
  using bgeot::base_node;

  using dal::dimension_error;
  using dal::internal_error;
  using dal::not_linear_error;
  using dal::to_be_done_error;
  using dal::failure_error;

}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_CONFIG_H  */
