/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm.h : generic algorithms on linear algebra                 */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2003  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GMM++                                            */
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

#ifndef __GMM_H
#define __GMM_H

#include <dal_std.h>

#if !defined(NOGMM_VERIFY) && defined(GETFEM_VERIFY)
#   define GMM_VERIFY
#endif

#include <gmm_def.h>
#include <gmm_sub_index.h>
#include <gmm_interface.h>
#include <gmm_scaled.h>
#include <gmm_conjugated.h>
#include <gmm_transposed.h>
#include <gmm_blas.h>
#include <gmm_sub_vector.h>
#include <gmm_sub_matrix.h>
#include <gmm_vector.h>
#include <gmm_matrix.h>
#include <gmm_dense_Householder.h>
#include <gmm_dense_lu.h>
#include <gmm_dense_qr.h>
#include <gmm_solvers.h>
#include <gmm_condest.h>
#include <gmm_opt.h>
#include <gmm_lapack_interface.h>

#endif //  __GMM_H
