/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_mesh.h : meshes with points.                           */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
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


#ifndef __BGEOT_MESH_H
#define __BGEOT_MESH_H

#include <bgeot_mesh_structure.h>
#include <bgeot_convex.h>

namespace bgeot
{

  template<class PT> struct mesh_point_ct
          : public dal::dynamic_tree_sorted<PT,
                   dal::lexicographical_less<PT,
                   dal::approx_less<typename PT::value_type> > > { };

  template<class PT, class PT_TAB = mesh_point_ct<PT> > class mesh
    : public mesh_structure
  {
    public :

      typedef PT point_type;
      typedef typename PT::vector_type vector_type;

      typedef dal::tab_ref_index_ref< typename PT_TAB::const_iterator,
                     ref_mesh_point_ind_ct::const_iterator> ref_mesh_pt_ct;
      typedef bgeot::convex<PT, ref_mesh_pt_ct> ref_convex;


    protected :

      dim_type dimension;
      PT_TAB pts;
      
    public :

      const PT_TAB &points(void) const { return pts; }
      PT_TAB &points(void) { return pts; }
      void swap_points(int i, int j)
      { pts.swap(i,j); mesh_structure::swap_points(i,j); }
     
      /* acces aux points directeurs non fait : le faire d'abord dans      */
      /* bgeot_mesh_structure avec container adapte.                       */

      ref_mesh_pt_ct points_of_convex(size_type ic) const 
      {
	ref_mesh_point_ind_ct rct = ind_points_of_convex(ic);
	return ref_mesh_pt_ct(pts.begin(), rct.begin(), rct.end());
      } 

      PT dir_point_of_convex(size_type ic, size_type j) const
      {	return pts[ind_dir_point_of_convex(ic, j)]; }

      ref_convex convex(size_type ic) const
      { return ref_convex(structure_of_convex(ic), points_of_convex(ic)); }

      dim_type dim(void) const { return dimension; }
      void clear(void)
      { dimension = dim_type(-1); mesh_structure::clear(); pts.clear(); }

      mesh(void) { dimension = dim_type(-1); }

      int write_to_file(STD_NEEDED ostream &ost) const;
      int read_from_file(STD_NEEDED istream &ist);
 
  };

  template<class PT, class PT_TAB>
    int mesh<PT, PT_TAB>::read_from_file(STD_NEEDED istream &ist)
  {
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    char tmp[100];
    bool te = false, please_get = true;

    clear();
    ist.seekg(0);
    ftool::read_untill(ist, "BEGIN POINTS LIST");

    while (!te)
    {
      if (please_get) ftool::get_token(ist, tmp, 99); else please_get = true;

      if (!strcmp(tmp, "END"))
      { te = true; }
      else if (!strcmp(tmp, "POINT"))
      {
	ftool::get_token(ist, tmp, 99);
        size_type ip = atoi(tmp);
        dim_type d = 0;
	if (npt.is_in(ip))
	{
	  STD_NEEDED cerr << "BGEOT : Fatal Error, two points with"
	       << "the same index. loading aborted" << endl; return -1;
	}
	npt.add(ip);
	ftool::get_token(ist, tmp, 99);
	while (isdigit(tmp[0]) || tmp[0] == '-' || tmp[0] == '+'
	                       || tmp[0] == '.')
	{ tmpv[d++] = atof(tmp); ftool::get_token(ist, tmp, 99); }
	please_get = false;
	if (dimension == dim_type(-1)) dimension = d; else if (dimension != d)
	{ STD_NEEDED cerr << "BGEOT : Points of different dimensions\n"; return -1; }
	
	PT v(d);
	for (size_type i = 0; i < d; i++) v[i] = tmpv[i];
	points()[ip] = v;
      }
      else
      { STD_NEEDED cerr << "BGEOT : Syntax error in file\n"; return -1; }
    }
    return mesh_structure::read_from_file(ist);
  }

  template<class ITER> void _write_point_to_file(STD_NEEDED ostream &ost, ITER b, ITER e)
  { for ( ; b != e; ++b) ost << "  " << *b; ost << endl; }

  template<class PT, class PT_TAB>
    int mesh<PT, PT_TAB>::write_to_file(STD_NEEDED ostream &ost) const
  {
    ost << endl << "BEGIN POINTS LIST" << endl << endl;
    mesh_point_st_ct::const_iterator b = point_structures().begin();
    mesh_point_st_ct::const_iterator e = point_structures().end();
    for (size_type i = 0; b != e; ++b, ++i)
      if ( (*b).is_valid() )
      {
	ost << "  POINT  " << i;
	_write_point_to_file(ost, points()[i].begin(), points()[i].end());
      }
    ost << endl << "END POINTS LIST" << endl << endl;
    return mesh_structure::write_to_file(ost);
  }

}  /* end of namespace bgeot.                                              */


#endif /* __BGEOT_MESH_H                                                   */
