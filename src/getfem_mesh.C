/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_mesh.C : meshes for computations.                     */
/*     									   */
/*                                                                         */
/* Date : November 05, 1999.                                               */
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


#include <getfem_mesh.h>

namespace getfem
{
  getfem_mesh::getfem_mesh(dim_type NN)
  {
    dimension = NN; eps_p = 1.0E-10;
    pts.comparator() = dal::lexicographical_less<base_node,
                dal::approx_less<base_node::value_type> >(1.0E-10);
  }

  size_type getfem_mesh::add_point(const base_node &pt)
  {
    if (dimension == dim_type(-1)) dimension = pt.size();
    if (pt.size() != dimension)
      throw dimension_error("getfem_mesh::add_point : dimensions mismatch");
    
    bool present;
    size_type i = pts.add_norepeat(pt, false, &present);
    if (!present) lmsg_sender().send(MESH_ADD_POINT(i));
    return i;
  }

  void getfem_mesh::sup_point(size_type i)
  {
    if (!point_is_valid(i))
    { pts.sup(i); lmsg_sender().send(MESH_SUP_POINT(i)); }
  }

  void getfem_mesh::swap_points(size_type i, size_type j)
  {
    if (i != j)
    {
      bgeot::mesh<base_node>::swap_points(i,j);
      lmsg_sender().send(MESH_SWAP_POINT(i, j));
    }	 
  }

  void getfem_mesh::optimize_structure()
  {
    size_type i, j;
    j = nb_convex();
    for (i = 0; i < j; i++)
      if (!convex_tab.index_valid(i))
	swap_convex(i, convex_tab.ind_last());

    for (i = 0, j = (points_tab.end()-points_tab.begin())-1; i < j && j != ST_NIL; ++i, --j)
    {
      while (i < j && j != ST_NIL && points_tab[i].first != ST_NIL) ++i;
      while (i < j && j != ST_NIL && points_tab[j].first == ST_NIL) --j;
      if (i < j && j != ST_NIL ) swap_points(i, j);
    }
  }

  void getfem_mesh::translation(base_vector V)
  {
    dal::bit_vector nn = points().index();
    size_type i;
    for (i << nn; i != ST_NIL; i << nn)
      points()[i] += V;
    points().resort();
  }

  void getfem_mesh::transformation(base_matrix M)
  {
    dal::bit_vector nn = points().index();
    size_type i;
    for (i << nn; i != ST_NIL; i << nn)
      points()[i] *= M;
    points().resort();
  }

  void getfem_mesh::clear(void)
  { 
    bgeot::mesh<base_node>::clear();
    gtab.clear(); trans_exists.clear();
    lmsg_sender().send(MESH_CLEAR());
  }

  size_type getfem_mesh::add_segment(size_type a, size_type b)
  { 
    static bgeot::pgeometric_trans cs = bgeot::simplex_trans(1, 1);
    size_type ipt[2]; ipt[0] = a; ipt[1] = b;
    return add_convex(cs, &(ipt[0]));
  }
  
  size_type getfem_mesh::add_triangle(size_type a, 
					   size_type b, size_type c)
  {
    size_type ipt[3]; ipt[0] = a; ipt[1] = b; ipt[2] = c;
    return add_simplex(2, &(ipt[0]));
  }

  size_type getfem_mesh::add_triangle_by_points
    (const base_node &p1, const base_node &p2, const base_node &p3)
  { return add_triangle(add_point(p1), add_point(p2), add_point(p3)); }
  
  size_type getfem_mesh::add_tetrahedron(size_type a, 
				       size_type b, size_type c, size_type d)
  {
    size_type ipt[4]; ipt[0] = a; ipt[1] = b; ipt[2] = c; ipt[3] = d;
    return add_simplex(3, &(ipt[0]));
  }

  size_type getfem_mesh::add_tetrahedron_by_points
    (const base_node &p1, const base_node &p2, const base_node &p3, const base_node &p4)
  {
    return add_tetrahedron(add_point(p1), add_point(p2),
			   add_point(p3), add_point(p4));
  }

  void getfem_mesh::sup_convex(size_type ic)
  {
    bgeot::mesh<base_node>::sup_convex(ic);
    trans_exists[ic] = false;
    lmsg_sender().send(MESH_SUP_CONVEX(ic));
  }

  void getfem_mesh::swap_convex(size_type i, size_type j) {
    if (i != j) {
      bgeot::mesh<base_node>::swap_convex(i,j);
      trans_exists.swap(i, j);
      gtab.swap(i,j);
      lmsg_sender().send(MESH_SWAP_CONVEX(i, j));
    }
  }

  base_vector getfem_mesh::normal_of_face_of_convex(size_type ic, short_type f,
						 const base_node &pt) const {
    size_type N  = dim();
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix S(N,pgt->nb_points());
    
    for (size_type i=0; i < pgt->nb_points(); i++) {
      std::copy(points_of_convex(ic)[i].begin(), 
		points_of_convex(ic)[i].end(), S.begin()+i*N);
    }
    
    return bgeot::compute_normal(S, f, pgt, pt);
  }


  base_vector getfem_mesh::normal_of_face_of_convex(size_type ic, short_type f,
						    size_type n) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_node pt = pgt->geometric_nodes()
      [pgt->structure()->ind_points_of_face(f)[n]];
    return normal_of_face_of_convex(ic, f, pt);
  }

  int getfem_mesh::read_from_file(std::istream &ist) {
   
    int r = bgeot::mesh<base_node>::read_from_file(ist);
    dal::bit_vector nn = convex_index();
    size_type i;
    for (i << nn; i != size_type(-1); i << nn)
      if (!(trans_exists[i])) {
	gtab[i] = bgeot::associated_trans(structure_of_convex(i));
	trans_exists[i] = true;
      }

    lmsg_sender().send(MESH_READ_FROM_FILE(ist));
    return r;
  }

  int getfem_mesh::read_from_file(const std::string &name)
  { 
    std::ifstream o(name.data());
    if (!o) DAL_THROW(std::invalid_argument, "Mesh file does not exist");
    return read_from_file(o);
  }

  int getfem_mesh::write_to_file(std::ostream &ost) const
  {
    bgeot::mesh<base_node>::write_to_file(ost);
    lmsg_sender().send(MESH_WRITE_TO_FILE(ost));
    return 0;
  }

  int getfem_mesh::write_to_file(const std::string &name) const
  {
    std::ofstream o(name.data());
    if (o)
    {
      o << "% GETFEM MESH FILE " << endl;
      o << "% GETFEM VERSION " << __GETFEM_VERSION << "."
	                      << __GETFEM_REVISION << endl;
      o << "% Yves.Renard@gmm.insa-tlse.fr " << endl
	<< endl << endl;
    
      return write_to_file(o);
    }
    return -1;
  }

}  /* end of namespace getfem.                                             */
