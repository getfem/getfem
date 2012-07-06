/*===========================================================================
 
 Copyright (C) 2005-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

#include "getfem/getfem_mesh_im.h"


namespace getfem {
  
  void mesh_im::update_from_context(void) const {
    for (dal::bv_visitor i(im_convexes); !i.finished(); ++i) {
      if (linked_mesh_->convex_index().is_in(i)) {
	if (v_num_update < linked_mesh_->convex_version_number(i))
	  const_cast<mesh_im *>(this)
	    ->set_integration_method(i, auto_add_elt_pim);
      }
      else const_cast<mesh_im *>(this)->set_integration_method(i, 0);
    }
    for (dal::bv_visitor i(linked_mesh_->convex_index());
	 !i.finished(); ++i) {
      if (!im_convexes.is_in(i)
	  && v_num_update < linked_mesh_->convex_version_number(i)) {
	if (auto_add_elt_pim != 0)
	  const_cast<mesh_im *>(this)
	    ->set_integration_method(i, auto_add_elt_pim);
      }
    }
    v_num_update = v_num = act_counter();
  }


  void mesh_im::set_integration_method(size_type cv, pintegration_method pim) {
    GMM_ASSERT1(linked_mesh_ != 0, "Uninitialized mesh_im");
    context_check();
    if (pim == NULL)
      { if (im_convexes.is_in(cv))
	  { im_convexes.sup(cv); touch(); v_num = act_counter(); } }
    else if (!im_convexes.is_in(cv) || ims[cv] != pim) {
      GMM_ASSERT1
	(pim->type() == IM_NONE || 
	 linked_mesh_->structure_of_convex(cv)->basic_structure()
	 == pim->structure(),
	 "Incompatibility between integration method "
	 << getfem::name_of_int_method(pim) << " and mesh element " <<
	 bgeot::name_of_geometric_trans(linked_mesh_->trans_of_convex(cv)));
      im_convexes.add(cv);
      ims[cv] = pim;
      touch(); v_num = act_counter();
    }
  }

  void mesh_im::set_integration_method(const dal::bit_vector &cvs, 
				       pintegration_method pim) { 
    for (dal::bv_visitor cv(cvs); !cv.finished(); ++cv)
      set_integration_method(cv, pim);
  }

  void mesh_im::set_integration_method(pintegration_method pim) { 
    set_integration_method(linked_mesh().convex_index(), pim);
    set_auto_add(pim);
  }
  
  void mesh_im::set_integration_method(const dal::bit_vector &cvs, 
				       dim_type im_degree) {
    GMM_ASSERT1(im_degree != dim_type(-1), "im_degree==-1");
    for (dal::bv_visitor cv(cvs); !cv.finished(); ++cv) {
      pintegration_method pim = 
	getfem::classical_approx_im(linked_mesh().trans_of_convex(cv), im_degree);
      set_integration_method(cv, pim);
    }
  }

  void mesh_im::clear(void) {
    ims.clear(); im_convexes.clear();
    touch(); v_num = act_counter();
  }


  void mesh_im::init_with_mesh(mesh &me) {
    GMM_ASSERT1(linked_mesh_ == 0, "Mesh level set already initialized");
    linked_mesh_ = &me;
    this->add_dependency(me);
    auto_add_elt_pim = 0;
    v_num_update = v_num = act_counter();
  }
  
  mesh_im::mesh_im(void) {
    linked_mesh_ = 0; auto_add_elt_pim = 0;
    v_num_update = v_num = act_counter();
  }

  mesh_im::mesh_im(const mesh_im &mim) {
    GMM_ASSERT1(mim.linked_mesh_ == 0,
		"Copy constructor is not allowed for non void mesh_im");
    linked_mesh_ = 0; auto_add_elt_pim = 0;
    v_num_update = v_num = act_counter();
  }

  mesh_im & mesh_im::operator=(const mesh_im &mim) {
    GMM_ASSERT1(linked_mesh_ == 0 && mim.linked_mesh_ == 0,
		"Copy operator is not allowed for non void mesh_im");
    return *this;
  }

  mesh_im::mesh_im(mesh &me)
  { linked_mesh_ = 0; init_with_mesh(me); }

  mesh_im::~mesh_im() {}


  void mesh_im::read_from_file(std::istream &ist) {
    GMM_ASSERT1(linked_mesh_ != 0, "Uninitialized mesh_im");
    gmm::stream_standard_locale sl(ist);
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    std::string tmp;
    ist.precision(16);
    clear();
    ist.seekg(0);ist.clear();
    bgeot::read_until(ist, "BEGIN MESH_IM");

    while (true)
    {
      ist >> std::ws; bgeot::get_token(ist, tmp);
      if (bgeot::casecmp(tmp, "END")==0) {
	break;
      } else if (bgeot::casecmp(tmp, "CONVEX")==0) {
	bgeot::get_token(ist, tmp);
	size_type ic = atoi(tmp.c_str());
	GMM_ASSERT1(linked_mesh().convex_index().is_in(ic), "Convex " << ic <<
		    " does not exist, are you sure "
		    "that the mesh attached to this object is right one ?");
	
	int rgt = bgeot::get_token(ist, tmp);
	if (rgt != 3) { // for backward compatibility with version 1.7
	  char c; ist.get(c);
	  while (!isspace(c)) { tmp.push_back(c); ist.get(c); }
	}
	getfem::pintegration_method pfi = getfem::int_method_descriptor(tmp);
	GMM_ASSERT1(pfi, "could not create the integration method '"
		    << tmp << "'");
	
	set_integration_method(ic, pfi);
      } else if (tmp.size()) {
	GMM_ASSERT1(false, "Unexpected token '" << tmp <<
		  "' [pos=" << std::streamoff(ist.tellg()) << "]");
      } else if (ist.eof()) {
	GMM_ASSERT1(false, "Unexpected end of stream "
		    << "(missing BEGIN MESH_IM/END MESH_IM ?)");	
      }
    }
  }

  void mesh_im::read_from_file(const std::string &name)
  { 
    std::ifstream o(name.c_str());
    GMM_ASSERT1(o, "mesh_im file '" << name << "' does not exist");
    read_from_file(o);
    o.close();
  }

  void mesh_im::write_to_file(std::ostream &ost) const {
    context_check();
    gmm::stream_standard_locale sl(ost);
    ost << '\n' << "BEGIN MESH_IM" << '\n' << '\n';
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      ost << " CONVEX " << cv;
      ost << " \'" << name_of_int_method(int_method_of_element(cv));
      ost << "\'\n";
    }

    ost << "END MESH_IM" << '\n';
  }

  void mesh_im::write_to_file(const std::string &name, bool with_mesh) const
  {
    std::ofstream o(name.c_str());
    GMM_ASSERT1(o, "impossible to open file '" << name << "'");
    o << "% GETFEM MESH_IM FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    if (with_mesh) linked_mesh().write_to_file(o);
    write_to_file(o);
    o.close();
  }
}  /* end of namespace getfem.                                             */



