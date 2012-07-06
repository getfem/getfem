/*===========================================================================
 
 Copyright (C) 2007-2012 Julien Pommier.
 
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
#include <getfemint_levelset.h>
#include <getfemint_levelset.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_arch_config.h>

#if GETFEM_HAVE_MUPARSER_MUPARSER_H
#include <muParser/muParser.h>
#elif GETFEM_HAVE_MUPARSER_H
#include <muParser.h>
#endif

namespace getfemint {
  getfemint_levelset* getfemint_levelset::get_from(getfem::level_set *ls,
                                                   int flags) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(ls));
    getfemint_levelset *gls = 0;
    if (!o) {
      const getfem::mesh &m = ls->get_mesh_fem().linked_mesh();
      getfemint_mesh *gm =
        getfemint_mesh::get_from(const_cast<getfem::mesh*>(&m),
                                 flags);
      gls = new getfemint_levelset();
      gls->ls = ls;
      gls->ikey = getfem_object::internal_key_type(ls);
      gls->set_flags(flags);
      getfemint::workspace().push_object(gls);
      getfemint::workspace().set_dependance(gls, gm);
    } else gls = dynamic_cast<getfemint_levelset*>(o);
    assert(gls);
    return gls;

  }

  void getfemint_levelset::values_from_poly(unsigned idx,
                                            const std::string &s) {
    const getfem::mesh_fem &mf = levelset().get_mesh_fem();
    assert(!mf.is_reduced());
    bgeot::base_poly p = bgeot::read_base_poly(mf.linked_mesh().dim(), s);
    ls->values(idx).resize(mf.nb_dof());
    for (unsigned i=0; i < mf.nb_dof(); ++i) {
      const getfem::base_node x = mf.point_of_basic_dof(i);
      ls->values(idx)[i] = p.eval(x.begin());
    }
  }
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
  void getfemint_levelset::values_from_func(unsigned idx,
                                            const std::string &s) {
    const getfem::mesh_fem &mf = levelset().get_mesh_fem();
    assert(!mf.is_reduced());

    double* x = (double *)calloc(mf.linked_mesh().dim(), sizeof(double));

    mu::Parser p;
    try {
      switch(mf.linked_mesh().dim()) {
        case 1:
         p.DefineVar("x",&x[0]); break;
        case 2:
         p.DefineVar("x",&x[0]);
         p.DefineVar("y",&x[1]); break;
        case 3:
         p.DefineVar("x",&x[0]);
         p.DefineVar("y",&x[1]);
         p.DefineVar("z",&x[2]); break;
      }
      p.SetExpr(s);
    } catch (mu::Parser::exception_type &e) {
      std::cout << "Message  : " << e.GetMsg() << std::endl;
      std::cout << "Formula  : " << e.GetExpr() << std::endl;
      std::cout << "Token    : " << e.GetToken() << std::endl;
      std::cout << "Position : " << e.GetPos() << std::endl;
      std::cout << "Errc     : " << e.GetCode() << std::endl;
    }

    ls->values(idx).resize(mf.nb_dof());
    bool is_set = 0;
    for (unsigned i=0; i < mf.nb_dof(); ++i) {
      is_set = 0;
      switch(mf.linked_mesh().dim()) {
        case 1 :
         x[0] = mf.point_of_basic_dof(i)[0];
         is_set = 1;
         break;
        case 2 :
         x[0] = mf.point_of_basic_dof(i)[0];
         x[1] = mf.point_of_basic_dof(i)[1];
         is_set = 1;
         break;
        case 3 :
         x[0] = mf.point_of_basic_dof(i)[0];
         x[1] = mf.point_of_basic_dof(i)[1];
         x[2] = mf.point_of_basic_dof(i)[2];
         is_set = 1;
         break;
      }
      try {
        if (is_set) ls->values(idx)[i] = p.Eval();
      } catch (mu::Parser::exception_type &e) {
        std::cout << "Message  : " << e.GetMsg() << std::endl;
        std::cout << "Formula  : " << e.GetExpr() << std::endl;
        std::cout << "Token    : " << e.GetToken() << std::endl;
        std::cout << "Position : " << e.GetPos() << std::endl;
        std::cout << "Errc     : " << e.GetCode() << std::endl;
      }
    }
    free(x);
  }
#endif
}
