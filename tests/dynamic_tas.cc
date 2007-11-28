// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================
#include "getfem/dal_tas.h"

using namespace std;

int main(void)
{
  try {
    dal::dynamic_tas<int> t;
    t.add(6);
    // cout << "index = " << t.index() << endl;
    // cout << "first free place : " << t.index().first_false() << endl;
    if (t.index().first_false() != 1)
      throw gmm::internal_error("dynamic_tas.C : structure error 1");
    
    t.add(4);
    t.add(2);
    
    // cout << "index = " << t.index() << endl;
    
    t.sup(1);
    // cout << "index = " << t.index() << endl;
    // cout << "first free place : " << t.index().first_false() << endl;
    if (t.index().first_false() != 1)
      throw gmm::internal_error("dynamic_tas.C : structure error 2");
    
    t.add(3);
    
    dal::dynamic_tas<int>::iterator it = t.begin(), end = t.end();
    // cout << " card = " << t.card() << endl;
    // cout << "index = " << t.index() << endl;
    // cout << "index.first() = " << t.index().first() << endl;
    
    for ( ; it != end; ++it) { cout << " : " << *it; }
    cout << endl;
    
    return 0;
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;
}
