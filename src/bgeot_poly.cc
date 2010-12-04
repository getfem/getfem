// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2008 Yves Renard
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



#include "getfem/bgeot_poly.h"
#include "getfem/bgeot_vector.h"
#include "getfem/bgeot_ftool.h"

namespace bgeot {

#define STORED 150
  static gmm::dense_matrix<size_type> alpha_M_(STORED, STORED);
  static void alpha_init_() {
    static bool init = false;
    if (!init) {
      for (short_type i = 0; i < STORED; ++i) {
	alpha_M_(i, 0) = alpha_M_(0, i) = 1;
	for (short_type j = 1; j <= i; ++j)
	  alpha_M_(i,j) = alpha_M_(j,i) = (alpha_M_(i, j-1) * (i+j)) / j;
      }
      init = true;
    }
  }
  static inline size_type alpha_(short_type n, short_type d)
  { return alpha_M_(d,n); }

  size_type alpha(short_type n, short_type d) {
    alpha_init_();
    GMM_ASSERT1(n < STORED && d < STORED,
		"alpha called with n = " << n << " and d = " << d);
    return alpha_(n,d);
  }

  const power_index &power_index::operator ++() {
    short_type n = short_type(size()), l;
    if (n > 0) {
      size_type g_idx = global_index_; short_type deg = degree_;
      reverse_iterator it = rbegin() + 1;
      for (l = short_type(n-2); l != short_type(-1); --l, ++it)
	if (*it != 0) break;
      short_type a = (*this)[n-1]; (*this)[n-1] = 0;
      (*this)[short_type(l+1)] = short_type(a + 1);
      if (l != short_type(-1)) ((*this)[l])--;
      else if (short_type(deg+1)) degree_ = short_type(deg+1);
      if (g_idx+1) global_index_ = g_idx+1;
      //degree_ = short_type(-1);
    }
    return *this;
  }
  
  const power_index &power_index::operator --() {
    short_type n = short_type(size()), l;
    if (n > 0) {
      size_type g_idx = global_index_; short_type deg = degree_;
      reverse_iterator it = rbegin();
      for (l = short_type(n-1); l != short_type(-1); --l, ++it) {
	if (*it != 0) break;
      }
      if (l != short_type(-1)) {
	short_type a = (*this)[l];
	(*this)[l] = 0; (*this)[n-1] = short_type(a - 1);
	if (l > 0) ((*this)[l-1])++; 
        else if (short_type(deg+1)) degree_ = short_type(deg-1);
      }
      if (g_idx+1) global_index_ = g_idx-1;
    }
    return *this;
  }
  
  short_type power_index::degree() const {
    if (degree_ != short_type(-1)) return degree_;
    degree_ = short_type(std::accumulate(begin(), end(), 0)); 
    return degree_;
  }

  size_type power_index::global_index(void) const {
    if (global_index_ != size_type(-1)) return global_index_;
    short_type d = degree(), n = short_type(size());
    global_index_ = 0;
    const_iterator it = begin(), ite = end();
    for ( ; it != ite && d > 0; ++it)
      { global_index_ += alpha_(n,short_type(d-1)); d=short_type(d-*it); --n; }
    return global_index_;
  }
 
  power_index::power_index(short_type nn) : v(nn), degree_(0), global_index_(0)
  { std::fill(begin(), end(), short_type(0)); alpha_init_(); }


  // functions to read a polynomial on a stream

  static void parse_error(int i)
  { GMM_ASSERT1(false, "Syntax error reading a polynomial " << i); }

  static std::string stored_s;
  int stored_tokent;

  static void unget_token(int i, std::string s)
  { stored_s = s; stored_tokent = i; }

  static int get_next_token(std::string &s, std::istream &f) {
    if (stored_s.size() == 0)
      return get_token(f, s, true, false, false);
    else { s = stored_s; stored_s.clear(); return stored_tokent; }
  }

  static base_poly read_expression(short_type n, std::istream &f) {
    gmm::stream_standard_locale sl(f);
    gmm::standard_locale sll;
    
    base_poly result(n,0);
    std::string s;
    int i = get_next_token(s, f), j;
    switch (i) {
    case 2 : result.one();
      result *= opt_long_scalar_type(::strtod(s.c_str(), 0));
      break;
    case 4 :
      if (s == "x") result = base_poly(n, 1, 0);
      else if (s == "y" && n > 1) result = base_poly(n, 1, 1);
      else if (s == "z" && n > 2) result = base_poly(n, 1, 2);
      else if (s == "w" && n > 3) result = base_poly(n, 1, 3);
      else if (s == "v" && n > 4) result = base_poly(n, 1, 4);
      else if (s == "u" && n > 5) result = base_poly(n, 1, 5);
      else if (s == "t" && n > 6) result = base_poly(n, 1, 6);
      else if (s == "sqrt") {
	base_poly p = read_expression(n, f);
	if (p.degree() > 0) parse_error(1);
	result.one();  result *= sqrt(p[0]);
      }
      else parse_error(2);
      break;
    case 5 :
      switch (s[0]) {
      case '(' :
	result = read_base_poly(n, f);
	j = get_next_token(s, f);
	if (j != 5 || s[0] != ')') parse_error(3);
	break;
	default : parse_error(4);
      }
      break;
    default : parse_error(5);
    }
    return result;
  }

  static void operator_priority_(int i, char c, int &prior, int &op) {
    if (i == 5)
      switch (c) {
      case '*' : prior = 2; op = 1; return;
      case '/' : prior = 2; op = 2; return;
      case '+' : prior = 3; op = 3; return;
      case '-' : prior = 3; op = 4; return;
      case '^' : prior = 1; op = 5; return;
      }
    prior = op = 0;
  }

  void do_bin_op(std::vector<base_poly> &value_list,
		 std::vector<int> &op_list,
		 std::vector<int> &prior_list) {
    base_poly &p2 = *(value_list.rbegin());
    if (op_list.back() != 6) {
      assert(value_list.size()>1);
      base_poly &p1 = *(value_list.rbegin()+1);    
      switch (op_list.back()) {
        case 1  : p1 *= p2; break;
        case 2  : if (p2.degree() > 0) parse_error(6); p1 /= p2[0]; break;
        case 3  : p1 += p2; break;
        case 4  : p1 -= p2; break;
        case 5  : 
	  {
	    if (p2.degree() > 0) parse_error(7);
	    int pow = int(p2[0]);
	    if (p2[0] !=  opt_long_scalar_type(pow) || pow < 0) parse_error(8);
	    base_poly p = p1; p1.one();
	    for (int i = 0; i < pow; ++i) p1 *= p;
	  }
	  break;
        default: assert(0);
      }
      value_list.pop_back(); 
    } else {
      p2 *= opt_long_scalar_type(-1);
    }
    op_list.pop_back(); prior_list.pop_back();
  }

  base_poly read_base_poly(short_type n, std::istream &f) {
    std::vector<base_poly> value_list;
    std::string s;
    std::vector<int> op_list, prior_list;
    
    int i = get_next_token(s, f), prior, op;
    if (i == 5 && s[0] == '-')
      { op_list.push_back(6); prior_list.push_back(2); }
    else if (i == 5 && s[0] == '+') ;
    else unget_token(i, s);

    value_list.push_back(read_expression(n, f));
    i = get_next_token(s, f);
    operator_priority_(i, s[0], prior, op);
    while (op) {
      while (!prior_list.empty() && prior_list.back() <= prior)
	do_bin_op(value_list, op_list, prior_list);

      value_list.push_back(read_expression(n, f));
      op_list.push_back(op);
      prior_list.push_back(prior);
      
      i = get_next_token(s, f);
      operator_priority_(i, s[0], prior, op);
    }
    
    if (i == 5 && s[0] == ')') { f.putback(')'); }
    else if (i != 0 && (i != 5 || s[0] != ';')) {
      cout << "s = " << s << endl;
      parse_error(9);
    }

    while (!prior_list.empty()) do_bin_op(value_list, op_list, prior_list);

    return value_list[0];
  }

  base_poly read_base_poly(short_type n, const std::string &s)
  { std::stringstream f(s); return read_base_poly(n, f); }


}  /* end of namespace bgeot.                                             */
