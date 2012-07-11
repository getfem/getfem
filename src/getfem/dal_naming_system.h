/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard
 
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

#ifndef DAL_NAMING_SYSTEM_H
#define DAL_NAMING_SYSTEM_H

#include <deque>
#include <map>
#include "dal_static_stored_objects.h"


namespace dal {

  /** @file dal_naming_system.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date August 17, 2002.
      @brief Naming system.
  */


  /** Associate a name to a method descriptor and store method descriptors.
   *
   * Methods may have parameters such as integer or other methods.
   *  The class METHOD have to derive from dal::static_stored_object
   */
  template <class METHOD> class naming_system {

  public :

    typedef boost::intrusive_ptr<const METHOD> pmethod;

    struct parameter {
      int type_; // 0 = numeric value, 1 = pointer on another method.
      double num_;
      pmethod pm_;

      pmethod method(void) const { return pm_; }
      double num(void) const { return num_; }
      int type(void) const { return type_; }
      parameter(double e) : type_(0), num_(e), pm_(0) {}
      parameter(pmethod p) : type_(1), num_(0.), pm_(p) {}
    };

    typedef std::deque<parameter> param_list;
    typedef pmethod (* pfunction)(param_list &,
				  std::vector<pstatic_stored_object> &);
    typedef pmethod (* pgenfunction)(std::string,
				     std::vector<pstatic_stored_object> &);
    typedef size_t size_type;

  protected :

    std::string prefix;
    std::map<std::string, size_type> suffixes;
    std::vector<pfunction> functions;
    std::vector<pgenfunction> genfunctions;
    std::map<std::string, std::string> shorter_names;
    std::map<std::string, std::string> aliases;

    struct method_key : virtual public static_stored_object_key {
      std::string name;

      virtual bool compare(const static_stored_object_key &oo) const {
	const method_key &o = dynamic_cast<const method_key &>(oo);
	if (name < o.name) return true; else return false;
      }
      method_key(const std::string &name_) : name(name_) {}
    };

    int mns_lexem(std::string s, size_type i, size_type &lenght);
    pmethod method_(std::string name, size_type &i, bool throw_if_not_found);




  public :

    void add_suffix(std::string name, pfunction pf);
    void add_generic_function(pgenfunction pf);
    std::string normative_name_of_method(pmethod pm) const;
    std::string shorter_name_of_method(pmethod pm) const;
    pmethod method(std::string name, size_type &i,
		   bool throw_if_not_found = true)
    { gmm::standard_locale sl; return method_(name, i, throw_if_not_found); }
    naming_system(std::string pr) : prefix(pr) {}
    bool delete_method(std::string name);
  };

  template <class METHOD>
  void naming_system<METHOD>::add_suffix(std::string name,
		       typename naming_system<METHOD>::pfunction pf) {
    std::string tname = prefix + '_' + name;
    if (suffixes.find(tname) != suffixes.end()) {
      functions[suffixes[tname]] = pf;
    } else {
      suffixes[tname] = functions.size();
      functions.push_back(pf);
    }

  }

  template <class METHOD>
  void naming_system<METHOD>::add_generic_function(pgenfunction pf) {
    genfunctions.push_back(pf);
  }

  template <class METHOD>
  std::string naming_system<METHOD>::normative_name_of_method(typename
			         naming_system<METHOD>::pmethod pm)  const {
    pstatic_stored_object_key k = key_of_stored_object(pm);
    const method_key *p;
    if (!k || !(p = dynamic_cast<const method_key *>(k)))
      return prefix + "_UNKNOWN";
    return p->name;
  }

  template <class METHOD> std::string
  naming_system<METHOD>::shorter_name_of_method(typename
			        naming_system<METHOD>::pmethod pm)  const {
    pstatic_stored_object_key k = key_of_stored_object(pm);
    const method_key *p;
    if (!k || !(p = dynamic_cast<const method_key *>(k)))
      return prefix + "_UNKNOWN";
    const std::string &name(p->name);
    std::map<std::string, std::string>::const_iterator
      it = shorter_names.find(name);
    if (it != shorter_names.end()) return it->second;
    return name;
  }

  /* 0 = end of the string
     1 = espace
     2 = method name
     3 = number
     4 = '('
     5 = ')'
     6 = ','
  */
  template <class METHOD>
  int naming_system<METHOD>::mns_lexem(std::string s, size_type i,
				       size_type &lenght) {
    lenght = 1;
    if (i >= s.size()) return 0;
    char c = s[i];
    if (isspace(c)) return 1;
    if (isalpha(c) || c == '_') {
      for (c = s[++i] ; isalpha(c) || c == '_' || isdigit(c); c = s[++i])
	++lenght;
      return 2;
    }
    if (isdigit(c) || c == '-' || c == '+') {
      for (c = s[++i] ; isdigit(c) || c == 'e' || c == 'E' ||
	     c == '.' || c == '-' || c == '+' ; c = s[++i]) ++lenght;
      return 3;
    }
    if (c == '(') return 4;
    if (c == ')') return 5;
    if (c == ',') return 6;
    GMM_ASSERT1(false, "Invalid character on position " << i
		<< " of the string : " << s);
  }


  template <class METHOD>
  typename naming_system<METHOD>::pmethod
  naming_system<METHOD>::method_(std::string name, size_type &i,
				 bool throw_if_not_found) {
    int state = 0;
    bool error = false;
    bool isend = false;
    pmethod pm = 0;
    size_type ind_suff = size_type(-1);
    size_type l;
    param_list params;
    std::string suff;

    for(;;) {
      int lex = mns_lexem(name, i, l);
      switch (state) {
      case 0 :
	switch (lex) {
	case 1  : i += l; break;
	case 2  :
	  suff = name.substr(i, l);
	  if (suffixes.find(suff) != suffixes.end())
	    ind_suff = suffixes[suff];
	  state = 1; i += l; break;
	default : error = true;
	}
	break;
      case 1 :
	switch (lex) {
	case 4  : state = 2; i += l; break;
	default : isend = true; break;
	}
	break;
      case 2 :
	switch (lex) {
	case 1  : i += l; break;
	case 2  :
	  pm = method_(name, i, throw_if_not_found);
	  if (!pm) return pm;
	  params.push_back(parameter(pm));
	  state = 3; break;
	case 3  : {
	  char *p;
	  gmm::standard_locale sl;
	  params.push_back(parameter(strtod(&(name[i]), &p)));
	  i += l; if (p < &(name[i])) error = true;
	  state = 3; break;
	}
	case 5  : i += l; isend = true; break;
	default : error = true;
	}
	break;
      case 3 :
	switch (lex) {
	case 1  : i += l; break;
	case 5  : i += l; isend = true; break;
	case 6  : i += l; state = 2; break;
	default : error = true;
	}
	break;
      }
      GMM_ASSERT1(!error, "Syntax error on position " << i
		  << " of the string : " << name);
      if (isend) {
	std::stringstream norm_name; norm_name.imbue(std::locale("C"));
	norm_name << suff;
	if (params.size() > 0) {
	  norm_name << '(';
	  typename param_list::const_iterator it = params.begin(),
	    ite = params.end();
	  for (; it != ite; ++it) {
	    if ((*it).type() == 0) norm_name << (*it).num();
	    if ((*it).type() == 1)
	      norm_name << normative_name_of_method((*it).method());
	    if (it+1 != ite) norm_name << ',';
	  }
	  norm_name << ')';
	}
	method_key nname(norm_name.str());
	if (aliases.find(norm_name.str()) != aliases.end())
	  nname.name = aliases[norm_name.str()];
	pstatic_stored_object o = search_stored_object(nname);
	if (o) return stored_cast<METHOD>(o);
	pm = 0;
	std::vector<pstatic_stored_object> dependencies;
	for (size_type k = 0; k < genfunctions.size() && pm == 0; ++k) {
	  pm = (*(genfunctions[k]))(nname.name, dependencies);
	}
	if (!pm) {
	  if (ind_suff == size_type(-1)) {
	    GMM_ASSERT1(!throw_if_not_found, "Unknown method: " << nname.name);
	    return 0;
	  }
	  pm = (*(functions[ind_suff]))(params, dependencies);
	}
	
	pstatic_stored_object_key k = key_of_stored_object(pm);
	if (!k) {
	  add_stored_object(new method_key(nname), pm,
			    dal::PERMANENT_STATIC_OBJECT);
	  for (size_type j = 0; j < dependencies.size(); ++j)
	    add_dependency(pm, dependencies[j]);
	}
	else {
	  std::string normname((dynamic_cast<const method_key *>(k))->name);
	  aliases[nname.name] = normname;
	  if (nname.name.size() < normname.size()) {
	    if (shorter_names.find(normname) != shorter_names.end()) {
	      if (nname.name.size() < shorter_names[normname].size())
		shorter_names[normname] = nname.name;
	    }
	    else shorter_names[normname] = nname.name;
	  }
	}
	return pm;
      }
    }

  }

  /**deletion of static_stored_object in the naming system*/
  template <class METHOD>
  bool naming_system<METHOD>::delete_method(std::string name) 
  {

    pmethod pm = 0;

	method_key nname(name);
	pstatic_stored_object o = search_stored_object(nname);

	if (!o) return false;
	pm = stored_cast<METHOD>(o);
	pstatic_stored_object_key k = key_of_stored_object(pm);
	dal::del_stored_object(pm, false);
	return true;
  }

}
#endif
