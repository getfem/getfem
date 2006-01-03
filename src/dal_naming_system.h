// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : File and string TOOL (ftool)
// File    : ftool_naming.h : 
//           
// Date    : August 17, 2002.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2002-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

#ifndef DAL_NAMING_SYSTEM_H
#define DAL_NAMING_SYSTEM_H

#include <stdio.h>
#include <ctype.h>
#include <deque>
#include <map>
#include <dal_static_stored_objects.h>


namespace dal {

  /** @file dal_naming_system.h
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
    dynamic_tree_sorted<std::string> suffixes;
    dynamic_array<pfunction> functions;
    dynamic_array<pgenfunction> genfunctions;
    std::map<std::string, std::string> shorter_names;
    std::map<std::string, std::string> aliases;
    int nb_genfunctions;
    
    struct method_key : virtual public static_stored_object_key {
      std::string name;
      
      virtual bool compare(const static_stored_object_key &oo) const {
	const method_key &o = dynamic_cast<const method_key &>(oo);
	if (name < o.name) return true; else return false;
      }
      method_key(const std::string &name_) : name(name_) {}
    };
    
    int mns_lexem(std::string s, size_type i, size_type &lenght);
    
  public :
    
    void add_suffix(std::string name, pfunction pf);
    void add_generic_function(pgenfunction pf);
    std::string normative_name_of_method(pmethod pm) const;
    std::string shorter_name_of_method(pmethod pm) const;
    pmethod method(std::string name, size_type &i);
    naming_system(std::string pr) : prefix(pr) { nb_genfunctions = 0; }
    
  };
  
  template <class METHOD>
  void naming_system<METHOD>::add_suffix(std::string name,
		       typename naming_system<METHOD>::pfunction pf) {
    size_type i = suffixes.add(prefix + '_' + name);
    functions[i] = pf;
  }

  template <class METHOD>
  void naming_system<METHOD>::add_generic_function(pgenfunction pf) {
    genfunctions[nb_genfunctions++] = pf;
  }
  
  template <class METHOD>
  std::string naming_system<METHOD>::normative_name_of_method(typename 
			         naming_system<METHOD>::pmethod pm)  const {
    pstatic_stored_object_key k = key_of_stored_object(pm);
    if (!k) DAL_THROW(failure_error, "Unknown method");
    return (dynamic_cast<const method_key *>(k))->name;
  }
  
  template <class METHOD> std::string
  naming_system<METHOD>::shorter_name_of_method(typename
			        naming_system<METHOD>::pmethod pm)  const { 
    pstatic_stored_object_key k = key_of_stored_object(pm);
    if (!k)  return "UNKNOWN";
    const std::string &name((dynamic_cast<const method_key *>(k))->name);
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
    DAL_THROW(failure_error, "Invalid character on position " << i
	      << " of the string : " << s);
  }
  
  
  template <class METHOD>
  typename naming_system<METHOD>::pmethod
  naming_system<METHOD>::method(std::string name, size_type &i) {
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
	  ind_suff = suffixes.search(suff);
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
	  pm = method(name, i);
	  params.push_back(parameter(pm));
	  state = 3; break;
	case 3  :
	  char *p;
	  params.push_back(parameter(strtod(&(name[i]), &p)));
	  i += l; if (p < &(name[i])) error = true;
	  state = 3; break;
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
      if (error) 
	DAL_THROW(failure_error,
	     "Syntax error on position " << i << " of the string : " << name);
      if (isend) {
	std::stringstream norm_name;
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
	for (int k = 0; k < nb_genfunctions && pm == 0; ++k) {
	  pm = (*(genfunctions[k]))(nname.name, dependencies);
	}
	if (!pm) {
	  if (ind_suff == size_type(-1))
	    DAL_THROW(failure_error, "Unknown method : " << nname.name);
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
  
}
#endif
