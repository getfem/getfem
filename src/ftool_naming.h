/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  File and string TOOL (ftool)                                 */
/* File    :  ftool_naming.h :                                             */
/*     									   */
/*                                                                         */
/* Date : August 17, 2002.                                                 */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
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

#include <stdio.h>
#include <dal_basic.h>
#include <ctype.h>
#include <deque>


namespace ftool
{
  /* ********************************************************************* */
  /*                                                                       */
  /*   Naming system                                                       */
  /*                                                                       */
  /* ********************************************************************* */

  /** This class associates a name to a method descriptor and store
   *  method descriptors. Methods may have parameters such as integer or
   *  other methods.
   */
  template <class METHOD> class naming_system {

  public :

    typedef const METHOD *pmethod;
    
    struct parameter {
      int _type; // 0 = numerique, 1 = pointeur sur autre methode.
      double _num;
      pmethod _pm;
      
      pmethod method(void) const { return _pm; }
      double num(void) const { return _num; }
      int type(void) const { return _type; }
      parameter(double e) : _type(0), _num(e) {}
      parameter(pmethod p) : _type(1), _pm(p) {}
    };
    
    
    typedef std::deque<parameter> param_list;
    typedef pmethod (* pfunction)(param_list &);
    typedef pmethod (* pgenfunction)(std::string);
    typedef size_t size_type;
    
  protected :
    
    std::string prefix;
    dal::dynamic_tree_sorted<std::string> suffixes;
    dal::dynamic_array<pfunction> functions;
    dal::dynamic_array<pgenfunction> genfunctions;
    int nb_genfunctions;
    
    struct meth_sto {
      pmethod pm;
      std::string name;
      
      bool operator < (const meth_sto &l) const {
	if (pm < l.pm) return true; return false;
      }
      meth_sto(void) {}
      meth_sto(pmethod p, std::string n) : pm(p), name(n) {} 
    };
    
    dal::dynamic_tree_sorted<meth_sto> meth_tab;
    
    struct meth_sto_bn {
      pmethod pm;
      std::string name;
      
      bool operator < (const meth_sto_bn &l) const {
	if (name < l.name) return true; return false;
      }
      meth_sto_bn(void) {}
      meth_sto_bn(pmethod p, std::string n) : pm(p), name(n) {} 
    };
    
    dal::dynamic_tree_sorted<meth_sto_bn> meth_tab_bn;
    
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
    typename dal::dynamic_tree_sorted<meth_sto>::const_sorted_iterator
      it = meth_tab.sorted_ge(meth_sto(pm, ""));
    const std::string *p = 0;
    if (it.index() == size_type(-1))
      DAL_THROW(dal::failure_error, "Unknown method");
    while (it.index() != size_type(-1) && (*it).pm == pm)
      { p = &((*it).name); ++it; }
    return *p;
  }
  
  template <class METHOD> std::string
  naming_system<METHOD>::shorter_name_of_method(typename
			        naming_system<METHOD>::pmethod pm)  const {
    typename dal::dynamic_tree_sorted<meth_sto>::const_sorted_iterator
      it = meth_tab.sorted_ge(meth_sto(pm, "")), it2 = it;
    const std::string *p = 0;
    size_type s = size_type(-1);
    if (it.index() == size_type(-1))
      DAL_THROW(dal::failure_error, "Unknown method");
    while (it.index() != size_type(-1) && (*it).pm == pm) {
      if (((*it).name).size() < s) {
	s = ((*it).name).size();
	p = &((*it).name); 
      }
      ++it;
    }
    
    // La boucle qui suit est à supprimer normalement
    while (it2.index() != size_type(-1) && (*it2).pm == pm) {
      if (((*it2).name).size() < s) {
	s = ((*it2).name).size();
	p = &((*it2).name);
	DAL_THROW(dal::internal_error, "This loop is not to be suppressed !!");
      }
      --it2;
    }
    
    return *p;
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
    DAL_THROW(dal::failure_error, "Invalid character on position " << i
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
	  //cerr << "params.size()=" << params.size() << ", i=" << i
	  // << ", name=" << name << ", pm=" << pm << "state="<< state << endl;
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
	DAL_THROW(dal::failure_error,
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
	std::string nname = norm_name.str();
	size_type j = meth_tab_bn.search(meth_sto_bn(0, nname));
	if (j == size_type(-1)) {
	  pm = 0;
	  for (int k = 0; k < nb_genfunctions && pm == 0; ++k) {
	    pm = (*(genfunctions[k]))(nname);
	  }
	  if (!pm) {
	    if (ind_suff == size_type(-1))
	      DAL_THROW(dal::failure_error, "Unknown method : " << nname);
	    pm = (*(functions[ind_suff]))(params);
	  }
	  j = meth_tab_bn.add(meth_sto_bn(pm, nname));
	  meth_tab.add(meth_sto(pm, nname));
	}
	return meth_tab_bn[j].pm;
      }
    }
    
  }
  
}
