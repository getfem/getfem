/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_except.h : Exceptions.                                   */
/*     									   */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, Julien.pommier@gmm.insa-tlse.fr                */
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

#ifndef __DAL_EXCEPT_H
#define __DAL_EXCEPT_H

#include <dal_std.h>

/* *********************************************************************** */
/*	Getfem++ generic errors.                     			   */
/* *********************************************************************** */

  /* errors definitions  */

  class dimension_error : public std::logic_error {
  public:
    dimension_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class file_not_found_error : public std::logic_error {
  public:
    file_not_found_error(const std::string& what_arg)
      : std::logic_error (what_arg) { }
  };

  class internal_error : public std::logic_error {
  public:
    internal_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class failure_error : public std::logic_error {
  public:
    failure_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class not_linear_error : public std::logic_error {
  public:
    not_linear_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class to_be_done_error : public std::logic_error {
  public:
    to_be_done_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  #define DAL_STANDARD_CATCH_ERROR   catch(std::logic_error e) \
  { \
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::runtime_error e)\
  {\
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::bad_alloc) {\
    cerr << "============================================\n";\
    cerr << "|  A bad allocation has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  }\
  catch(std::bad_typeid) { \
    cerr << "============================================\n";\
    cerr << "|  A bad typeid     has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(std::bad_exception) { \
    cerr << "============================================\n";\
    cerr << "|  A bad exception  has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(std::bad_cast) { \
    cerr << "============================================\n";\
    cerr << "|    A bad cast  has been detected !!!     |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(...) {\
    cerr << "============================================\n";\
    cerr << "|  An unknown error has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  }
//   catch(ios_base::failure) { 
//     cerr << "============================================\n";
//     cerr << "| A ios_base::failure has been detected !!!|\n";
//     cerr << "============================================\n";
//     exit(1);
//   } 

  class exception_callback {
  public:
    virtual void callback(const std::string& msg) {};
  };
  
  void do_exception_callback(const std::string &msg);
  void set_exception_callback(exception_callback *e);

#define DAL_THROW(type, thestr) {                \
    std::stringstream msg;                       \
    msg << "in "__FILE__ << ", line "            \
        << __LINE__ << ": \n" << thestr << ends; \
    dal::do_exception_callback(msg.str());       \
    throw (type)(msg.str());                     \
  }
  
} /* end of namespace dal.                                                */



#endif /* __DAL_EXCEPT_H */
