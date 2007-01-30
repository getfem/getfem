// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/** @file gmm_except.h 
    @author  Yves Renard <Yves.Renard@insa-toulouse.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
    @date September 01, 2002.
    @brief Definition of basic exceptions.
*/

#ifndef GMM_EXCEPT_H__
#define GMM_EXCEPT_H__

#include "gmm_std.h"


namespace gmm {

/* *********************************************************************** */
/*	Getfem++ generic errors.                     			   */
/* *********************************************************************** */

  /* errors definitions  */

  using std::invalid_argument;

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

  #define GMM_STANDARD_CATCH_ERROR   catch(std::logic_error e) \
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

  /* callback handler for gmm/getfem exceptions */
  struct exception_callback {
    virtual void callback(const std::string&) = 0; //{};

    static exception_callback *which_except(exception_callback *p = 0) {
      static exception_callback *exc_cback = 0;
      if (p != 0) exc_cback = p;
      return exc_cback;
    }

    static void do_exception_callback(const std::string &msg)
      { if (which_except()) which_except()->callback(msg); }
  
    static void set_exception_callback(exception_callback *e)
      { which_except(e); }
    virtual ~exception_callback() {}
  };

  /* crashing callback for debug mode */
  struct exception_callback_debug : public gmm::exception_callback  {
    virtual void callback(const std::string& msg)
    { cerr << msg << endl; *(int *)(0) = 0; }
    virtual ~exception_callback_debug() {}
  };

#define GMM_SET_EXCEPTION_DEBUG {                                                       \
    gmm::exception_callback::set_exception_callback(new gmm::exception_callback_debug);	\
  }

  /** user function for changing the default exception callback */ 
  inline void set_exception_callback(exception_callback *e)
  { exception_callback::which_except(e); }

#ifdef GETFEM_HAVE_PRETTY_FUNCTION
#  define GMM_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else 
#  define GMM_PRETTY_FUNCTION ""
#endif

#define GMM_THROW(type, thestr) {					\
    std::stringstream msg;						\
    msg << "Error in "__FILE__ << ", line "				\
        << __LINE__ << " " << GMM_PRETTY_FUNCTION << ": \n" << thestr << ends; \
    gmm::exception_callback::do_exception_callback(msg.str());		\
    throw (type)(msg.str());						\
  }

#ifdef DEBUG_MODE
#  define GMM_INTERNAL_ERROR(thestr) { \
  cerr << "Internal error: " << GMM_PRETTY_FUNCTION << " " << thestr << endl; \
   ::abort(); \
   }
#else
#  define GMM_INTERNAL_ERROR(thestr) GMM_THROW(gmm::internal_error, "Internal error: " << thestr)
#endif

/* *********************************************************************** */
/*	Getfem++ warnings.                         			   */
/* *********************************************************************** */

  // This allows to dynamically hide warnings
  struct warning_level {
    static int level(int l = -2)
    { static int level_ = 3; return (l != -2) ? (level_ = l) : level_; }
  };

  inline void set_warning_level(int l) { warning_level::level(std::max(0,l)); }

  // This allow not too compile some Warnings
#ifndef GMM_WARNING_LEVEL
# define GMM_WARNING_LEVEL 4
#endif

  // Warning levels : 0 always printed
  //                  1 very important : specify a possible error in the code.
  //                  2 important : specify a default of optimization for inst.
  //                  3 remark
  //                  4 ignored by default.

#define GMM_WARNING_MSG(level_, thestr)  {			       \
      std::stringstream msg;                                           \
      msg << "Level " << level_ << " Warning in "__FILE__ << ", line " \
          << __LINE__ << ": " << thestr << ends;		       \
       std::cerr << msg.str() << std::endl;                            \
    }

#define GMM_WARNING0(thestr) GMM_WARNING_MSG(0, thestr)

#if GMM_WARNING_LEVEL > 0
# define GMM_WARNING1(thestr)                                           \
  { if (1 <= gmm::warning_level::level()) GMM_WARNING_MSG(1, thestr) }
#else
# define GMM_WARNING1(thestr) {}
#endif

#if GMM_WARNING_LEVEL > 1
# define GMM_WARNING2(thestr)                                           \
  { if (2 <= gmm::warning_level::level()) GMM_WARNING_MSG(2, thestr) } 
#else
# define GMM_WARNING1(thestr) {}
#endif

#if GMM_WARNING_LEVEL > 2
# define GMM_WARNING3(thestr)                                           \
  { if (3 <= gmm::warning_level::level()) GMM_WARNING_MSG(3, thestr) } 
#else
# define GMM_WARNING1(thestr) {}
#endif

#if GMM_WARNING_LEVEL > 3
# define GMM_WARNING4(thestr)                                           \
  { if (4 <= gmm::warning_level::level()) GMM_WARNING_MSG(4, thestr) } 
#else
# define GMM_WARNING1(thestr) {}
#endif

/* *********************************************************************** */
/*	Getfem++ traces.                         			   */
/* *********************************************************************** */

  // This allows to dynamically hide traces
  struct traces_level {
    static int level(int l = -2)
    { static int level_ = 3; return (l != -2) ? (level_ = l) : level_; }
  };

  inline void set_traces_level(int l) { traces_level::level(std::max(0,l)); }

  // This allow not too compile some Warnings
#ifndef GMM_TRACES_LEVEL
# define GMM_TRACES_LEVEL 3
#endif

  // Traces levels : 0 always printed
  //                 1 Susceptible to occur once in a program.
  //                 2 Susceptible to occur occasionnaly in a program (10).
  //                 3 Susceptible to occur often (100).
  //                 4 Susceptible to occur very often (>1000).

#define GMM_TRACE_MSG_MPI     // for Parallelized version
#define GMM_TRACE_MSG(level_, thestr)  {			       \
   GMM_TRACE_MSG_MPI {                                                  \
      std::stringstream msg;                                           \
      msg << "Trace " << level_ << " in "__FILE__ << ", line "         \
          << __LINE__ << ": " << thestr  \
          << ends;						       \
       std::cout << msg.str() << std::endl;                            \
    }                                                                  \
  }        

#define GMM_TRACE0(thestr) GMM_TRACE_MSG(0, thestr)

#if GMM_TRACES_LEVEL > 0
# define GMM_TRACE1(thestr)                                           \
  { if (1 <= gmm::traces_level::level()) GMM_TRACE_MSG(1, thestr) }
#else
# define GMM_TRACE1(thestr) {}
#endif

#if GMM_TRACES_LEVEL > 1
# define GMM_TRACE2(thestr)                                           \
  { if (2 <= gmm::traces_level::level()) GMM_TRACE_MSG(2, thestr) } 
#else
# define GMM_TRACE2(thestr) {}
#endif

#if GMM_TRACES_LEVEL > 2
# define GMM_TRACE3(thestr)                                           \
  { if (3 <= gmm::traces_level::level()) GMM_TRACE_MSG(3, thestr) } 
#else
# define GMM_TRACE3(thestr) {}
#endif

#if GMM_TRACES_LEVEL > 3
# define GMM_TRACE4(thestr)                                           \
  { if (4 <= gmm::traces_level::level()) GMM_TRACE_MSG(4, thestr) } 
#else
# define GMM_TRACE4(thestr) {}
#endif

 
}


#endif /* GMM_EXCEPT_H__ */
