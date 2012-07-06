/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2000-2012 Yves Renard
 
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

/**@file bgeot_ftool.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date March 09, 2000.
   @brief "File Tools"
*/

#ifndef BGEOT_FTOOL_H
#define BGEOT_FTOOL_H

#include <stdio.h>
#include "dal_basic.h"
#include <map>

namespace bgeot
{
  /* ********************************************************************* */
  /*       Read a char string.                                             */
  /* ********************************************************************* */

  bool read_until(std::istream &ist, const char *st);

  /** Very simple lexical analysis of general interest for reading small
   *  langages with a "MATLAB like" syntax :
   *  spaces are ignored, '%' indicates a commentary until the end of the line,
   *  '...' indicates that the instruction continue on the next line
   *  (or separate two sub part of the same character string).
   *
   *  The function returns
   *    0 if there is nothing else to read in the file
   *    1 if an end line has been found (st is empty in this case)
   *    2 if a number as been read
   *    3 if a string has been read
   *    4 if a alphanumeric name has been read
   *    5 for a one character symbol
   *    6 for a two characters symbol (<=, >=, ==, !=, &&, ||)
   *
   *  Note that when the end of line is not significant the option ignore_cr
   *  allows to consider the carriage return as a space character.
   */
  int get_token(std::istream &ist, std::string &st,
		bool ignore_cr = true, bool to_up = true, 
		bool read_un_pm = true, int *linenb = 0);

  struct skip {
    const char *s;
    skip(const char *s_) : s(s_) {}
  };
  std::istream& operator>>(std::istream& is, const skip& t);

  /* ********************************************************************* */
  /*       Case-insensitive string operations                              */
  /* ********************************************************************* */
  int casecmp(const char *a, const char *b, unsigned n=unsigned(-1));

  inline int casecmp(const std::string& a, const char *b,
		     unsigned n=unsigned(-1))
  { return casecmp(a.c_str(),b,n); }

  inline int casecmp(const std::string& a, const std::string& b,
		     unsigned n=unsigned(-1))
  { return casecmp(a.c_str(), b.c_str(),n); }

  inline int casecmp(char a, char b)
  { return toupper(a)<toupper(b) ? -1 : (toupper(a) == toupper(b) ? 0 : +1); }

  /* ********************************************************************* */
  /*       Read a parameter file.                                          */
  /* ********************************************************************* */
 
  // The associated langage has approximatively the following grammar:
  //
  // 'instruction' := 'parameter_name' '=' 'expression';
  //              or 'if' 'expression'
  //                 'instruction_list'
  //                 [ 'else' 'instruction_list' ]
  //                 'end'          
  // 'expression' := '[' 'expression' ',' ... ']'
  //              or 'parameter_name'
  //              or 'numeric_value'
  //              or 'character_string'
  //              or '(' 'expression' ')'
  //              or '+' 'expression'
  //              or '-' 'expression'
  //              or '!' 'expression'
  //              or 'expression' '+' 'expression'
  //              or 'expression' '-' 'expression'
  //              or 'expression' '*' 'expression'
  //              or 'expression' '/' 'expression'
  //              or 'expression' '||' 'expression'
  //              or 'expression' '&&' 'expression'
  //              or 'expression' '==' 'expression'
  //              or 'expression' '!=' 'expression'
  //              or 'expression' '<' 'expression'
  //              or 'expression' '>' 'expression'
  //              or 'expression' '<=' 'expression'
  //              or 'expression' '>=' 'expression'

  class md_param {
  public :

    typedef enum { REAL_VALUE, STRING_VALUE, ARRAY_VALUE } param_type;

    class param_value {
      param_type pt;
      double real_value;
      std::string string_value;
      std::vector<param_value> array_value;
      
    public :
      param_type type_of_param(void) const { return pt; }
      double &real(void) { return real_value; }
      double real(void) const { return real_value; }
      std::string &string(void) { return string_value; }
      const std::string &string(void) const { return string_value; }
      std::vector<param_value> &array(void) { return array_value; }
      const std::vector<param_value> &array(void) const { return array_value; }
      param_value(double e = 0.0) : pt(REAL_VALUE), real_value(e) {}
      param_value(std::string s) : pt(STRING_VALUE), real_value(0.0),
				   string_value(s) {}
      param_value(char *s) :  pt(STRING_VALUE), real_value(0.0),
			      string_value(s) {}
      param_value(param_type p): pt(p), real_value(0.0) {}
    };
    
  protected :

    std::map<std::string, param_value> parameters;
    bool token_is_valid;
    int current_line;
    std::string current_file;

    int get_next_token(std::istream &f);
    void valid_token(void);
    std::string temp_string;
    param_value read_expression_list(std::istream &f, bool skipped);
    param_value read_expression(std::istream &f, bool skipped);
    int read_instruction_list(std::istream &f, bool sk = false);
    int read_instruction(std::istream &f, bool sk = false);
    void parse_error(const std::string &t);
    void syntax_error(const std::string &t);
    void do_bin_op(std::vector<md_param::param_value> &value_list,
		   std::vector<int> &op_list, std::vector<int> &prior_list);
    
  public :
      

    double real_value(const std::string &name, const char *comment = 0);
    long int_value(const std::string &name, const char *comment = 0);
    const std::string &string_value(const std::string &name,
				    const char *comment = 0);
    const std::vector<param_value> &array_value(const std::string &name,
						const char *comment = 0);
    void add_int_param(const std::string &name, long e)
    { parameters[name] = param_value(double(e)); }
    void add_real_param(const std::string &name, double e)
    { parameters[name] = param_value(e); }
    void add_string_param(const std::string &name, const std::string &s)
    { parameters[name] = param_value(s); }
 
    // todo : add the access to the arrays
    
    void read_param_file(std::istream &f);

    /** Read the parameters on the command line. If a name is found the 
     *	corresponding .param file is searched. If a -dNOM=VALUE is found
     *	(or -d NOM=VALUE), it is evaluated.
     */
    void read_command_line(int argc, char *argv[]);
  };

}


namespace ftool {
 
  // For compatibility with Getfem 2.0

  using bgeot::md_param;

}

#endif
