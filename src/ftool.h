// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : File and string TOOL (ftool)
// File    : ftool.h : 
//           
// Date    : March 09, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
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

/**@file ftool.h
   @brief "File Tools"
*/
#ifndef FTOOL_H
#define FTOOL_H

#include <stdio.h>
#include <dal_basic.h>

namespace ftool
{
  /* ********************************************************************* */
  /*       Read a char string.                                             */
  /* ********************************************************************* */


  bool read_until(std::istream &ist, const char *st);
  bool get_token(std::istream &ist, char *st, int nb);

  struct skip {
    const char *s;
    skip(const char *s_) : s(s_) {}
  };
  std::istream& operator>>(std::istream& is, const skip& t);

  /* ********************************************************************* */
  /*       Case-insensitive string operations                              */
  /* ********************************************************************* */
  int casecmp(const char *a, const char *b, unsigned n=unsigned(-1));
  inline int casecmp(const std::string& a, const char *b, unsigned n=unsigned(-1)) {
    return casecmp(a.c_str(),b,n);
  }
  inline int casecmp(const std::string& a, const std::string& b, unsigned n=unsigned(-1)) {
    return casecmp(a.c_str(), b.c_str(),n);
  }
  inline int casecmp(char a, char b) { 
    return toupper(a) < toupper(b) ? -1 : (toupper(a) == toupper(b) ? 0 : +1); 
  }

  /* ********************************************************************* */
  /*       Clock functions.                                                */
  /* ********************************************************************* */

  /** get CPU time, in seconds. */
# ifdef HAVE_SYS_TIMES
  inline double uclock_sec(void) {
    static double ttclk = 0.;
    if (ttclk == 0.) ttclk = sysconf(_SC_CLK_TCK);
    tms t; times(&t); return double(t.tms_utime) / ttclk;
  }
# else
  inline double uclock_sec(void)
  { return double(clock())/double(CLOCKS_PER_SEC); }
# endif

  /* ********************************************************************* */
  /*       Read a parameter file.                                          */
  /* ********************************************************************* */
 
  class md_param {
    public :
      struct t_value {
	int type;  /* 0: nothing, 1 : real, 2 : integer, 3 : string,       */
	/* 4 : array.                                                      */
	union {
	  const char *v_string;
	  double v_real;
	  long v_int;
	  struct {
	    dal::dynamic_array<t_value> *list;
	    int nb;
	  } v_list;
	} value;
	
	t_value(void) { type = 0; }

	void  clear_value(void) {
	  switch(type) {
	  case 3:
	    delete[] value.v_string;
	    break;
	  case 4: delete value.v_list.list;
	  }
	  type = 0;
	}
	~t_value(void)
	  { clear_value(); }
      };

    struct t_param {
      const char *name;   /* parameter name.                               */
      const char *comment;
      t_value value;
      
      void param(void) { value.type = 0; name = comment = 0; }
      
      double &real_value(void) { return value.value.v_real; }
      long &int_value(void) { return value.value.v_int; }
      const char * &string_value(void) { return value.value.v_string; }
      int &type_value(void) { return value.type; }
      int &nb_sub_param(void) { return value.value.v_list.nb; }
      dal::dynamic_array<t_value> * &sub_param(void)
	{ return value.value.v_list.list; }
      
      void clear_value(void)
	{ value.clear_value(); }
      
      ~t_param(void) {
	if (name != 0) delete[] name;
	if (comment != 0)  delete[] comment;
      }
      t_param() : name(0), comment(0) {}
    };
    
    protected :
      dal::dynamic_array<t_param> param_list;
    int nb_param;
    
    int search_param(const char *name);
    
    int state, nbcharread;
    char string_read[255];
    char name_read[255];
    int lex_read_char(char c);
    int read_char(char c);
    
    int blk1, blk2, blk3, blk4, dts;
    bool is_text;
    FILE *fid;
    int f_status; /* 1 = read, 2 = write.                                  */
    int lblk_count;
    char *temp_blk;
    int temp_blk_size;
    int clk; int flushtime;
    
    public :
      
      md_param(void) { 
      nb_param = 0; state = 0; nbcharread = 0; blk1 = blk2 = blk3 = blk4 = 1;
      dts = 4; is_text = true; f_status = 0; temp_blk_size = 0; flushtime = 20;
      }
    /***********************************************************************/
    /* Table of parameters management.                                     */
    /***********************************************************************/

    const char * &param_name(int i) { return param_list[i].name; }
    double &real_value(int i) { return param_list[i].real_value(); }
    long &int_value(int i) { return param_list[i].int_value(); }
    const char * &string_value(int i) { return param_list[i].string_value(); }
    int &param_type(int i) { return param_list[i].type_value(); }
    void clear_value(int i) { param_list[i].clear_value(); }
    const char * &param_comment(int i) { return param_list[i].comment; }
    int &nb_sub_param(int i) { return param_list[i].nb_sub_param(); }
    dal::dynamic_array<t_value> * &sub_param(int i)
      { return param_list[i].sub_param();}
    double real_value(const char *name);
    double real_value(const char *name, const char *comment);
    long int_value(const char *name);
    long int_value(const char *name, const char *comment);
    const char *string_value(const char *name);
    const char *string_value(const char *name, const char *comment);
    int nb_sub_param(const char *name);
    int nb_sub_param(const char *name, const char *comment);
    long sub_int_value(const char *name, int n);
    const char *sub_string_value(const char *name, int n);
    double sub_real_value(const char *name, int n);
    int add_real_param(const char *name, double e);
    int add_int_param(const char *name, long e);
    int add_string_param(const char *name, const char *st);
    int add_string_param_to_param(const char *name, const char *value);
    
    /********************************************************************/
    /* Reading and writting parameters.                                 */
    /********************************************************************/

    int read_string(const char *st);
  
    void read_param_file(const char *fn); /* read parameters in a file. */
    void read_command_line(int argc, char *argv[]);

    void write_param_file(const char *fn); /*write parameters in a file.*/
    void write_param_file(void);

    /********************************************************************/
    /* Data file management.                                            */
    /********************************************************************/
    
    void blk_size(int a, int b, int c, int d)
      { blk1 = a;  blk2 = b; blk3 = c; blk4 = d; }
    void set_to_text(void) { add_string_param("DATA_TYPE", "TEXT");}
    void set_to_bin(void) { add_string_param("DATA_TYPE", "BIN"); }
    void ropen_data_file(void);
    void wopen_data_file(void);
    void close_data_file(void);
    void write_in_data_file(int nb, double *array);
    void load_from_data_file(int nb, double *array);
  };

}

#endif
