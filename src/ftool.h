/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  File and string TOOL (ftool)                                 */
/* File    :  ftool.h :                                                    */
/*     									   */
/*                                                                         */
/* Date : March 09, 2000.                                                  */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
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

#include <stdio.h>
#include <dal_basic.h>

namespace ftool
{
  /* ********************************************************************* */
  /*       Read a char string.                                             */
  /* ********************************************************************* */


  bool read_untill(std::istream &ist, const char *st);
  bool get_token(std::istream &ist, char *st, int nb);

  /* ********************************************************************* */
  /*       Clock functions.                                                */
  /* ********************************************************************* */

  double uclock_sec(void);

  /* ********************************************************************* */
  /*       Read a parameter file.                                          */
  /* ********************************************************************* */
 
  class md_param {
    public :
      struct t_value {
	int type;  /* 0: vide, 1 : reel, 2 : entier, 3 : chaine,           */
	/* 4 : tableau.                                                    */
	union {
	  char *v_string;
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
	  case 3: // cout << "suppression de " << value.v_string << endl;
	    delete[] value.v_string; // cout << "faite ... /n";
	    break;
	  case 4: delete value.v_list.list;
	  }
	  type = 0;
	}
	~t_value(void)
	  { clear_value(); }
      };

    struct t_param {
      char *name;   /* Nom du parametre.                                  */
      char *comment;
      t_value value;
      
      void param(void) { value.type = 0; name = comment = NULL; }
      
      double &real_value(void) { return value.value.v_real; }
      long &int_value(void) { return value.value.v_int; }
      char * &string_value(void) { return value.value.v_string; }
      int &type_value(void) { return value.type; }
      int &nb_sub_param(void) { return value.value.v_list.nb; }
      dal::dynamic_array<t_value> * &sub_param(void)
	{ return value.value.v_list.list; }
      
      void clear_value(void)
	{ value.clear_value(); }
      
      ~t_param(void) {
	if (name != NULL) delete[] name;
	if (comment != NULL)  delete[] comment;
      }
    };
    
    protected :
      dal::dynamic_array<t_param> param_list;
    int nb_param;
    
    int search_param(char *name);
    
    int state, nbcharread;
    char string_read[255];
    char name_read[255];
    int lex_read_char(char c);
    int read_char(char c);
    
    int blk1, blk2, blk3, blk4, dts;
    bool is_text;
    FILE *fid;
    int f_status; /* 1 = ouvert en lecture, 2 en ecriture.              */
    int lblk_count;
    char *temp_blk;
    int temp_blk_size;
    int clk; int flushtime;
    
    public :
      
      md_param(void) { 
      nb_param = 0; state = 0; nbcharread = 0; blk1 = blk2 = blk3 = blk4 = 1;
      dts = 4; is_text = true; f_status = 0; temp_blk_size = 0; flushtime = 20;
      }
    /********************************************************************/
    /* Table of parameters management.                                  */
    /********************************************************************/

    char * &param_name(int i) { return param_list[i].name; }
    double &real_value(int i) { return param_list[i].real_value(); }
    long &int_value(int i) { return param_list[i].int_value(); }
    char * &string_value(int i) { return param_list[i].string_value(); }
    int &param_type(int i) { return param_list[i].type_value(); }
    void clear_value(int i) { param_list[i].clear_value(); }
    char * &param_comment(int i) { return param_list[i].comment; }
    int &nb_sub_param(int i) { return param_list[i].nb_sub_param(); }
    dal::dynamic_array<t_value> * &sub_param(int i)
      { return param_list[i].sub_param();}
    double real_value(char *name); /* regarde si le param existe,       */
    /* retourne la valeur ou 0.0.                                       */
    double real_value(char *name, char *comment); /* idem mais demande  */
    /* a l'utilisateur en cas de non existence du parametre.            */
    long int_value(char *name);
    long int_value(char *name, char *comment);
    char *string_value(char *name);
    char *string_value(char *name, char *comment);
    int nb_sub_param(char *name);
    int nb_sub_param(char *name, char *comment);
    long sub_int_value(char *name, int n);
    char *sub_string_value(char *name, int n);
    double sub_real_value(char *name, int n);
    int add_real_param(char *name, double e);
    int add_int_param(char *name, long e);
    int add_string_param(char *name, char *st);
    int add_string_param_to_param(char *name, char *value);
    
    /********************************************************************/
    /* Reading and writting parameters.                                 */
    /********************************************************************/

    int read_string(char *st); /* lit des parametres dans une chaine    */
           /* de caractere.                                             */
           /* rend 0 si tout c'est bien passe.                          */
  
    void read_param_file(char *fn); /* read parameters in a file.       */
    void read_command_line(int argc, char *argv[]);

    void write_param_file(char *fn); /* write all parameters in a file. */
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
