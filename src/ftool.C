/* *********************************************************************** */
/*                                                                         */
/* Library :  File and string TOOL (ftool)                                 */
/* File    :  ftool.C :                                                    */
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


#include <ftool.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>

namespace ftool
{

#ifndef CLK_TCK
#define TTCLK (double(CLOCKS_PER_SEC) / 10000.0)
#else
#define TTCLK double(CLK_TCK)
#endif

  double uclock_sec(void)
  { tms t; times(&t); return double(t.tms_utime) / TTCLK; }

  bool read_untill(STD_NEEDED istream &ist, const char *st)
  {
    int i = 0, l = strlen(st); char c;
    while (!ist.eof() && i < l)
      { ist.get(c); if (toupper(c) == st[i]) i++; else i = 0; }
    if (ist.eof()) return false; else return true;
  }

  bool get_token(STD_NEEDED istream &ist, char *st, int nb)
  {
    char c;
    int i = 0;
    bool te = false;
    st[0] = 0;
    if (ist.eof()) return false;
    ist.get(c);
    while (!te)
    {
      while (isspace(c)) { if (ist.eof()) return true; ist.get(c); }
      if (c == '%')
	{ while (c != '\n') { if (ist.eof()) return true; ist.get(c); } }
      else { te = true; }
    }
    while (i < nb-1 && !isspace(c) && !iscntrl(c))
      { st[i++] = toupper(c); if (ist.eof()) break; ist.get(c); }
    st[i] = 0;
    return true;
  }


  static char temp_string[512];
  
  int md_param::search_param(char *name) {
    for (int i = 0; i < nb_param; i++)
      if (!strcmp(name, param_name(i))) return i;
    return -1;
  }


  /*************************************************************************/
  /* type des parametres :                                                 */
  /*     0 : inconnu.                                                      */
  /*     1 : reel (ou entier).                                             */
  /*     2 : entier.                                                       */
  /*     3 : chaine de caractere.                                          */
  /*     4 : liste.                                                        */
  /*************************************************************************/


  /*************************************************************************/
  /* L'analyse syntaxique est la suivante.                                 */
  /* classe de caracteres :                                                */
  /*     0 : espace ou ';' ou TAB.                                         */
  /*     1 : % indiquant debut de commentaire.                             */
  /*     2 : " ou '                                                        */
  /*     3 : =                                                             */
  /*     4 : Alphabetique + '_'                                            */
  /*     5 : Numerique + '.' + '+' + '-'                                   */
  /*     6 : fin de ligne.                                                 */ 
  /*     7 : Autre.                                                        */
  /*     8 : '['                                                           */
  /*     9 : ']'                                                           */
  /* Les transitions sont :                                                */
  /* ---------------------------------------------------------------       */
  /* |car\etat | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12|       */
  /* ---------------------------------------------------------------       */
  /* |  0      | 0 | 2 | 2 | 3 | 0 | 5 | 0 | 7 | 8 | 9 | 10| 8 | 8 |       */
  /* |  1      | 7 | -1| -1| -1| 7 | 5 | 7 | 7 | 9 | 9 | 10| 9 | 9 |       */
  /* |  2      | -1| -1| -1| 5 | -1| 6 | 5 | 7 | 10| 9 | 12| -1| 10|       */
  /* |  3      | -1| 3 | 3 | -1| -1| 5 | -1| 7 | -1| 9 | 10| -1| -1|       */
  /* |  4      | 1 | 1 | -1| -1| 4 | 5 | 1 | 7 | -1| 9 | 10| 11| -1|       */
  /* |  5      | -1| 1 | -1| 4 | 4 | 5 | -1| 7 | 11| 9 | 10| 11| -1|       */
  /* |  6      | 0 | 2 | 2 | 3 | 0 | -1| 0 | 0 | 8 | 8 | -1| 8 | 8 |       */
  /* |  7      | -1| -1| -1| -1| -1| 5 | -1| 7 | -1| 9 | 10| -1| -1|       */
  /* |  8      | -1| -1| -1| 8 | -1| 5 | -1| 7 | -1| 9 | 10| -1| -1|       */
  /* |  9      | -1| -1| -1| -1| -1| 5 | -1| 7 | 0 | 9 | 10| 0 | 0 |       */
  /* ---------------------------------------------------------------       */
  /*                                                                       */
  /* ou les etats sont :                                                   */
  /* 0 : attente de la lecture d'un parametre.                             */
  /* 1 : lecture du nom du parametre.                                      */ 
  /* 2 : attente de la lecture du =                                        */
  /* 3 : lecture du =                                                      */
  /* 4 : lecture de la valeur numerique.                                   */
  /* 5 : lecture d'une chaine de caractere.                                */
  /* 6 : lecture du " fermant d'une chaine de caractere ou double ".       */
  /* 7 : lecture d'un commentaire.                                         */
  /* 8 : debut lecture d'un tableau, attente premiere valeur.              */
  /* 9 : lecture d'un commentaire dans un tableau.                         */
  /* 10: lecture d'une chaine de caractere dans un tableau.                */
  /* 11: lecture d'une valeur numerique dans un tableau.                   */
  /* 12: lecture du " fermant d'une chaine dans un tableau.                */
  /*************************************************************************/
  
  
  static int __automat[10][13]=
  {  0,  2,  2,  3,  0,  5,  0,  7,  8,  9,  10,  8,  8,
     7, -1, -1, -1,  7,  5,  7,  7,  9,  9,  10,  9,  9, 
     -1, -1, -1,  5, -1,  6,  5,  7, 10,  9,  12, -1, 10,
     -1,  3,  3, -1, -1,  5, -1,  7, -1,  9,  10, -1, -1,
     1,  1, -1, -1,  4,  5,  1,  7, -1,  9,  10, 11, -1,
     -1,  1, -1,  4,  4,  5, -1,  7, 11,  9,  10, 11, -1,
     0,  2,  2,  3,  0, -1,  0,  0,  8,  8,  -1,  8,  8,
     -1, -1, -1, -1, -1,  5, -1,  7, -1,  9,  10, -1, -1, 
     -1, -1, -1,  8, -1,  5, -1,  7, -1,  9,  10, -1, -1,
     -1, -1, -1, -1, -1,  5, -1,  7,  0,  9,  10,  0,  0
  };
  
  static int __car_type_amp(char c) {
    switch (c) {
    case ' ' :  case ';' : case 9 : return 0;
    case 0 : case 13 : case 10 : return 6;
    case '%' : return 1;
    case '"' : case '\'' : return 2;
    case '=' : return 3;
    case '_' : return 4;
    case '+' : case '-' : case '.' : return 5;
    case '[' : return 8;
    case ']' : return 9;
    }
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) return 4;
    if ((c >= '0' && c <= '9')) return 5;
    return 7;
  }
  
  
  int md_param::read_char(char c) {
    int crt = __car_type_amp(c);
    int newstate = __automat[crt][state];
    // STD_NEEDED cout << "car : " << c << " type : " << crt << " old state : "
    //     << state << " new state : " << newstate << endl; getchar();
    if (newstate == 1 || newstate == 4 || newstate == 5
	|| newstate == 10 || newstate == 11 ) { 
      if (nbcharread >= 254) STD_NEEDED cerr << "String too long\n";
      else string_read[nbcharread++] = c;
    }
    
    if (state != newstate) {
      switch (newstate) {
      case 2 : case 3 :
	if (state == 1) {
	  string_read[nbcharread] = 0;
	  strcpy(name_read, string_read);
	  nbcharread = 0;
	}
	break;
      case 7 : case 0 : case 1 : case 8 : case 9 :
	if (state == 6 || state == 4) {
	  string_read[nbcharread] = 0;
	  add_string_param(name_read, string_read);
	  nbcharread = 0;
	}
	if (state == 12 || state == 11) {
	  add_string_param_to_param(name_read, string_read);
	  nbcharread = 0;
	}
	break;
      case 5 : nbcharread = 0; break;
      }
    }
    
    state = newstate; if (newstate == -1) return -1; else return 0;
  }
  
  int md_param::read_string(char *st) {
    for (int i = 0; i < strlen(st); i++) {
      if (read_char(st[i]) == -1) return -1;
    }
    if (read_char(13) == -1) return -1;
    return 0;
  }
  
  void md_param::read_param_file(char *fn) {
    FILE *F = fopen(fn, "r");
    if (F == NULL) { STD_NEEDED cout << "MODEL_PARAM : file not found.\n"; exit(1); }
    for(;;) {
      char c = getc(F); if (feof(F)) break;
      if (read_char(c) == -1)
	{ STD_NEEDED cout << "MODEL_PARAM : syntax error in file " << fn << endl; exit(1);}
    }
    if (read_char(13) == -1) { STD_NEEDED cout << "MODEL_PARAM : syntax error"; exit(1); }
    if (state != 0) { STD_NEEDED cout << "MODEL_PARAM : incorrect end of file"; exit(1); }
    fclose(F);
  }
  
  void md_param::read_command_line(int argc, char *argv[]) /* lit les        */
  { /* parametres sur la ligne de commande. Si un nom est trouve, on cherche */
    /* le fichier .param correspondant, si une option -dNOMP=VALUE est       */
    /* trouvee, elle est evaluee. Cela laisse la place pour d'autres options */
    /* pour le programme lui meme.                                           */
    
    for (int aa = 1; aa < argc; aa++) {
      if (argv[aa][0] != '-')
	{
	  strcpy(temp_string, argv[aa]);
	  FILE *F = fopen(temp_string, "r");
	  if (F == NULL)
	    { strcat(temp_string, ".param"); F = fopen(temp_string, "r"); }
	  if (F == NULL)
	    {
	      sprintf(temp_string,"%s%s", argv[aa], ".m");
	      F = fopen(temp_string, "r");
	    }
	  if (F != NULL) { fclose(F); read_param_file(temp_string); }
	}
      else if (argv[aa][1] == 'd')
	{
	  if (strlen(argv[aa]) == 2)
	    { if (aa < argc - 1) read_string(argv[++aa]); }
	  else
	    read_string(&(argv[aa][2]));
	}
    }
  }
  
  static void __write_value_file(FILE *F, md_param::t_value &v) {
    switch(v.type) {
    case 1 : fprintf(F, "%g", v.value.v_real); break;
    case 2 : fprintf(F, "%ld", v.value.v_int); break;
    case 3 : fprintf(F, "'%s'", v.value.v_string); break;
    case 4 : fprintf(F, "[ ");
      for (int j = 0; j < v.value.v_list.nb; j++ ) {
	__write_value_file(F, (*v.value.v_list.list)[j]);
	fprintf(F, " ");
      }
      fprintf(F, "]");
      break;
    }
  }
  
  
  void md_param::write_param_file(char *fn) {
    FILE *F = fopen(fn, "w");
    add_int_param("BLKSIZE1", blk1); add_int_param("BLKSIZE2", blk2);
    add_int_param("BLKSIZE3", blk3); add_int_param("BLKSIZE4", blk4);
    if (is_text)
      add_string_param("DATA_TYPE", "TEXT");
    else
      add_string_param("DATA_TYPE", "BIN");
    if (dts == 4)
      add_string_param("DATA_FORMAT", "float32");
    else
      add_string_param("DATA_FORMAT", "float64");
    
    for (int i = 0; i < nb_param; i++) {
      fprintf(F, "%s = ", param_name(i));
      __write_value_file(F, param_list[i].value);
      
      if (param_comment(i) !=NULL)
	fprintf(F, ";\t\t%% %s\n", param_comment(i));
      else
	fprintf(F, ";\n");
    }
    fclose(F);
  }
  
  void md_param::write_param_file(void) {
    sprintf(temp_string, "%s%s",
	    string_value("ROOTFILENAME", "Nom du fichier de donnees"), ".m");
    write_param_file(temp_string);
  }
  
  
  double md_param::real_value(char *name) { 
    int i = search_param(name);
    double e;
    if (i == -1) return 0.0;
    
    switch(param_type(i)) {
    case 1 : return real_value(i);
    case 2 : return double(int_value(i));
    case 3 : e = atof(string_value(i));
      clear_value(i);
      add_real_param(name,e);
      return e;
    }
    return 0.0;
  }
  
  double md_param::real_value(char *name, char *comment) {
    int i = search_param(name);
    if (i == -1) {
      double f;
      STD_NEEDED cout << "No parameter " << name << " found, please enter its value\n";
      STD_NEEDED cout << comment << " : "; STD_NEEDED cin >> f; sprintf(string_read, "%g", f);
      i = add_real_param(name, atof(string_read));
    }
    if (param_comment(i) != NULL) delete param_comment(i);
    param_comment(i) = new char[strlen(comment)+1];
    strcpy(param_comment(i), comment);
    return real_value(name);
  }
  
  long md_param::int_value(char *name) { 
    int i = search_param(name);
    long e;
    if (i == -1) return 0;
    
    switch(param_type(i)) {
    case 1 : return long(real_value(i));
    case 2 : return int_value(i);
    case 3 : e = atol(string_value(i));
      clear_value(i);
      add_int_param(name,e);
      return e;
    }
    return 0;
  }
  
  long md_param::int_value(char *name, char *comment) {
    int i = search_param(name);
    if (i == -1) {
      long f;
      STD_NEEDED cout << "\nNo parameter " << name << " found, please enter its value\n";
      STD_NEEDED cout << comment << " : "; STD_NEEDED cin >> f; sprintf(string_read, "%ld", f);
      i = add_int_param(name, atol(string_read));
    }
    if (param_comment(i) != NULL) delete param_comment(i);
    param_comment(i) = new char[strlen(comment)+1];
    strcpy(param_comment(i), comment);
    return int_value(name);
  }
  
  char *md_param::string_value(char *name) {
    int i = search_param(name); if (i==-1) return NULL; return string_value(i);
  }
  
  char *md_param::string_value(char *name, char *comment) {
    int i = search_param(name);
    if (i == -1) {
      STD_NEEDED cout << "No parameter " << name << " found, please enter its value\n";
      STD_NEEDED cout << comment << " : "; STD_NEEDED cin >> string_read;
      i = add_string_param(name, string_read);
    }
    if (param_comment(i) != NULL) delete param_comment(i);
    param_comment(i) = new char[strlen(comment)+1];
    strcpy(param_comment(i), comment);
    return string_value(name);
  }
  
  int md_param::nb_sub_param(char *name) {
    int i = search_param(name); if (i == -1) return 0;
    if (param_type(i) != 4) return 0;
    return nb_sub_param(i);
  }
  
  int md_param::nb_sub_param(char *name, char *comment) {
    int i = search_param(name), nb;
    
    if (i == -1 || param_type(i) != 4) {
      
      if (i == -1) {
	i = nb_param++;
	param_name(i) = new char[strlen(name)+1];
	strcpy(param_name(i), name);
	param_type(i) = 4; nb_sub_param(i) = 0;
	sub_param(i) = new dal::dynamic_array<t_value>; 
      }
      
      if (param_type(i) != 4) {
	clear_value(i);  param_type(i) = 4;  nb_sub_param(i) = 0;
	sub_param(i) = new dal::dynamic_array<t_value>;
      }
      STD_NEEDED cout << "No list parameter "<< name << " found, please enter its value\n"; 
      STD_NEEDED cout << comment << " : \n";
      STD_NEEDED cout << "Number of sub parameters : "; STD_NEEDED cin >> nb;
      for (int k = 0; k < nb; k++)
	{
	  STD_NEEDED cout << "Value " << k << " : "; STD_NEEDED cin >> string_read;
	  add_string_param_to_param(name, string_read);
	}
    }
    if (param_comment(i) != NULL) delete param_comment(i);
    param_comment(i) = new char[strlen(comment)+1];
    strcpy(param_comment(i), comment);
    return nb_sub_param(i);
  }
  
  long md_param::sub_int_value(char *name, int n) {
    int i = search_param(name);
    if (i == -1) return 0;
    if (param_type(i) != 4) return 0;
    if (nb_sub_param(i) <= n) return 0;
    
    switch((* sub_param(i))[n].type) {
    case 1 : return long( (* sub_param(i))[n].value.v_real);
    case 2 : return  (* sub_param(i))[n].value.v_int;
    case 3 : return atol( (* sub_param(i))[n].value.v_string);
    }
    return 0;
  }
  
  char *md_param::sub_string_value(char *name, int n) {
    int i = search_param(name);
    if (i == -1) return 0;
    if (param_type(i) != 4) return 0;
    if (nb_sub_param(i) <= n) return 0;
    
    switch((* sub_param(i))[n].type) {
    case 1 : sprintf(string_read, "%g", (* sub_param(i))[n].value.v_real);
      return string_read;
    case 2 : sprintf(string_read, "%ld", (* sub_param(i))[n].value.v_int);
      return string_read;
    case 3 : return  (* sub_param(i))[n].value.v_string;
    }
    
    return 0;
  }
  
  double md_param::sub_real_value(char *name, int n) {
    int i = search_param(name);
    if (i == -1) return 0;
    if (param_type(i) != 4) return 0;
    if (nb_sub_param(i) <= n) return 0;
    
    switch((* sub_param(i))[n].type) {
    case 1 : return  (* sub_param(i))[n].value.v_real;
    case 2 : return  double((* sub_param(i))[n].value.v_int);
    case 3 : return atof( (* sub_param(i))[n].value.v_string);
    }
    
    return 0.0;
  }
  
  int md_param::add_string_param(char *name, char *value) {
    int i = search_param(name);
    if (!strcmp(name, "DATA_TYPE")) {
      if (!strcmp(value, "BIN")) is_text = false; else is_text = true;
    }
    if (!strcmp(name, "DATA_FORMAT")) {
      if (!strcmp(value, "float64")) dts = 8; else dts = 4;
    }
    if (i == -1) {
      i = nb_param++;
      param_name(i) = new char[strlen(name)+1];
      strcpy(param_name(i), name);
    }
    else {
      clear_value(i);    
    }
    string_value(i) = new char[strlen(value)+1];
    strcpy(string_value(i), value);
    param_type(i) = 3;
    return i;
  }

  int md_param::add_string_param_to_param(char *name, char *value) {
    int i = search_param(name);
    if (i == -1)
      {
	i = nb_param++;
	param_name(i) = new char[strlen(name)+1];
	strcpy(param_name(i), name);
	param_type(i) = 4; nb_sub_param(i) = 0;
	sub_param(i) = new dal::dynamic_array<t_value>;
      }
    else
      {
	if (param_type(i) != 4)
	  {
	    clear_value(i);  param_type(i) = 4;  nb_sub_param(i) = 0;
	    sub_param(i) = new dal::dynamic_array<t_value>;
	  }  
      }
    (* sub_param(i))[nb_sub_param(i)].value.v_string=new char[strlen(value)+1];
    strcpy((* sub_param(i))[nb_sub_param(i)].value.v_string, value);
    (* sub_param(i))[nb_sub_param(i)].type = 3;
    nb_sub_param(i)++;
    
    return i;
  }
  

  int md_param::add_real_param(char *name, double e) {
    int i = search_param(name);
    if (i == -1)
      {
	i = nb_param++;
	param_name(i) = new char[strlen(name)+1];
	strcpy(param_name(i), name);
      }
    else
      {
	clear_value(i);    
      }
    real_value(i) = e;
    param_type(i) = 1;
    return i;
  }

  int md_param::add_int_param(char *name, long e) {
    int i = search_param(name);
    if (i == -1)
      {
	i = nb_param++;
	param_name(i) = new char[strlen(name)+1];
	strcpy(param_name(i), name);
      }
    else
      {
	clear_value(i);    
      }
    int_value(i) = e;
    param_type(i) = 2;
    return i;
  }
  
  void md_param::wopen_data_file(void) {
    if (f_status != 0) close_data_file();
    temp_blk_size = 100; temp_blk = new char[temp_blk_size*dts];
    flushtime = int_value("DATA_FLUSH_TIME");
    if (flushtime == 0) flushtime = 20;
    write_param_file();
    sprintf(temp_string, "%s%s", string_value("ROOTFILENAME"), ".data");
    fid = fopen(temp_string, "w");
    f_status = 2;
    lblk_count = 0;
    clk = clock() / CLOCKS_PER_SEC;
  }
  
  void md_param::ropen_data_file(void) {
    if (f_status != 0) close_data_file();
    temp_blk_size = 100; temp_blk = new char[temp_blk_size*dts];
    write_param_file();
    sprintf(temp_string, "%s%s", string_value("ROOTFILENAME"), ".data");
    fid = fopen(temp_string, "r");
    f_status = 1;
  }
  
  void md_param::close_data_file(void) {
    fclose(fid); f_status = 0;
    if (temp_blk_size > 0) { delete[] temp_blk; temp_blk_size = 0; }
  }
  
  void  md_param::write_in_data_file(int nb, double *array) {
    if (f_status != 2) wopen_data_file();
    if (is_text)
      {
	for (int i = 0; i < nb; i++)
	  {
	    fprintf(fid, "%g", array[i]);
	    lblk_count++;
	    if (lblk_count == blk1) { lblk_count = 0; fprintf(fid, "\n"); }
	    else fprintf(fid, "\t");
	  }
      }
    else
      {
	if (dts == 4)
	  {
	    if (nb > temp_blk_size)
	      {
		delete[] temp_blk; temp_blk_size=nb;
		temp_blk = new char[temp_blk_size*dts];
	      }
	    for (int i = 0; i < nb; i++)
	      ((float *)temp_blk)[i] = array[i];
	    fwrite(temp_blk, dts, nb, fid);
	  }
	else
	  fwrite(array, dts, nb, fid);
      }
    if ( dal::abs(int(clk - clock() / CLOCKS_PER_SEC)) > flushtime )
      { fflush(fid); clk = clock() / CLOCKS_PER_SEC; }
  }
  
  void  md_param::load_from_data_file(int nb, double *array) {
    if (f_status != 1) ropen_data_file();
    if (is_text)
      {
	for (int i = 0; i < nb; i++)
	  fscanf(fid, "%lg", &(array[i]));
      }
    else
      {
	if (dts == 4)
	  {
	    if (nb > temp_blk_size)
	      {
		delete[] temp_blk; temp_blk_size=nb;
		temp_blk = new char[temp_blk_size*dts];
	      }
	    fread(temp_blk, dts, nb, fid);
	    for (int i = 0; i < nb; i++)
	      array[i] = ((float *)temp_blk)[i];
	  }
	else
	  fread(array, dts, nb, fid);
      }
  }
  
}
