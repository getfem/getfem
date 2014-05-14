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
 
===========================================================================*/


#include "getfem/bgeot_config.h"
#include "getfem/bgeot_ftool.h"
#include <ctype.h>
#include <limits.h>
#ifndef _WIN32
#  include <unistd.h>
#endif
#include <fstream>

namespace bgeot {

  bool read_until(std::istream &ist, const char *st) {
    int i = 0, l = int(strlen(st)); char c;
    while (!ist.eof() && i < l)
      { ist.get(c); if (toupper(c) == toupper(st[i])) i++; else i = 0; }
    if (ist.eof()) return false; else return true;
  }
  
#define get_c__(r, c) {	ist.get(c); if (ist.eof()) return r;  \
    if (to_up) c = char(toupper(c)); }

#define sdouble__(c, e) {  st.push_back(c); get_c__(5, d); \
    if (d == e) { st.push_back(e); return 6; }		   \
    else { ist.putback(d); return 5; } }		   \

  int get_token(std::istream &ist, std::string &st,
		bool ignore_cr, bool to_up, bool read_un_pm, int *linenb) {
    st.resize(0);
    char c = char(-1), d, e;
   
    get_c__(0, c);

    for(;;) { // Go through spaces, commentaries and '...'
      if (!ignore_cr && c == '\n') { if (linenb) (*linenb)++; return 1; }
      if (isspace(c)) { while (isspace(c)) get_c__(0, c); }
      else if (c == '%') { while (c != '\n') get_c__(0, c); }
      else if (c == '.') {
	if (ist.eof()) break; else {
	  get_c__(0, d);
	  if (d == '.'  && !ist.eof()) {
	    get_c__(0, e);
	    if (e == '.') {
	      while (c != '\n') get_c__(0, c);
	      if (linenb) (*linenb)++; 
	      get_c__(0, c);
	    }
	    else { ist.putback(e); ist.putback(d); break; }
	  }
	  else { ist.putback(d); break; }
	}
      }
      else break;
    }

    if (read_un_pm)
      if (c == '-' || c == '+') { // reading a number beginning with '+' or '-'
	get_c__(2, d);
	if (isdigit(d) || d == '.') { st.push_back(c); c = d; }
	else ist.putback(d);
      }

    if (isdigit(c) || c == '.') { // reading a number
      while (isdigit(c) || c == '.' || c == 'e'  || c == 'E') {
	st.push_back(c); 
	if (c == 'e' || c == 'E') {
	  get_c__(2, c);
	  if (c == '+' || c == '-') st.push_back(c);
	  else ist.putback(c);
	} 
	get_c__(2, c);
      }
      ist.putback(c);
      return 2;
    }

    if (c == '\"') { // reading a string
      get_c__(3, c);
      while (true) {
	if (c == '\"' || c == '\n') return 3;
	if (c == '\\') { st.push_back(c); get_c__(3, c); }
	st.push_back(c);
        get_c__(3, c);
      }
      return 3;
    }

    if (c == '\'') { // reading a string
      get_c__(3, c);
      while (true) {
	if (c == '\'' || c == '\n') return 3;
	if (c == '\\') { st.push_back(c); get_c__(3, c); }
	st.push_back(c);
        get_c__(3, c);
      }
      return 3;
    }

    if (isalpha(c) || c == '_') { // reading a name
      while (isalnum(c) || c == '_') { st.push_back(c); get_c__(4,c); }
      ist.putback(c);
      return 4;
    }

    if (c == '|') sdouble__(c, '|');
    if (c == '&') sdouble__(c, '&');
    if (c == '=') sdouble__(c, '=');
    if (c == '~') sdouble__(c, '=');
    if (c == '<') sdouble__(c, '=');
    if (c == '>') sdouble__(c, '=');   

    st.push_back(c); return 5; // return the symbol read.
  }

  std::istream& operator>>(std::istream& is, const skip& t) {
    char c;
    int i = 0;
    while (!is.get(c).eof() && isspace(c)) /*continue*/;    
    for (i=0; t.s[i]; ++i) {
      if (i) is.get(c);
      GMM_ASSERT1(toupper(c) == toupper(t.s[i]) && !is.eof(),
		  "expected token '"<<t.s<<"' not found");
    }
    return is;
  }

  int casecmp(const char *a, const char *b, unsigned n) {
    unsigned i;
    for (i=0; i < n && a[i] && b[i]; ++i) {
      if (toupper(a[i]) < toupper(b[i])) return -1;
      else if (toupper(a[i]) > toupper(b[i])) return -1;
    }
    if (a[i]) return +1;
    else if (b[i]) return -1;
    else return 0;
  }
  
  void md_param::parse_error(const std::string &t) {
    GMM_ASSERT1(false, "Parse error reading "
		<< current_file << " line " << current_line << " near " << t);
  }

  void md_param::syntax_error(const std::string &t) {
    GMM_ASSERT1(false, "Error reading "
		<< current_file << " line " << current_line << " : " << t);
  }

  int md_param::get_next_token(std::istream &f) {
    static int token_type = 0;
    if (!token_is_valid)
      token_type = get_token(f, temp_string, false, false, false,
			     &current_line);
    token_is_valid = false;
    return token_type;
  }

  void md_param::valid_token(void) { token_is_valid = true; }

  std::ostream &operator <<(std::ostream &o, const md_param::param_value& p) {
    switch (p.type_of_param()) {
    case md_param::REAL_VALUE : o << p.real(); break;
    case md_param::STRING_VALUE : o << '\'' << p.string() << '\''; break;
    case md_param::ARRAY_VALUE : 
      o << "[";
      if (p.array().size()) o << p.array()[0];
      for (unsigned i = 1; i < p.array().size(); ++i)
	o << ", " << p.array()[i];
      o << "]";
    }
    return o;
  }

  md_param::param_value md_param::read_expression(std::istream &f,
						  bool skipped) {
    param_value result;
    int i = get_next_token(f);
    if (i == 2) { // a number
      result = param_value(::strtod(temp_string.c_str(), 0));
    }
    else if (i == 3) { // a string
      result = param_value(temp_string);
      int j = get_next_token(f);
      while (j == 3) {
	result.string() += temp_string;
	j = get_next_token(f);
      }
      valid_token();
    }
    else if (i == 4) { // a parameter name
      std::string name(temp_string);
      if (parameters.find(name) != parameters.end())
	result = parameters[name];
      else if (!skipped) {
	std::stringstream s; s << "Parameter " << name << " not found";
	syntax_error(s.str());
      }
    }
    else if (i == 5) { // unary operators, parentheses and arrays
      switch (temp_string[0]) {
      case '(' :
	{
	  result = read_expression_list(f, skipped);
	  int j = get_next_token(f);
	  if (j != 5 || temp_string[0] != ')') parse_error(temp_string);
	}
	break;
      case '+' :
	result = read_expression(f, skipped);
	if (result.type_of_param() != REAL_VALUE)
	  syntax_error("Sorry, unary + does not support string "
		       "or array values");
	break;
      case '-' :
	result = read_expression(f, skipped);
	if (result.type_of_param() != REAL_VALUE)
	  syntax_error("Sorry, unary - does not support string "
			 "or array values");
	result.real() *= -1.0;
	break;
      case '~' : 
	result = read_expression(f, skipped);
	if (result.type_of_param() != REAL_VALUE)
	  syntax_error("Sorry, unary ! does not support string "
			 "or array values");
	result.real() = !(result.real());
	break;
      case '[' :
	{
	  bool first = true;
	  result = param_value(ARRAY_VALUE);
	  while (true) {
	    int j = get_next_token(f);
	    if (j == 5 && temp_string[0] == ']') break;
	    if (!first && temp_string[0] != ',') parse_error(temp_string);
	    if (first) valid_token();
	    result.array().push_back(read_expression_list(f, skipped));
	    first = false;
	  }
	}
	break;
      default : parse_error(temp_string);
      }
    }
    else parse_error(temp_string);

    return result;
  }

  static void operator_priority_ftool(int i, char c, int &prior, int &op) {
    if (i == 5)
      switch (c) {
      case '*' : prior = 1; op = 1; return;
      case '/' : prior = 1; op = 2; return;
      case '+' : prior = 2; op = 3; return;
      case '-' : prior = 2; op = 4; return;
      case '<' : prior = 3; op = 5; return;
      case '>' : prior = 3; op = 6; return;
      }
    if (i == 6)
      switch (c) {
      case '<' : prior = 3; op =  7; return; // <= 
      case '>' : prior = 3; op =  8; return; // >= 
      case '=' : prior = 3; op =  9; return; // == 
      case '~' : prior = 3; op = 10; return; // != 
      case '&' : prior = 4; op = 11; return; // && 
      case '|' : prior = 4; op = 12; return; // ||
      }
    prior = op = 0;
  }

  void md_param::do_bin_op(std::vector<md_param::param_value> &value_list,
			std::vector<int> &op_list,
			std::vector<int> &prior_list) {
    param_value &p1(*(value_list.end() - 2));
    param_value &p2(*(value_list.end() - 1));
    if (p1.type_of_param() != REAL_VALUE || p2.type_of_param() != REAL_VALUE)
      syntax_error("Sorry, binary operators does not support string "
		     "or array values");
    
    switch (op_list.back()) {
    case 1  : p1.real() *= p2.real(); break;
    case 2  : p1.real() /= p2.real(); break;
    case 3  : p1.real() += p2.real(); break;
    case 4  : p1.real() -= p2.real(); break;
    case 5  : p1.real() = (p1.real() < p2.real()); break;
    case 6  : p1.real() = (p1.real() > p2.real()); break;
    case 7  : p1.real() = (p1.real() <= p2.real()); break;
    case 8  : p1.real() = (p1.real() >= p2.real()); break;
    case 9  : p1.real() = (p1.real() == p2.real()); break;
    case 10 : p1.real() = (p1.real() != p2.real()); break;
    case 11 : p1.real() = ((p1.real() != 0.0) && (p2.real() != 0.0)); break;
    case 12 : p1.real() = ((p1.real() != 0.0) || (p2.real() != 0.0)); break;
    }
    value_list.pop_back(); op_list.pop_back(); prior_list.pop_back();
  }


  md_param::param_value md_param::read_expression_list(std::istream &f,
						       bool skipped) {
    std::vector<param_value> value_list;
    value_list.push_back(read_expression(f, skipped));
    std::vector<int> op_list, prior_list;
    int i = get_next_token(f), prior, op;
    operator_priority_ftool(i, temp_string[0], prior, op);
    while (op) {
      while (!prior_list.empty() && prior_list.back() <= prior)
	do_bin_op(value_list, op_list, prior_list);

      value_list.push_back(read_expression(f, skipped));
      op_list.push_back(op);
      prior_list.push_back(prior);

      i = get_next_token(f);
      operator_priority_ftool(i, temp_string[0], prior, op);
    }
    valid_token();

    while (!prior_list.empty()) do_bin_op(value_list, op_list, prior_list);

    return value_list[0];
  }

  int md_param::read_instruction(std::istream &f, bool skipped) {
    int i = 1;
    while (i == 1 || (i == 5 && temp_string[0] == ';')) i = get_next_token(f);
    if (i == 0) return 1;
    if (i != 4) parse_error(temp_string);
    if (temp_string == "end") return 1;
    if (temp_string == "else") return 2;
    if (temp_string == "elseif") return 3;
    if (temp_string == "if") {
      param_value p = read_expression_list(f, skipped);
      if (p.type_of_param() != REAL_VALUE)
	syntax_error("if instruction needs a condition");
      bool b = (p.real() != 0.0);
      int j = read_instruction_list(f, !b || skipped);
      if (j == 0) syntax_error("Unterminated if");
      if (j == 2) {
	int k = read_instruction_list(f, b || skipped);
	if (k != 1) syntax_error("Unterminated else");
      }
      if (j == 3) {
	int k = 0;
	do {
	  if (b) skipped = true;
	  p = read_expression_list(f, skipped);
	  if (p.type_of_param() != REAL_VALUE)
	    syntax_error("elseif instruction needs a condition");
	  b = (p.real() != 0.0);
	  k = read_instruction_list(f, !b || skipped);
	  if (k == 2) {
	    k = read_instruction_list(f, b || skipped);
	    break;
	  }
	} while (k == 3);
	if (k != 1) syntax_error("Unterminated elseif");
      }     
      return 0;
    }
    if (temp_string == "error") {
      param_value p = read_expression_list(f, skipped);
      GMM_ASSERT1(skipped, "Error in parameter file: " << p);
      return 0;
    }
    std::string name(temp_string);
    i = get_next_token(f);
    if (i != 5 || temp_string[0] != '=') parse_error(temp_string);
    param_value result = read_expression_list(f, skipped);
    i = get_next_token(f);
    if (i != 0 && i != 1 && (i != 5 || temp_string[0] != ';'))
      parse_error(temp_string);
    if (!skipped) parameters[name]=result;
    return 0;
  }

  int md_param::read_instruction_list(std::istream &f, bool skipped) {
    int i; while (!(i = read_instruction(f, skipped))) { }
    return i;
  }

  void md_param::read_param_file(std::istream &f) {
    gmm::standard_locale sl;
    token_is_valid = false; current_line = 1;
    if (read_instruction_list(f) > 1)
      syntax_error("Parameter file terminated by an else");
  }
  
  void md_param::read_command_line(int argc, char *argv[]) {
    gmm::standard_locale sl;
    for (int aa = 1; aa < argc; aa++) {
      if (argv[aa][0] != '-') {
	current_file = std::string(argv[aa]);
	std::ifstream f1(current_file.c_str());
	if (f1) { read_param_file(f1); f1.close(); }
	else {
	  std::string r = current_file;
	  current_file += ".param";
	  std::ifstream f2(current_file.c_str());
	  if (f2) { read_param_file(f2); f2.close(); }
	  else GMM_ASSERT1(false,  "Parameter file " << r << "not found");
	}
      }
      else if (argv[aa][1] == 'd') {
	current_file = "command line";
	if (strlen(argv[aa]) == 2)
	  { std::stringstream ss(argv[++aa]); read_param_file(ss); }
	else 
	  { std::stringstream ss(&(argv[aa][2])); read_param_file(ss); }
      }
    }
  }
  
  double md_param::real_value(const std::string &name, const char *comment) {
    if (parameters.find(name) == parameters.end()) {
      if (comment == 0) return 0.0;
      else {
	double f;
	gmm::standard_locale sl;
	cout << "No parameter " << name << " found, please enter its value\n";
	cout << comment << " : "; cin >> f;
	parameters[name] = param_value(f);
      }
    }
    param_value &p(parameters[name]);
    GMM_ASSERT1(p.type_of_param() == REAL_VALUE,
		"Parameter " << name << " is not real");
    return p.real();
  }
  
  long md_param::int_value(const std::string &name, const char *comment) {
    if (parameters.find(name) == parameters.end()) {
      if (comment == 0) return 0;
      else {
	long f;
	gmm::standard_locale sl;
	cout << "No parameter " << name << " found, please enter its value\n";
	cout << comment << " : "; cin >> f;
	parameters[name] = param_value(double(f));
      }
    }
    param_value &p(parameters[name]);
    GMM_ASSERT1(p.type_of_param() == REAL_VALUE,
		"Parameter " << name << " is not real");
    return long(p.real());
  }
  
  const std::string &md_param::string_value(const std::string &name,
				     const char *comment) {
    static const std::string empty_string;
    if (parameters.find(name) == parameters.end()) {
      if (comment == 0) return empty_string;
      else {
	std::string s;
	gmm::standard_locale sl;
	cout << "No parameter " << name << " found, please enter its value\n";
	cout << comment << " : "; cin >> s;
	parameters[name] = param_value(s);
      }
    }
    param_value &p(parameters[name]);
    GMM_ASSERT1(p.type_of_param() == STRING_VALUE, "Parameter " << name
		<< " is not a character string");
    return p.string();
  }
  
  const std::vector<md_param::param_value> &
  md_param::array_value(const std::string &name, const char *comment) {

    static std::vector<md_param::param_value> empty_array;
    if (parameters.find(name) == parameters.end()) {
      if (comment == 0) return empty_array;
      else {
	std::string s;
	gmm::standard_locale sl;
	cout << "No parameter " << name << " found, please enter its value\n";
	cout << comment << " : "; cin >> s;
	parameters[name] = param_value(s);
      }
    }
    param_value &p(parameters[name]);
    GMM_ASSERT1(p.type_of_param() == ARRAY_VALUE, "Parameter " << name
		<< " is not an array");
    return p.array();
  }
}
