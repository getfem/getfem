/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_except.C                                                 */
/*     									   */
/*                                                                         */
/* Date : September 01, 2002.                                              */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, pommier@gmm.insa-tlse.fr                       */
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
#include <dal_except.h>

namespace dal {
  static exception_callback *exc_cback = 0;
  
  void do_exception_callback(const std::string &msg) {
    if (exc_cback) {
      exc_cback->callback(msg);
    }
  }
  
  void set_exception_callback(exception_callback *e) {
    exc_cback = e;
  }

  /* juste pour forcer l'instanciation dans libgetfem de tous les
     type qu'on est suceptible afficher. Si il en manque un (c'est le cas pour unsigned et float)
     dans libgetfem.so qui est utilise dans libgetfemint.so alors on retrouve un warning
     comme quoi le symbole __T_Q13std8ios_base et __T_Q13std11logic_error sont dupliqués
     visiblement il reinstancie un peu trop de trucs .. 

     pour s'assurer que le probleme n'est pas de retour, faire
     nm libgetfem.so | grep 'td::basic_ostream<char, std::char_traits<char> >::operator <<'
     et faire de même dans tous les .o de matlabint, et s'assurer que tous ceux de matlabint
     sont bien inclus dans libgetfem.so

     enfin c'est pas tout a fait ça puisque tous ces trucs sont deja instancies dans libcxxstd.a, et c'est le noeud du probleme j'imagine:
     nm /usr/lib/cmplrs/cxx/V6.3-008/libcxxstd.a | grep 'td::basic_ostream<char, std::char_traits<char> >::operator <<' 
     ---> std::basic_ostream<char, std::char_traits<char> >::operator <<(long long) | 0000000000006368 | T | 0000000000000008
          std::basic_ostream<char, std::char_traits<char> >::operator <<(const void*) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(std::basic_ostream<char, std::char_traits<char> >& (*)(std::basic_ostream<char, std::char_traits<char> >&)) | 0000000000000000 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(unsigned long long) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(unsigned int) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(unsigned long) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(unsigned short) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(double) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(float) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(int) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(bool) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(long double128) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(long) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(short) | 0000000000006368 | T | 0000000000000008
	  std::basic_ostream<char, std::char_traits<char> >::operator <<(long double64) | 0000000000006368 | T | 0000000000000008
  */
  void just_for_the_fine_cxx() {
    std::stringstream s; s << int(1) << double(2.0) 
    << "hello" << std::string("hello") << unsigned(1) << float(2) 
    << char('a') << (unsigned char)('b') << short(1) << (unsigned short)(2) 
			   << long(1) << (unsigned long)(2) << (const void*)NULL;
  }
}
