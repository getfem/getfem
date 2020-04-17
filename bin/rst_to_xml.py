#!/usr/bin/env python
# -*- python -*-
#
# Copyright (C) 2010-2020 Yves Renard.
#                                                       
# This file is a part of GetFEM                                         
#                                                                         
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#
############################################################################
"""  Transform a rst file into a xml one.

  xml2rst is used for the text part and tralics for the math formulaes.

  $Id: extract_doc 3304 2009-11-03 13:17:46Z renard $
"""
import os
import string
import sys


class ParseError(Exception):
    def __init__(self, value):
      self.value = value
    def __str__(self):
      return repr(self.value)


if (len(sys.argv) != 2):
    raise SystemExit, 'Format : rst_to_xml filename'

filename = sys.argv[1]

fl = open(filename)
temprst = open(filename+'_temp.rst', 'w')
in_math_mode = 0
count_math_f = 0
ntab = 0
math_forms = []


# read the file and detect the ..math:: and :math: replace it by some tags
# and store the formulaes
l = fl.readline()
while(len(l)):
    ll = l.strip()
    if (in_math_mode == 0):
        if (ll[0:2] == '..' and ll[2:].strip()[0:6] == 'math::'):
            in_math_mode = 1
            math_form = ''
        elif (ll.find(':math:') != -1):
            orgl = l
            j = l.find(':math:')
            while (j != -1):
                temprst.write(l[:j])
                l = l[j+6:].strip()
                if (l[0] != '`'): raise ParseError, orgl
                l = l[1:].strip()
                j = l.find('`')
                math_form = ''
                while (j == -1):
                    math_form += ' ' + l
                    l = fl.readline()
                    if (not len(l)): raise ParseError, 'Reach end of file'
                    l = l.strip()
                    j = l.find('`')
                math_form += ' ' + l[:j]
                math_forms.append('$'+math_form+'$')
                count_math_f += 1
                temprst.write("MATHZFORMULE%06d" % count_math_f)
                l = l[j+1:]
                j = l.find(':math:')
            temprst.write(l+'\n')
            l = ''
    elif (in_math_mode == 1 and ll != ''):
        math_form += ll
        for i in range(len(l)):
            ntab = i;
            if (not l[i].isspace()): break
        in_math_mode = 2
    elif (in_math_mode == 2 and ll == ''):
        if (math_form != ''):
            count_math_f += 1
            temprst.write("MATHZFORMULE%06d" % count_math_f)
            math_forms.append('$$'+math_form+'$$')
        math_form = ''
    elif (in_math_mode == 2 and ll != ''):
        for i in range(len(l)):
            nntab = i;
            if (not l[i].isspace()): break
        if (nntab == ntab):
           math_form += ll
        else:
           in_math_mode = 0
           if (math_form != ''):
               count_math_f += 1
               temprst.write("MATHZFORMULE%06d\n" % count_math_f)
               math_forms.append('$$'+math_form+'$$')
 
    if (in_math_mode == 0):
        temprst.write(l)
    l = fl.readline()
    

temprst.close()
fl.close()


math_forms_trans = []

for iform in range(count_math_f):
    temprst = open(filename+'_temp_f.tex', 'w')
    math_form = math_forms[iform];
    math_form = math_form.replace('\\mathscr', '\\cal')
    print(math_form)
    if (math_form.count('&')):
        temprst.write('\\begin{eqnarray*}\n')
        temprst.write(math_form[2:len(math_form)-2] + '\n')
        temprst.write('\\end{eqnarray*}\n')
    else:
        temprst.write(math_form)
    temprst.close()
    if (os.system('tralics ' + filename+'_temp_f.tex')): exit(1)
    fl = open(filename+'_temp_f.xml')
    for l in fl:
        if (l[:13] == '<formula type'):
            math_forms_trans.append(l)
            print(("Formule %d : " % iform) + l)
            break
        if (l[:16] == '<p><formula type'):
            math_forms_trans.append(l[3:])
            print(("Formule %d : " % iform) + l)
            break
    fl.close()


fl = os.popen('rst2xml ' + filename+'_temp.rst')
rfl = open(filename+'.xml', 'w')
for l in fl:
    if (l.find("MATHZFORMULE") != -1):
        j = l.find("MATHZFORMULE")
        while (j != -1):
          r = l[j+12:j+18]
          print(r)
          nf = int(r)
          print(nf)
          print("MATHZFORMULE%06d" % nf)
          print(math_forms_trans[nf-1])
          l = string.replace(l, ("MATHZFORMULE%06d" % nf), math_forms_trans[nf-1])
          print(l)
          j = l.find("MATHZFORMULE")
    rfl.write(l)
rfl.close()
    


print("there were ", count_math_f, " formulaes")
