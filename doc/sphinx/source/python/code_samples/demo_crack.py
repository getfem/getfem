#!/usr/bin/env python
# -*- coding: utf8 -*-
"""  Linear Elastostatic problem with a crack.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

"""

##########################################################################
#  Exact solution.                                                       #
##########################################################################

tol = 0.0001

def sint2(x,y):
  """.
  returns sin(theta/2) where theta is the angle of 0-(x,y) with the axis Ox
  """
  r = sqrt(x*x+y*y)
  if r < tol:
    return 0
  elif y<0:
    return -sqrt(abs(r-x)/(2*r))
  return sqrt(abs(r-x)/(2*r))

def cost2(x,y):
  """.
  returns cos(theta/2) where theta is the angle of 0-(x,y) with the axis 
  Ox
  """
  r = sqrt(x*x+y*y)
  if r < tol:
    return 0
  return sqrt(abs(r+x)/(2*r))

#
# analytical solution for a semi-infinite crack [-inf,a] in an
# infinite plane submitted to +sigma above the crack
# and -sigma under the crack. (The crack is directed along the x axis).
#
# nu and E are the poisson ratio and young modulus
#
# solution taken from "an extended finite elt method with high order
# elts for curved cracks", Stazi, Budyn,Chessa, Belytschko
#

def elasticite2lame(young_modulus, poisson_ratio):
  """.
  returns lamé coeficients (lambda, mu) 
  Ox
  """
  mu = young_modulus/(2*(1+poisson_ratio))
  la = 2*mu*poisson_ratio/(1-poisson_ratio)
  return (la,mu)
