#!/usr/bin/python
import re
import glob

def doxrename(f):
    latexmacro = re.sub('[._:]','',f)
    latexmacro = re.sub('0','zero',latexmacro)
    latexmacro = re.sub('1','one',latexmacro)
    latexmacro = re.sub('2','two',latexmacro)

    doxname = re.sub('_','__',f)
    doxname = re.sub('\.','_8',doxname)
    doxname = re.sub(':','_1',doxname)

    escapedname = re.sub('_', '\\_', f)
    
    return (latexmacro,doxname,escapedname)

flist=glob.glob1('../../src/', '*.h') + glob.glob1('../../src/', '*.cc')

out=file('doxygenlinks.tex','wt')

for f in flist:
    n = doxrename(f)
    print "doing file %s" % (n,)
    out.write('\\newcommand{\\%s}{\\doxfilename{%s}{%s}}\\xspace\n' % (n[0],n[2], n[1]))

classes="""
dal::bit_vector
dal::bv_visitor
bgeot::convex_structure/bgeot::pconvex_structure
bgeot::convex_ref/bgeot::pconvex_ref
bgeot::geometric_trans/bgeot::pgeometric_trans
getfem::virtual_fem/getfem::pfem
getfem::mesh
getfem::mesh_region
getfem::mr_visitor
bgeot::mesh_structure
getfem::mesh_im
getfem::mesh_fem
getfem::stored_mesh_slice
getfem::dx_export
getfem::vtk_export
"""

for c in classes.split():
    ll=c.split('/')
    n = doxrename(ll[0])

    print "doing class %s" % (n,)

    out.write('\\newcommand{\\%s}{\\doxref{%s}{class%s}}\\xspace\n' % (n[0],n[2], n[1]))

    if (len(ll)>1):
        m = doxrename(ll[1])
        print "doing alias %s" % (ll[1],)
        out.write('\\newcommand{\\%s}{\\doxref{%s}{class%s}}\\xspace\n' % (m[0],m[2], n[1]))
