#transforme le contenu du rep CBLAS de superlu en un seul fichier c
import sys
import re

f = open('BLAS_c','w')
f.write('#include "f2c_lite.h"\n')

for fname in sys.argv[1:]:
	cf = open(fname);
	defines = []
	for l in cf.readlines():
		if (l.startswith('#include')):
			continue
		
		if (l.find('#define') != -1):
			m = re.search('#define *([A-Za-z0-9_]*)',l)
			if (m):
				print l,
				defines += [m.group(1)]
		f.write(l)		
	for d in defines:
		l = '#undef ' + d + '\n'
		print l,
		f.write(l)
