
INCLUDES_GETFEM = -I/mnt-home/gmmpc15/renard/source++/getfem++/include
GETFEM_LIB = /mnt-home/gmmpc15/renard/source++/getfem++/lib/libgetfem.a
COMP = g++-3.0
OPTIONS_STANDARD = -ftemplate-depth-40
OPTIONS_STRICT = -pedantic
OPTIONS_DEBUG    = -g
OPTIONS_OPTIMIZE = -O3
OPTIONS_BOTH     = -g -O3
INCLUDES_STANDARD = -I/mnt/opt/matlab/extern/include 
LIB_X11 = -L/usr/X11R6/lib -lX11
CPLUS = $(COMP) $(OPTIONS_STANDARD) $(OPTIONS_OPTIMIZE) 
CPLUSDEBUG = $(COMP) $(OPTIONS_STANDARD) $(OPTIONS_DEBUG) 
CPLUSOPTDEBUG = $(COMP) $(OPTIONS_STANDARD) $(OPTIONS_BOTH) 
INCLUDES = -I./include $(INCLUDES_STANDARD)  $(INCLUDES_GETFEM) 


LIB_X11 = -L/usr/X11R6/lib -lX11


EXECUTABLES = laplacien
# CPLUS = $(CPLUSOPTDEBUG)
# CPLUS = $(CPLUSDEBUG)
CPLUS = g++ -ainsi -ftemplate-depth-40
all : lib_outil $(EXECUTABLES)
INC = $(INCLUDES) -I../include $(INCLUDES)

lib_outil :
	cd .. ; $(MAKE)

libraries : $(PERSO_LIB)

test_mesh : test_mesh.C $(GETFEM_LIB)
	$(CPLUS) -g $(COPT) $(LIBDIR) $(INC) test_mesh.C ../sources/*.o $(libraries) $(LIB) $(LIB_X11) -o test_mesh

laplacien : laplacien.C $(GETFEM_LIB) $(GETFEM)
	$(CPLUS) -g $(COPT) $(INC) laplacien.C ../lib/libgetfem.a $(libraries) $(LIB) $(LIB_X11) -I/opt/mtl/include -o laplacien

elasto_statique : elasto_statique.C $(GETFEM_LIB) $(GETFEM)
	$(CPLUS) -g $(COPT) $(LIBDIR) $(INC) elasto_statique.C ../sources/*.o /opt/getfem/lib/libgetfem.a $(libraries) $(LIB) $(LIB_X11) -o elasto_statique

clean :
	/bin/rm -f *.o *.log *.dvi *% *.ps *.aux *.bak *~ *#
	/bin/rm -f $(EXECUTABLES) *.dat *.fig
