
include /opt/getfem/getfem++/defmakefile


EXECUTABLES = elasto_statique
CPLUS = $(CPLUSOPTDEBUG)
# CPLUS = $(CPLUSDEBUG)
# CPLUS = g++ -ainsi -ftemplate-depth-40
all : lib_outil $(EXECUTABLES)
INC += -I../include $(INCLUDES)

lib_outil :
	cd .. ; $(MAKE)

libraries : $(PERSO_LIB)

test_mesh : test_mesh.C $(GETFEM_LIB) $(GETFEM)
	$(CPLUS) -g $(COPT) $(LIBDIR) $(INC) test_mesh.C ../sources/*.o $(libraries) $(LIB) $(LIB_X11) -o test_mesh

laplacien : laplacien.C $(GETFEM_LIB) $(GETFEM)
	$(CPLUS) -g $(COPT) $(LIBDIR) $(INC) laplacien.C ../sources/*.o /opt/getfem/lib/libgetfem.a $(libraries) $(LIB) $(LIB_X11) -o laplacien

elasto_statique : elasto_statique.C $(GETFEM_LIB) $(GETFEM)
	$(CPLUS) -g $(COPT) $(LIBDIR) $(INC) elasto_statique.C ../sources/*.o /opt/getfem/lib/libgetfem.a $(libraries) $(LIB) $(LIB_X11) -o elasto_statique

clean :
	/bin/rm -f *.o *.log *.dvi *% *.ps *.aux *.bak *~ *#
	/bin/rm -f $(EXECUTABLES) *.dat *.fig
