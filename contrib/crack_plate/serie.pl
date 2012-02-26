# Copyright (C) 2001-2012 Jeremie Lasry
#
# This file is a part of GETFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

use strict ;

# This program can execute a serie of getfem calculus. The different steps are :
#  1. Definition of the .param file
#  2. Throwing calculs and get the result
#  3. write the result on a file.
# Only the meshes are put on a serie, but extensions are widely possible.

# Parameters to set
my $elem_type = "quad" ;
my $mixed_elements = 1 ; # important : set to 1 in order to split quadrangles cut by the crack.  
my $xfem_enrichment = 3 ;
my $mortar_type = 3 ;
my $fichier_entree = "serie.param" ;
my $fichier_sortie = "res.m" ;
my $sing_base_type = 0 ;  # set to 0 for singuls on 4 dofs, set to 1 for singuls on 2 dofs
my $mesh_files_or_nx = 0 ; # set to 0 for some getfem mesh files, set to 1 for some mesh files,
my $mesh_noised = 1 ;     # set to 1 if you want to "shake" the mesh 
my $mesh_directory = "/home-local/jlasry/RESULTATS_XFEM/KL/MAILLAGES_MATLAB/" ;
#my @mesh_names = ("10", "12", "15", "21", "30", "39", "60");
my @mesh_names = ("11", "15", "21", "31", "41", "51", "61", "71", "81", "101" ) ;
#my @mesh_names = ("9", "11") ;
my $radius_enrich = 0.2 ;

# Nothing to modify beyond this line.
# Last comment : printing in files must be "cleaned"
# --------------------------------------------------------------------------

my (@cond, @err_L2, @err_H1, @err_H2, @fic1, @fic2, @fic3, @fic4) ;


#initialisations
foreach (@mesh_names){
		@cond[scalar(@cond)] = "null" ;
		@err_L2[scalar(@cond)] = "null" ;
		@err_H1[scalar(@cond)] = "null" ;
		@err_H2[scalar(@cond)] = "null" ;
		@fic1[scalar(@fic1)] = "null" ;
		@fic2[scalar(@fic2)] = "null" ;
		@fic3[scalar(@fic3)] = "null" ;
		@fic4[scalar(@fic4)] = "null" ;
}

# print "err_L2 =" ;
# foreach(@err_L2) {
#    print $_ ;}	
# print "\n" ;

sub start_program {
  my $def   = $_[0];
  my $nb_fic = 0 ;

  print "def = $def\n";

  my @result = ("null", "null", "null","null","null","null","null","null");
  my $var = "null" ;
  open F, "./crack_bilaplacian $fichier_entree $def 2>&1|" or die("crack_bilaplacian not found");
  while (<F>) {
    if ($_ =~ /condition number/) {
      ($a, $var) = split(':', $_); 
      @result[0] = $var ;
    }
    if ($_ =~ /L2 ERROR/) {
      ($a, $var) = split(':', $_); 
      @result[1] = $var ;
    }
    if ($_ =~ /H1 ERROR/) {
      ($a, $var) = split(':', $_); 
      @result[2] = $var ;
    }
    if ($_ =~ /H2 ERROR/) {
      ($a, $var) = split(':', $_); 
      @result[3] = $var ;
    }
    if ($_ =~ /FIC/) {
      ($a, $var) = split(':', $_); 
      @result[4 + $nb_fic] = $var ;
      $nb_fic += 1 ;
    }

  }
  close(F);
  if ($?) {
    #`rm -f $tmp`; 
    print "./crack_bilaplacian failed\n";
    exit(1);
  }
  return @result;
}


# main -----------------------------------------------



# 1. Cr�ation du fichier .param

#my $i = 0 ;
#my $j = 0 ;
my $k = 0 ;

foreach (my $i = 0 ; $i< scalar(@mesh_names) ; $i++) {
        #$j = 0 ; -> etrangement, on dirait qu'on ne peut pas acceder 2 fois de suite � la m�me case d'un tableau.
		print "i = ", $i, "\n" ;
		print "mailage : ", @mesh_names[$i], "\n" ;
		
		open(PARAM, "+>".$fichier_entree) or die("opening".$fichier_entree.$!);
		# ouverture du fichier en lecture et ecriture 
			# (ecriture => ecrasement => peu subtil...)
			
		print PARAM  
		"%%%  PLATE PARAMETERS %%% \n",
		"TRANSLAT_X = 0.0 ;\n" ,
		"TRANSLAT_Y = 0.0 ;\n" ,
		"FT = 5.0 ;\n" ,
		"D = 1.0 ;\n" ,     
		"KL = 1 ;\n" ,
		"NU = 0.3 ;\n" ,
		"EPSILON=0.045 ;\n \n" ;
		
		# Setting the mesh parameter (same for tri or quad)
		my $quad = 0 ;
		if ($elem_type eq "quad") { $quad = 1 ; }
		print PARAM 
		"%%% MESH & FEMS PARAMETERS %%% \n" ;
		if ($mesh_files_or_nx eq 0){
		   print PARAM
		   "LX = 1.0; LY = LX; LZ = LX; \n",
		   "N = 2; \n",
		   "MESH_NOISED = ".$mesh_noised." ; \n",
		   "MIXED_ELEMENTS = ".$mixed_elements." ; \n",        
		   "NX = ".@mesh_names[$i]." ; \n". 
		   "QUAD = ".$quad." ; \n" ;
		}
		if ($mesh_files_or_nx eq 1){
			print PARAM
			"MESH_FILE = '".$mesh_directory."tri_".@mesh_names[$i].".mesh' ;\n" ;
		}
		if ($mixed_elements eq 0){
			# instructions that are different wether using triangles or quadrangles :
			if ($elem_type eq "tri"){
			print PARAM
			"MESH_TYPE = 'GT_PK(2,1)'; \n",   
			"DATA_FEM_TYPE = 'FEM_PK(2, 5)';\n",  
			"PARTITION_OF_UNITY_FEM_TYPE = 'FEM_REDUCED_HCT_TRIANGLE';\n",  
			"FEM_TYPE = 'FEM_REDUCED_HCT_TRIANGLE';\n",  
			"DIRICHLET_FEM_TYPE = 'FEM_PK(2,2)';\n",  
			"DIRICHLET_DER_FEM_TYPE = 'FEM_PK(2,1)';\n",  
			"INTEGRATION = 'IM_HCT_COMPOSITE(IM_TRIANGLE(13))';\n",  
			"MORTAR_FEM_TYPE = 'FEM_PK(2,2)';\n",  
			"MORTAR_DERIV_FEM_TYPE = 'FEM_PK(2,1)';\n \n";
			}
			elsif ($elem_type eq "quad"){
			print PARAM
			"MESH_TYPE = 'GT_QK(2,1)'; \n",   
			"DATA_FEM_TYPE = 'FEM_QK(2, 5)';\n",  
			"PARTITION_OF_UNITY_FEM_TYPE = 'FEM_REDUCED_QUADC1_COMPOSITE';\n",  
			"FEM_TYPE = 'FEM_REDUCED_QUADC1_COMPOSITE';\n",  
			"DIRICHLET_FEM_TYPE = 'FEM_QK(2,2)';\n",  
			"DIRICHLET_DER_FEM_TYPE = 'FEM_QK(2,1)';\n",  
			"INTEGRATION = 'IM_QUADC1_COMPOSITE(IM_TRIANGLE(13))';\n",  
			"MORTAR_FEM_TYPE ='FEM_QK(2,2)';\n",  
			"MORTAR_DERIV_FEM_TYPE = 'FEM_QK(2,1)';\n \n";
			} 
		}
		else{
		print PARAM
		"TRI_MESH_TYPE = 'GT_PK(2,1)'; \n",   
		"TRI_DATA_FEM_TYPE = 'FEM_PK(2, 5)';\n",  
		"TRI_PARTITION_OF_UNITY_FEM_TYPE = 'FEM_REDUCED_HCT_TRIANGLE';\n",  
		"TRI_FEM_TYPE = 'FEM_REDUCED_HCT_TRIANGLE';\n",  
		"TRI_DIRICHLET_FEM_TYPE = 'FEM_PK(2,2)';\n",  
		"TRI_DIRICHLET_DER_FEM_TYPE = 'FEM_PK(2,1)';\n",  
		"TRI_INTEGRATION = 'IM_HCT_COMPOSITE(IM_TRIANGLE(13))';\n",  
		"TRI_MORTAR_FEM_TYPE = 'FEM_PK(2,2)';\n",  
		"TRI_MORTAR_DERIV_FEM_TYPE = 'FEM_PK(2,1)';\n \n";
		print PARAM
		"QUAD_MESH_TYPE = 'GT_QK(2,1)'; \n",   
		"QUAD_DATA_FEM_TYPE = 'FEM_QK(2, 5)';\n",  
		"QUAD_PARTITION_OF_UNITY_FEM_TYPE = 'FEM_REDUCED_QUADC1_COMPOSITE';\n",  
		"QUAD_FEM_TYPE = 'FEM_REDUCED_QUADC1_COMPOSITE';\n",  
		"QUAD_DIRICHLET_FEM_TYPE = 'FEM_QK(2,2)';\n",  
		"QUAD_DIRICHLET_DER_FEM_TYPE = 'FEM_QK(2,1)';\n",  
		"QUAD_INTEGRATION = 'IM_QUADC1_COMPOSITE(IM_TRIANGLE(13))';\n",  
		"QUAD_MORTAR_FEM_TYPE ='FEM_QK(2,2)';\n",  
		"QUAD_MORTAR_DERIV_FEM_TYPE = 'FEM_QK(2,1)';\n \n";
		}
		# xfem parameters :
		print PARAM
		"%%%  XFEM PARAMETERS %%% \n",
		"ENRICHMENT_OPTION = ".$xfem_enrichment."; \n",		
		"RADIUS_ENR_AREA = ".$radius_enrich."; \n",
		"SING_BASE_TYPE = ".$sing_base_type." ; \n",
		"RESIDUAL = 1E-9 ; \n",
		"DIRICHLET_VERSION = 0; \n",
		"MORTAR_VERSION = 0 ; \n",
		"EPS_DIRICHLET_PENAL = 1E-12 ; \n",
		"NORM_EXACT = 0 ; \n",
		"RADIUS_SPLIT_DOMAIN = 0.0 ; \n",
		"ROOTFILENAME = 'serie_".$i."' ; \n",
		"VTK_EXPORT = 0; \n", 
		"MATLAB_EXPORT = 0; \n",
		"SIMPLEX_INTEGRATION = 'IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(13),3)';\n",
		"SINGULAR_INTEGRATION = 'IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2, 13), 9)';\n \n";
		
		# treating the mortar case as optionnal
		if ($xfem_enrichment == 3){
		   print PARAM
		   "%%% MORTAR PARAMETERS %%% \n",
		   "MULT_WITH_H = 1 ; \n",
		   "MORTAR_WITHOUT_SINGUL = 0 ; \n",
		   "MORTAR_TYPE = ".$mortar_type." ; \n",
		   "SEUIL = 1e-26 ; \n \n" ; 
		}
		
		
		close(PARAM) ;
		
		my @result = start_program("");
		print "le tableau :\n", @result, "\n" ;
		
		@cond[$k]   = @result[0] ;
		@err_L2[$k] = @result[1] ;
		@err_H1[$k] = @result[2] ;
		@err_H2[$k] = @result[3] ;
		@fic1[$k] = @result[4] ;
		@fic2[$k] = @result[5] ;
		@fic3[$k] = @result[6] ;
		@fic4[$k] = @result[7] ;
		$k += 1 ;
	#$i += 1 ;
}

# print "err_L2 =" ;
# foreach(@err_L2) {
#    print $_ ;}	
# print "\n" ;

# Ecriture des resultats dans un fichier .m
open(SORTIE, ">".$fichier_sortie);

print SORTIE "cond = [";
foreach (@cond){
   print SORTIE $_, " ";}
print SORTIE "]; \n" ;

print SORTIE "L2_error = [";
foreach (@err_L2){
   print SORTIE $_, " ";}
print SORTIE "]; \n" ;

print SORTIE "H1_error = [";
foreach (@err_H1){
   print SORTIE $_, " ";}
print SORTIE "]; \n" ;

print SORTIE "H2_error = [";
foreach (@err_H2){
   print SORTIE $_, " ";}   
print SORTIE "]; \n" ;

print SORTIE "FIC1 = [";
foreach (@fic1){
   print SORTIE $_, " ";}
print SORTIE "]; \n" ;

print SORTIE "FIC2 = [";
foreach (@fic2){
   print SORTIE $_, " ";}
print SORTIE "]; \n" ;

if ($sing_base_type eq 0){
   print SORTIE "FIC3 = [";
   foreach (@fic3){
      print SORTIE $_, " ";}
   print SORTIE "]; \n" ;

   print SORTIE "FIC4 = [";
   foreach (@fic4){
      print SORTIE $_, " ";}
   print SORTIE "]; \n" ;
}

# print informations relative to meshes :

my (@h, @nb_tri) ;
my $l = 0 ;

foreach (@mesh_names){
		@h[scalar(@h)] = 0. ;
}
if ($mesh_files_or_nx eq 1){
   @nb_tri = (312, 460, 758, 1498, 3086, 5182, 8500, 12254, 16746, 22222, 34006) ;
   foreach( @h) {
      @h[$l] = sqrt(2./@nb_tri[$l]) ;
      $l += 1 ;
   }
}
else{
   foreach( @h) {
      @h[$l] = 1./@mesh_names[$l] ;
      $l += 1 ;
   }
}

print SORTIE "h = [";
foreach (@h){
   print SORTIE $_, " ";}
print SORTIE "]; \n" ;








