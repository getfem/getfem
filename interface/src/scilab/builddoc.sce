// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM
// 
//  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.



Files = ['gf_asm', 'gf_global_function', ...
         'gf_mesh_get', 'gf_model_get', 'gf_spmat_get', ...
         'gf_compute', 'gf_global_function_get', ...
         'gf_mesh_im', 'gf_model_set', 'gf_spmat_set', ...
         'gf_cvstruct_get', 'gf_integ',  ...
         'gf_mesh_im_get', 'gf_poly', 'gf_undelete', ...
         'gf_delete', 'gf_integ_get', ...
         'gf_mesh_im_set', 'gf_precond', 'gf_util', ...
         'gf_eltm', 'gf_levelset', ...
         'gf_mesh_levelset', 'gf_precond_get', 'gf_workspace', ...
         'gf_fem', 'gf_levelset_get', 'gf_mesh', ...
         'gf_mesh_levelset_get', 'gf_slice', 'gf_fem_get', ...
         'gf_levelset_set', 'gf_mesh_fem', 'gf_mesh_levelset_set', ...
         'gf_slice_get', 'gf_geotrans', 'gf_linsolve', ...
         'gf_mesh_fem_get', 'gf_mesh_set', 'gf_slice_set', ...
         'gf_geotrans_get', 'gf_mesh_fem_set', ...
         'gf_model', 'gf_spmat'];

Files = ['gf_asm'];

Path = '../';

DocTokens = ['GFDOC','FUNC','MATLABEXT','MATLABFUNC','INIT',...
             'RDATTR','GET','SET','ARGS'];

TypesTable = ['@imat',   'imat'; ...
              '@ivec',   'ivec'; ...
              '@cvec',   'vec'; ...
              '@dcvec',  'vec'; ...
              '@dvec',   'vec'; ...
              '@vec',    'vec'; ...
              '@dmat',   'mat'; ...
              '@mat',    'mat'; ...
              '@str',    'string'; ...
              '@int',    'int'; ...
              '@bool',   'bool'; ...
              '@real',   'real'; ...
              '@scalar', 'scalar'; ...
              '@list',   'list'; ...
              '@tpoly',  'poly'; ...
              '@tmf',    'mesh_fem'; ...
              '@tgt',    'geotrans'; ...
              '@tgf',    'global_function'; ...
              '@tmls',   'mesh_levelset'; ...
              '@tmim',   'mesh_im'; ...
              '@tls',    'levelset'; ...
              '@tsl',    'slice'; ...
              '@tsp',    'spmat'; ...
              '@tpre',   'precond'; ...
              '@CELL',   ''];

for i=1:size(Files,'*')
  printf('Processing %s\n', Files(i));

  DocPage = mlist(['gf','gfdoc','func','matlabext','matlabfunc','init','set', 'get', 'args','rdattr'], ...
                        list([]), list([]),list([]),     list([]),      list([]),list([]),list([]),list([]),list([]));

  fid   = mopen(Path + Files(i) + '.cc', 'r');
  Lines = mgetl(fid, -1);
  mclose(fid);

  // Put block of comments with title DocTokens(j) in a mlist
  for j=1:size(DocTokens,'*')
    Index = grep(Lines,'/*@'+DocTokens(j));
    if ~isempty(Index) then
      for k=1:size(Index,'*')
        IndPos = Index(k);
        Result = isempty(strindex(Lines(IndPos),'@*/'));
        DocPage(convstr(DocTokens(j),'l'))(k) = [];
        while Result    
          DocPage(convstr(DocTokens(j),'l'))(k) = [DocPage(convstr(DocTokens(j),'l'))(k); Lines(IndPos)];
          Result = isempty(strindex(Lines(IndPos),'@*/'));
          IndPos = IndPos + 1;
        end
      end
    end
  end

  // Now, remove parsed tokens and types substitution
  for j=1:size(DocTokens,'*')
    // Tokens removal 
    for k=1:length(DocPage(convstr(DocTokens(j),'l')))
      // Starting token removal
      Index = grep(DocPage(convstr(DocTokens(j),'l'))(k), '/*@'+DocTokens(j));
      if ~isempty(Index) then
        for l=1:size(Index,'*')
          DocPage(convstr(DocTokens(j),'l'))(k)(Index(l)) = strsubst(DocPage(convstr(DocTokens(j),'l'))(k)(Index(l)),'/*@'+DocTokens(j),'');
        end
      end

      // Ending token removal
      Index = grep(DocPage(convstr(DocTokens(j),'l'))(k), '@*/');
      if ~isempty(Index) then
        for l=1:size(Index,'*')
          DocPage(convstr(DocTokens(j),'l'))(k)(Index(l)) = strsubst(DocPage(convstr(DocTokens(j),'l'))(k)(Index(l)),'@*/','');
        end
      end
    end

    // Types substitution
    for k=1:length(DocPage(convstr(DocTokens(j),'l')))
      for l=1:size(TypesTable,1)
        for m=1:size(DocPage(convstr(DocTokens(j),'l'))(k),'*')
          DocPage(convstr(DocTokens(j),'l'))(k)(m) = strsubst(DocPage(convstr(DocTokens(j),'l'))(k)(m),TypesTable(l,1),TypesTable(l,2));
        end
      end
    end

    // Strip blanks
    for k=1:length(DocPage(convstr(DocTokens(j),'l')))
      for l=1:size(DocPage(convstr(DocTokens(j),'l'))(k),'*')
        DocPage(convstr(DocTokens(j),'l'))(k)(l) = stripblanks(DocPage(convstr(DocTokens(j),'l'))(k)(l));
      end
    end

    // Remove empty lines
    for k=1:length(DocPage(convstr(DocTokens(j),'l')))
      tmp_str = [];
      for l=1:size(DocPage(convstr(DocTokens(j),'l'))(k),'*')
        if ~isempty(DocPage(convstr(DocTokens(j),'l'))(k)(l)) then
          tmp_str = [tmp_str; DocPage(convstr(DocTokens(j),'l'))(k)(l)];
        end
      end
      DocPage(convstr(DocTokens(j),'l'))(k) = tmp_str;
    end

    // Process math expressions
    for k=1:length(DocPage(convstr(DocTokens(j),'l')))
      for l=1:size(TypesTable,1)
        for m=1:size(DocPage(convstr(DocTokens(j),'l'))(k),'*')
          DocPage(convstr(DocTokens(j),'l'))(k)(m) = strsubst(DocPage(convstr(DocTokens(j),'l'))(k)(m),':math:`','<latex style=""text"">');
          DocPage(convstr(DocTokens(j),'l'))(k)(m) = strsubst(DocPage(convstr(DocTokens(j),'l'))(k)(m),'`','</latex>');
        end
      end
    end
  end

  // XML Processing
  fid = mopen('help/tmp/'+Files(i)+'.xml','w');
 
  // Header
  mfprintf(fid,"<?xml version=""1.0"" encoding=""UTF-8""?>\n");
  mfprintf(fid,"<refentry version=""5.0-subset Scilab"" xml:id=""%s"" xml:lang=""en""\n", Files(i));
  mfprintf(fid,"          xmlns=""http://docbook.org/ns/docbook""\n");
  mfprintf(fid,"          xmlns:xlink=""http://www.w3.org/1999/xlink""\n");
  mfprintf(fid,"          xmlns:xi=""http://www.w3.org/2001/XInclude""\n");
  mfprintf(fid,"          xmlns:svg=""http://www.w3.org/2000/svg""\n");
  mfprintf(fid,"          xmlns:mml=""http://www.w3.org/1998/Math/MathML""\n");
  mfprintf(fid,"          xmlns:html=""http://www.w3.org/1999/xhtml""\n");
  mfprintf(fid,"          xmlns:db=""http://docbook.org/ns/docbook"">\n");

  // Refnamediv
  mfprintf(fid,"<refnamediv>\n");
  mfprintf(fid,"  <refname>%s</refname>\n", Files(i));
  mfprintf(fid,"  <refpurpose>%s</refpurpose>\n",DocPage('gfdoc')(1)(1));
  mfprintf(fid,"</refnamediv>\n");

  // Synopsis
  mfprintf(fid,"<refsynopsisdiv>\n");
  mfprintf(fid,"  <title>Calling Sequence</title>\n");
  mfprintf(fid,"  <synopsis>\n");
  for j=1:length(DocPage('func'))
    mfprintf(fid, "%s\n", DocPage('func')(j)(1));
  end
  mfprintf(fid,"  </synopsis>\n");
  mfprintf(fid,"</refsynopsisdiv>\n");

  // Description
  mfprintf(fid,"<refsection>\n");
  mfprintf(fid,"  <title>Description</title>\n");
  mfprintf(fid,"  <para>\n");
  for j=1:size(DocPage('gfdoc')(1),'*')
    mfprintf(fid,"  %s\n", DocPage('gfdoc')(1)(j));
  end
  mfprintf(fid,"  </para>\n");

  mfprintf(fid,"  <itemizedlist>\n");
  for j=1:length(DocPage('func'))
    mfprintf(fid,"    <listitem>\n");
    mfprintf(fid,"      <para>\n");
    for k=1:size(DocPage('func')(j),"*")
      mfprintf(fid,"        %s\n", DocPage('func')(j)(k));
    end  
    mfprintf(fid,"      </para>\n");
    mfprintf(fid,"    </listitem>\n");
  end
  mfprintf(fid,"  </itemizedlist>\n");
  mfprintf(fid, "</refsection>\n");

  // See also section
  mfprintf(fid, "<refsection>\n");
  mfprintf(fid, "  <title>See Also</title>\n");
  mfprintf(fid, "    <simplelist type=""inline"">\n");
  mfprintf(fid, "    <member><link linkend=""gf_solve"">gf_solve</link></member>\n");
  mfprintf(fid, "  </simplelist>\n");
  mfprintf(fid, "</refsection>\n");

  // Author section
  mfprintf(fid, "<refsection>\n");
  mfprintf(fid, "  <title>Authors</title>\n");
  mfprintf(fid, "  <para>Y. Collette</para>\n");
  mfprintf(fid, "</refsection>\n");

  // Close the xml document
  mfprintf(fid, "</refentry>\n");

  mclose(fid);
end // End for
