readme.txt of the scilab getfem interface

To compile this interface, you will need to work with the scilab-master (the version in the git repository) because this interface uses the new interface api.

So, to compile this interface:
- launch scilab
- go to the getfem++/interface/src/scilab directory
- exec builder.sce; // The documentation is not yet built
- exec loader.sce;
You can try to launch a demo (be careful, there is a lot of work needed before the interface can be really useable).
- cd demos
- exec testyc.sce;

TODO:
- verifications of the data send / get to / from the getfem interface. Specially the hybrid types (list of list of doubles + strings ...)
- build an interface to meshash so as to access a sparse incomplete lu / cholesky factorisation (needed by gf_solve).
- test the scripts (specially the graphic ones).
- automatic generation of the documentation from the good interface functions ...

