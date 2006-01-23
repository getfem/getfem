




all : simplexification_refelt




simplexification_refelt.o : simplexification_refelt.cc
	g++ -O3 -I ../src -I ../gcc/src simplexification_refelt.cc -c

simplexification_refelt : simplexification_refelt.o
	g++ simplexification_refelt.o -o simplexification_refelt ../gcc/src/.libs/libgetfem.a -lqhull -lblas -lg2c

