# this Makefile produces the numerical library jbnumlib
# and the library libchiron.a
# it also allows you to produce a number of testing/example programs
# it is part of the CHIRON ChPT at two loops program collection
#
# Copyright (C) 2014-2015 Johan Bijnens, v1.02
# CHIRON is licenced under the GNU GPL version 2 or later,
# see COPYING for details.
# Please respect the Guidelines, see GUIDELINES for details.

# -std=c++0x or -std=c++11 might be needed
# OPTIONS:
# in both testing and library
CXX = g++-10 -O3  -I./include
#CXX = g++ -O3  -I./include -std=c++11
# added for the compilation of the libraries
CFLAGS = -Wall -Wextra -Wconversion
#added for the testing files
CFLAGSTEST =  -Wall -Wextra -Wconversion

OBJECTSJBNUMLIB = jbdgauss.o jbdcauch.o jbwgauss.o \
  jbdgauss2.o jbdquad15.o jbdquad21.o \
  jbdcauch2.o jbdsing15.o jbdsing21.o \
  jbwgauss2.o jbwquad15.o jbwquad21.o \
  jbdadmul.o jbdadmul2.o \
  jbdli2.o jbdlin.o \
  jbdbesik.o jbdtheta30.o jbdtheta30m1.o jbdtheta32.o jbdtheta34.o \
  jbdtheta2d0.o jbdtheta2d0m1.o jbdtheta2d02.o \
  jbdtheta3.o jbderivutheta3.o jbderiv2utheta3.o jbderiv3utheta3.o \
  jbdrteq3.o jbdzerox.o jbdYlm.o jbdadmulc2.o

OBJECTSCHIRON = inputs.o Li.o Ci.o inputsnf.o Linf.o Ki.o inputsnf2.o \
       linf2.o \
       massdecayvev.o getfpimeta.o \
       massdecayvevnf2.o \
       massdecayvevlo.o massdecayvevPQ.o \
       massdecayvevnf.o massdecayvevnfPQ.o \
       oneloopintegrals.o sunsetintegrals.o finitevolumeoneloopintegrals.o \
       finitevolumesunsetintegrals.o quenchedsunsetintegrals.o \
       fourvector.o finitevolumeonelooptwist.o \
       vectorformlo.o quenchedoneloopintegrals.o vectorformPQ.o vectorformPQS.o\
       massdecayvevTV.o

OBJECTSCHIRONTHETA =  massdecayvevVt.o massdecayvevloVt.o massdecayvevPQVt.o \
       massdecayvevnfVt.o massdecayvevnfPQVt.o massdecayvevnf2Vt.o
OBJECTSCHIRONBESSEL =  massdecayvevVb.o massdecayvevloVb.o massdecayvevPQVb.o \
       massdecayvevnfVb.o massdecayvevnfPQVb.o massdecayvevnf2Vb.o
OBJECTSCHIRON2 = $(OBJECTSCHIRONTHETA) $(OBJECTSCHIRONBESSEL)

INCLUDECHIRON =  include/inputs.h include/Li.h include/Ci.h\
       include/inputsnf.h include/Linf.h include/Ki.h \
       include/inputsnf2.h include/linf2.h \
       include/massdecayvev.h \
       include/massdecayvevnf2.h \
       include/massdecayvevlo.h \
       include/massdecayvevPQ.h\
       include/massdecayvevnf.h \
       include/massdecayvevnfPQ.h \
       include/oneloopintegrals.h \
       include/sunsetintegrals.h\
       include/finitevolumeoneloopintegrals.h\
       include/finitevolumesunsetintegrals.h\
       include/massdecayvevV.h \
       include/massdecayvevnf2V.h \
       include/massdecayvevloV.h \
       include/massdecayvevPQV.h \
       include/quenchedsunsetintegrals.h \
       include/fourvector.h \
       include/finitevolumeonelooptwist.h \
       include/vectorformlo.h \
       include/quenchedoneloopintegrals.h \
       include/vectorformPQ.h \
       include/vectorformPQS.h \
       include/massdecayvevTV.h \
       include/getfpimeta.h

TESTCHIRON = testintegralsreal testintegralsrealsingular \
       testintegralscomplex \
       testinputs testinputsnf testLi testCi testmassdecayvev \
       testinputsnf2 testlinf2 \
       testgetfpimeta testoneloopintegrals testsunsetintegrals \
       testfinitevolumeoneloopintegrals testfinitevolumesunsetintegrals \
       testmassdecayvevV testquenchedsunsetintegrals testLinf testKi \
       testmassdecayvevlo testmassdecayvevPQ testmassdecayvevPQV \
       testmassdecayvevloV testmassdecayvevnf testfourvector \
       testfinitevolumeonelooptwist testjbdrteq3 testvectorformlo \
       testquenchedoneloopintegrals testvectorformPQ testvectorformPQS \
       testmassdecayvevTV testmassdecayvevnf2 testmassdecayvevnf2V \
       testjbdzerox
all: libchiron.a libjbnumlib.a

libchiron.a: $(OBJECTSCHIRON) $(OBJECTSCHIRON2)
	ar r libchiron.a $(OBJECTSCHIRON) $(OBJECTSCHIRON2) $(INCLUDECHIRON) \
        ; cp libchiron.a lib

libjbnumlib.a: $(OBJECTSJBNUMLIB) include/jbnumlib.h
	ar r libjbnumlib.a $(OBJECTSJBNUMLIB) include/jbnumlib.h ; cp libjbnumlib.a lib

massdecayvevVt.o: src/massdecayvevV.cc include/massdecayvevV.h
	$(CXX) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevVt.o src/massdecayvevV.cc
massdecayvevloVt.o: src/massdecayvevloV.cc include/massdecayvevloV.h
	$(CXX) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevloVt.o src/massdecayvevloV.cc
massdecayvevPQVt.o: src/massdecayvevPQV.cc include/massdecayvevPQV.h
	$(CXX) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevPQVt.o src/massdecayvevPQV.cc

massdecayvevnfVt.o: src/massdecayvevnfV.cc include/massdecayvevnfV.h
	$(CXX) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevnfVt.o src/massdecayvevnfV.cc

massdecayvevnfPQVt.o: src/massdecayvevnfPQV.cc include/massdecayvevnfPQV.h
	$(CXX) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevnfPQVt.o src/massdecayvevnfPQV.cc

massdecayvevVb.o: src/massdecayvevV.cc include/massdecayvevV.h
	$(CXX) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevVb.o src/massdecayvevV.cc

massdecayvevloVb.o: src/massdecayvevloV.cc include/massdecayvevloV.h
	$(CXX) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevloVb.o src/massdecayvevloV.cc

massdecayvevPQVb.o: src/massdecayvevPQV.cc include/massdecayvevPQV.h
	$(CXX) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevPQVb.o src/massdecayvevPQV.cc

massdecayvevnfVb.o: src/massdecayvevnfV.cc include/massdecayvevnfV.h
	$(CXX) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevnfVb.o src/massdecayvevnfV.cc

massdecayvevnfPQVb.o: src/massdecayvevnfPQV.cc include/massdecayvevnfPQV.h
	$(CXX) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevnfPQVb.o src/massdecayvevnfPQV.cc

massdecayvevnf2Vt.o: src/massdecayvevnf2V.cc include/massdecayvevnf2V.h
	$(CXX) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevnf2Vt.o src/massdecayvevnf2V.cc

massdecayvevnf2Vb.o: src/massdecayvevnf2V.cc include/massdecayvevnf2V.h
	$(CXX) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevnf2Vb.o src/massdecayvevnf2V.cc

$(OBJECTSJBNUMLIB): %.o: src/%.cc
	$(CXX) -c $(CFLAGS) $< -o $@

$(OBJECTSCHIRON): %.o: src/%.cc include/%.h
	$(CXX) -c $(CFLAGS) $< -o $@

# making testing programs, output is always a.out

$(TESTCHIRON): %: test/%.cc  libjbnumlib.a libchiron.a
	$(CXX) $(CFLAGSTEST)  -o a.out $<  -lchiron -ljbnumlib -L./lib


# some installing and cleaning up bits
.PHONY: clean installjbnumlib installchiron doc

install: libjbnumlib.a libchiron.a
	cp libjbnumlib.a libchiron.a ~/lib ; cp include/jbnumlib.h $(INCLUDECHIRON) ~/include

installchiron: libjbnumlib.a libchiron.a
	cp libjbnumlib.a libchiron.a ~/lib ; cp include/jbnumlib.h $(INCLUDECHIRON) ~/include

installjbnumlib: libjbnumlib.a
	cp libjbnumlib.a ~/lib ; cp include/jbnumlib.h ~/include

clean:
	rm *.o *.a a.out lib/*.a \
          temp.dat test.dat *.log *.out *.toc *.aux \
          *.idx *.ilg *.ind *.pdf

doc:
	pdflatex doc/manual.tex ;pdflatex doc/manual.tex ; makeindex manual ;\
          pdflatex doc/manual.tex ; makeindex manual ; \
          pdflatex doc/manual.tex ; mv manual.pdf doc
