INCLUDE=-I/data/local/include -I${HOME}/primesieve-5.4.2/include
LIBRARY=-L/data/local/lib -L${HOME}/primesieve-5.4.2/.libs
INCLUDE_G2=-I${HOME}/local/include -I/data/local/include -I${HOME}/primesieve-5.4.2/include
LIBRARY_G2=-L${HOME}/local/lib -L/data/local/lib -L${HOME}/primesieve-5.4.2/.libs
#HWLPOLY=
HWLPOLY=-lhwlpoly

all: glfunc glfunc_g2 smalljac_test

glfunc:	glfunc.cpp main.cpp defines.h utils.h structs.h error.h acb_fft.h coeffs.h artin.h g2_g3_ell.h modform.h turing.h upsample.h zeros.h div_zeros.h output.h ell_qd.h ../andy_elliptic/pari.c
	echo "Building with includes/libs from /data/local"
	g++ -O2 glfunc.cpp  -o glfunc -finline-functions -fomit-frame-pointer -fopenmp ${INCLUDE} ${LIBRARY} -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari -lprimesieve ${HWLPOLY} -lmpfr

glfunc_g2:	glfunc.cpp main.cpp defines.h utils.h structs.h error.h acb_fft.h coeffs.h artin.h g2_g3_ell.h modform.h turing.h upsample.h zeros.h div_zeros.h output.h ../andy_elliptic/pari.c
	echo "Building with includes/libs from /${HOME}/local"
	g++ -O2 glfunc.cpp  -o glfunc_g2 -finline-functions -fomit-frame-pointer -fopenmp ${INCLUDE_G2} ${LIBRARY_G2} -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari -lprimesieve ${HWLPOLY} -lmpfr

smalljac_test: smalljac_test.cpp
	g++ -O2 smalljac_test.cpp  -o smalljac_test -finline-functions -fomit-frame-pointer -fopenmp ${INCLUDE} ${LIBRARY} -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari ${HWLPOLY}


clean: 
	rm -f glfunc
	rm -f glfunc_g2
	rm -f smalljac_test

