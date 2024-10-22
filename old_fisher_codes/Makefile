test:
	gcc test.c ps.c utils.c nr_utils.c fom.c -lm -Wall -O3\
          -DSWG36 -DGROWTH_GAMMA -DPLANCK_PAR -DGAUSS_CONVERGENCE
	gcc test-b.c ps.c utils.c nr_utils.c fom.c -lm -Wall -o test-b.x -O3\
          -DSWG36 -DGROWTH_GAMMA -DPLANCK_PAR -DGAUSS_CONVERGENCE -DRECONS_05
	gcc postprocess.c -lm -Wall -O3 -o getfom5.x -DWBIN_PRIOR=100.0 -DMAX_EIGEN_PRINT=10.0
#          -DSWG36 -DGROWTH_GAMMA -DGAUSS_CONVERGENCE -DPLANCK_PAR

