convexityCM10: convexityCM10.c gpc.o gpc.h
	cc -o convexityCM10 convexityCM10.c gpc.o -lm -O

skull3: skull3.c gpc.o gpc.h
	cc -o skull3 skull3.c gpc.o -lm -O

convex2: convex2.c
	cc     convex2.c   -o convex2 -lm

fit_triangle: fit_triangle.c
	cc -o fit_triangle fit_triangle.c -lm

pix_transf: pix_transf.c
	cc -o pix_transf pix_transf.c -lm

parallelogram: parallelogram.c
	cc -o parallelogram parallelogram.c -lm

NRCOBJ = powell.o nrutil.o linmin.o brent.o mnbrak.o f1dim.o
FITOBJ2 = superellipse_fit2.o ellipse_fit.o circle_fit.o moments.o ellipse_specific2.o

se2: se2.o $(FITOBJ2) $(NRCOBJ)
	cc -o se2 se2.o $(FITOBJ2) $(NRCOBJ) -lm
