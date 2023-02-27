PROFILE = -pg
PROFILE =
DEBUG = -g $(PROFILE)
DEBUG =
XLIB = 
XLIB = -lX11
XFILE = 
XFILE = ca-X.o
PLOTLIB = -l4014

ca:	ca-main.o ca-entropy.o ca-eval.o ca-parse.o ca-plot.o \
							ca-utilities.o $(XFILE)
	cc ca-main.o ca-entropy.o ca-eval.o ca-parse.o ca-utilities.o \
		$(XFILE) ca-plot.o $(DEBUG) $(XLIB) $(PLOTLIB) -lm \
		-o ca

ca-main.o:	ca-main.c ca.h
	cc -c  $(DEBUG) -O2 ca-main.c

ca-entropy.o:	ca-entropy.c ca.h
	cc -c  $(DEBUG) -O2 ca-entropy.c

ca-eval.o:	ca-eval.c ca.h
	cc -c  $(DEBUG) -O2 ca-eval.c

ca-parse.o:	ca-parse.c ca.h
	cc -c  $(DEBUG) -O2 ca-parse.c

ca-utilities.o:	ca-utilities.c ca.h
	cc -c  $(DEBUG) -O2 ca-utilities.c

ca-X.o:		ca-X.c ca.h ca_cursor.h ca_cursor_mask.h
	cc -c  $(DEBUG) -O2 ca-X.c

ca-plot.o:		ca-plot.c ca.h
	cc -c  $(DEBUG) -O2 ca-plot.c

# End of makefile
