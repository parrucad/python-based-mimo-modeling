# $Id: Makefile,v 1.5 2004/08/16 17:33:51 miguel Exp $

mimosim : mimosim.c source.c parser.c channel.c chestimator.c coder.c debug.c matrix.c memoryalloc.c mlreceiver.c modulator.c nr_functions.c support.c support_ml.c vblastreceiver.c vb-reduced.c vb-wubben.c
	gcc -O2 -finline-functions -lgsl -lgslcblas -lm -Wall -o mimosim mimosim.c
	mkdir -p results

debug :  mimosim.c source.c parser.c channel.c chestimator.c coder.c debug.c matrix.c memoryalloc.c mlreceiver.c modulator.c nr_functions.c support.c support_ml.c vblastreceiver.c nr_licensed.c vb-reduced.c vb-wubben.c
	gcc -g -lm -lgsl -lgslcblas -Wall -o mimosim mimosim.c
	mkdir -p results
