/* Copyright (c) 2004 Miguel Bazdresch

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Please see the file LICENSE for more details. */

/* $Id: source.c,v 1.2 2004/07/31 00:04:25 miguel Exp $ */

/* Miguel Bazdresch */

/* MIMO Simulation */

/* Source function
   This function generates a stream of random, uniformly
   distributed bits.
   The function returns a 0 if everything went alright, and
   returns different than zero if there was an error.
   For random numbers, this function uses ran2 from NR (p. 282) */

int source (void) {

    /* function prototypes */
    double nr_ran2(long *idum);

    /* local variables */
    unsigned int j;

    switch(cSourceType) {
        case 0: /* All zeroes */
            for(j=1;j<=iNumberInfoBitsFrame;j++) {
                pSourceBits[j] = 0;
            }
            break;
        case 1: /* All ones */
            for(j=1;j<=iNumberInfoBitsFrame;j++) {
                pSourceBits[j] = 1;
            }
            break;
        case 2: /* random */
            for(j=1;j<=iNumberInfoBitsFrame;j++) {
                if(debug && ((cDebug & 512) == 512)) { /* read data from datafile */
                    int foo;
                    fscanf(datafile, "%d", &foo);
                    pSourceBits[j] = foo;
                }
                else {
                    pSourceBits[j] = gsl_rng_uniform(source_rng) > 0.5 ? 1 : 0;
                }
            }
            break;
	default:
	    return 1;
    }
    return 0;
}
