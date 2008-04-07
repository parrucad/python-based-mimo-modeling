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

/* $Id: chestimator.c,v 1.1 2004/05/02 22:06:13 miguel Exp $ */

/* MIMO simulation
   Miguel Bazdresch

   Function chestimator() performs channel estimation. Estimated
   channel matrix is stored in pH */

int chestimator (void) {

    /* prototypes */
    void debug2(int a1);

    /* local variables */
    int Hrow, Hcol;
    double a;  

    if (debug && ((cDebug & 2) == 2)) debug2(0); /* before estimation, print H matrix to file */
        /* Estimate H. How to do this depends on variable cMethodChannelEstimation
           The estimate of H is saved in the same variables than the original H */
    if (cMethodChannelEstimation == 0) { /* use genie knowledge */
            /* if the program is reading H values from a file (debug
               == 512) then the H values must not be regenerated here,
               since it would require re-reading data from the
               noise.out file. */
        switch (cMatrixHComponentGen) {
            case 0: /* random values */
                if(!debug) {
                    for(Hrow=1;Hrow<=N;Hrow++) {
                        for(Hcol=1;Hcol<=M;Hcol++) {
                            /* Hvar determines the variance of the elements of H */
                            a = gsl_ran_gaussian(hrx_rng,Hvar);
                            pH[Hrow][Hcol] = a;
                            pH[Hrow+N][Hcol+M] = a;
                        }
                        for (Hcol=1;Hcol<=M;Hcol++) {
                            a = gsl_ran_gaussian(hrx_rng,Hvar);
                            pH[Hrow][Hcol+M] = -a;
                            pH[Hrow+N][Hcol] = a;
                        }
                    }
                } else {
                    if ((cDebug & 512) != 512) { /* NOT reading H values from file */
                        for(Hrow=1;Hrow<=N;Hrow++) {
                            for(Hcol=1;Hcol<=M;Hcol++) {
                                /* Hvar determines the variance of the elements of H */
                                a = gsl_ran_gaussian(hrx_rng,Hvar);
                                pH[Hrow][Hcol] = a;
                                pH[Hrow+N][Hcol+M] = a;
                            }
                            for (Hcol=1;Hcol<=M;Hcol++) {
                                a = gsl_ran_gaussian(hrx_rng,Hvar);
                                pH[Hrow][Hcol+M] = -a;
                                pH[Hrow+N][Hcol] = a;
                            }
                        }
                    }
                }
                break;
            case 1: /* H is identity matrix - channel has no effect - only for square systems */
                for (Hrow=1;Hrow<=M2;Hrow++) {
                    for (Hcol=1;Hcol<=M2;Hcol++) {
                        if (Hrow==Hcol) {
                            pH[Hrow][Hrow] = 1.0;
                        }
                        else {
                            pH[Hrow][Hcol] = 0.0;
                        }
                    }
                }
                break;
	    default:
		return 1;
        } /* end switch for H estimation */
    }
    if (debug && ((cDebug & 2) == 2)) debug2(1); /* print estimated H matrices to file */
    
    return 0;
}
