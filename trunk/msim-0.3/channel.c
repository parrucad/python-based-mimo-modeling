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

/* $Id: channel.c,v 1.2 2004/07/30 23:57:21 miguel Exp $ */

/* Miguel Bazdresch

   MIMO simulator

   Channel function

   This function calculates vector rl = Ha + v, where
     rl is the received n-vector.
     H is a complex nxm matrix. hij is the transfer function
       from transmitter j to receiver i. m<=n. m is the number of
       transmitters, and n is the number of receivers.
       a is the transmitted m-vector. a is taken from the rows
       of the transmitted matrix, pModulatedReal + j*pModulatedImag.
     v is an m-vector with gaussian components. Its amplitude
       depends on noise power, fNoisePower.
   This function is repeated for each row in the transmitted
   matrix, and the results are stored in matrix
   pReceivedReal + j*pReceivedImag. The size of these matrices is
   the same as for the transmitted matrices, that is,
   cTransmitAntennas rows and iTxColumns columns.
   The matrix H is different for each frame. */

int channel(void) {
    
    /* function prototypes */
    double nr_gasdevnoise(long *idum);
    double nr_gasdevHch(long *idum);
    void debug4(void);
    void debug8(int a1, double d1, double d2);
    
    /* local variables */
    double MatMult, a, noise;
    unsigned int Hcol, Hrow, column, vector, i, j;
    
    column = 1;
    /* calculate H
       H is such that r = aH + v
       H dimensions are 2M*2N
       How to calculate H depends on the value of global variable cMatrixHComponentGen:
         0: random values
    */
    switch (cMatrixHComponentGen) {
        case 0:  /* random values */
            for(Hrow=1;Hrow<=N;Hrow++) {
                for(Hcol=1;Hcol<=M;Hcol++) {
                    /* Hvar determines the variance of the elements of H */
                    if (debug && ((cDebug & 512) == 512)) { /* read from file */
                        double foo;
                        fscanf(noisefile, "%lf", &foo);
                        a = Hvar*foo;
                    }
                    else {
                        a = gsl_ran_gaussian(hchannel_rng, Hvar);
                    }
                    pH[Hrow][Hcol] = a;
                    pH[Hrow+N][Hcol+M] = a;
                }
                for (Hcol=1;Hcol<=M;Hcol++) {
                    if (debug && ((cDebug & 512) == 512)) { /* read from file */
                        double foo;
                        fscanf(noisefile, "%lf", &foo);
                        a = Hvar*foo;
                    }
                    else {
                        a = gsl_ran_gaussian(hchannel_rng, Hvar);
                    }
                    pH[Hrow][Hcol+M] = -a;
                    pH[Hrow+N][Hcol] = a;
                }
            }
            break;
        case 1: /* identity matrix - channel doesn't alter signal - only for square systems */
            for (Hrow=1;Hrow<=N2;Hrow++) {
                for (Hcol=1;Hcol<=N2;Hcol++) {
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
    } /* end H matrix generation switch */
    /* Now multiply each vector in the frame by H and add noise */
    for (vector=1;vector<=(Lt+L);vector++) {
        for (i=1;i<=N2;i++) {
	    MatMult = 0.0;
	    for (j=1;j<=M2;j++) {
		MatMult += pModulated[j][column]*pH[i][j];
	    }
	    if( Nvar != 0 ) {
		if (debug && ((cDebug & 512) == 512)) { /* read from file */
		    double foo;
		    fscanf(noisefile, "%lf", &foo);
		    noise = foo*Nvar;
		}
		else {
		    noise = gsl_ran_gaussian(noise_rng,Nvar);
		}
	    }
	    else noise = 0;
	    pReceived[i][column] = MatMult + noise;
	    noisematrix[i][column] = noise;
	    if (debug && ((cDebug & 8) == 8))
		debug8(0, noise, MatMult); /* estimates of noise and received power */
	}
        column++;
    }
    if (debug && ((cDebug & 4) == 4)) debug4(); /* print received matrix to outfile */
    return 0;
}
