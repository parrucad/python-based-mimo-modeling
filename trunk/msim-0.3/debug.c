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

/* $Id: debug.c,v 1.4 2004/07/30 23:57:29 miguel Exp $ */

/* MIMO Simulation */

/* Miguel Bazdresch */

/* Debug function.
These functions print to a file the variables or warnings specified by the variable cDebug */

void debug1(int iterations, double anorm) {
    /* warn no convergence found */
    fprintf(svdfile, "\nWarning:");
    fprintf(svdfile, "no convergence in %u dsvdcmp iterations\n", iterations);
    fprintf(svdfile, "anorm = %1.20f\n", anorm);
}

void debug2(int a1) {
    int Hrow, Hcol;
    if (a1 == 1) {
        /* print estimated H matrices*/
        fprintf(hfile, "%% Channel H Estimate\n");
        fprintf(hfile, "frame = %u\n", iFrameCounter);
        fprintf(hfile, "EstH = [");
        for(Hrow=1;Hrow<=N2;Hrow++) {
            fprintf(hfile, "[");
            for(Hcol=1;Hcol<=M2;Hcol++)
                fprintf(hfile, "%1.5f ", pH[Hrow][Hcol]);
                fprintf(hfile, "];\n      ");
        }
        fprintf(hfile, "];\n\n");
    }
    else {
        /* print H matrix to file */
        fprintf(hfile, "%%%% Channel H Matrix.\n");
        fprintf(hfile, "frame = %u\n", iFrameCounter);
        fprintf(hfile, "H = [");
        for(Hrow=1;Hrow<=N2;Hrow++) {
            fprintf(hfile, "[");
            for(Hcol=1;Hcol<=M2;Hcol++)
                fprintf(hfile, "%1.5f ", pH[Hrow][Hcol]);
                fprintf(hfile, "];\n      ");
        }
        fprintf(hfile, "];\n");
    }
}

void debug4(void) {
    int i, j;
    /* print received matrix to outfile */
    fprintf(outfile, "frame = %d\n", iFrameCounter);
    fprintf(outfile, "\n%%%% Received Matrix.\n");
    fprintf(outfile, "RxMx = [");
    for(i=1;i<=N;i++) {
        fprintf(outfile, "[");
        for(j=1;j<=(Lt+L);j++) {
            fprintf(outfile, "%1.5f + %1.5f*j, ", pReceived[i][j], pReceived[i+N][j]);
        }
        fprintf(outfile, "];\n        ");        
    }
    fprintf(outfile, "];\n");
}

void debug8(int a1, double noise, double MatMult) {
    int i, j;
    if (a1 == 1) {
        /* calculate sent power (modulated matrix) */
        for(i=1;i<=M;i++) {
            for(j=1;j<=(Lt+L)*iBitsPerSymbol;j++) {
                EstSentPwr += SQR(pModulated[i][j]);
                EstSentPwr += SQR(pModulated[i+M][j]);
            }
        }
    }
    else {
        /* calculate received and noise power */
        EstNoisePwr += SQR(noise);
        EstRxPwr += SQR(MatMult);
    }
}

void debug16(void) {
    int row, column;
    /* output coded matrix to out file */
    fprintf(outfile, "\n%%%% Transmit (coded) matrix.\n");
    fprintf(outfile, "frame = %d\n", iFrameCounter);
    fprintf(outfile, "TxMx = [");
    for(row=1;row<=M;row++) {
        fprintf(outfile, "[");
        for(column=1;column<=(Lt+L)*iBitsPerSymbol;column++)
            fprintf(outfile, "%u, ", pTransmitMatrix[row][column]);
            fprintf(outfile, "];\n       ");
    }
    fprintf(outfile, "];\n");
}

void debug32(void) {
    int row, col;
    /* print modulated matrix if debug is enabled*/
    fprintf(outfile, "\n%%%% Modulated Matrix.\n");
    fprintf(outfile, "MdMx = [");
    for(row=1;row<=M;row++){
        fprintf(outfile, "[");
        for(col=1;col<=(Lt+L);col++)
            fprintf(outfile, "%1.5f + %1.5f*j, ", pModulated[row][col], pModulated[row+M][col]);
            fprintf(outfile, "];\n        ");
    }
    fprintf(outfile, "];\n");
}

void debug64(int a1, int a2) {
    int index, Hrow, Hcol;
    /* print VBLAST variables */
    fprintf(outfile, "-----------------------\n");
    fprintf(outfile, "Frame = %d, column = %d\n\n", iFrameCounter, a1);
    /* print pH*/
    fprintf(outfile, "H = [ ");
    for(Hrow=1;Hrow<=N2;Hrow++) {
        fprintf(outfile, "[");
        for(Hcol=1;Hcol<=M2;Hcol++)
            fprintf(outfile, "%1.5f, ", pH[Hrow][Hcol]);
	fprintf(outfile, "];\n       ");
    }
    fprintf(outfile, "];\n");
    /* print pseudoinverse*/
    fprintf(outfile, "G = [ ");
    for(Hrow=1;Hrow<=M2;Hrow++) {
        fprintf(outfile, "[");
        for(Hcol=1;Hcol<=N2;Hcol++)
            fprintf(outfile, "%1.5f, ", storedpinvs[a2][Hrow][Hcol]);
            fprintf(outfile, "];\n       ");
    }
    fprintf(outfile, "];\n");
    /* print original r vector*/
    fprintf(outfile, "r_original = [");
    for(index=1;index<=N;index++) {
        fprintf(outfile, "[%1.3f + %1.3f*j]; ", pReceived[index][a1], pReceived[index+N][a1]);
    }
    fprintf(outfile, "];\n");
    /* print receiver variables */
    fprintf(outfile, "i = %u; k[%u] = %u\n", a2, a2, k[a2]);
    /* print received, modified vector */
    fprintf(outfile, "r_current = [");
    for(index=1;index<=N;index++) {
        fprintf(outfile, "[%1.3f + %1.3f*j]; ", r[index], r[index+N]);
    }
    fprintf(outfile, "];\n");
    /* print weight vector */
    fprintf(outfile, "w%u = [", a2);
    for(index=1;index<=N;index++) {
        fprintf(outfile, "%1.5f + %1.5f*j, ", wg[index], wg[index+N]);
    }
    fprintf(outfile, "];\n");
    /* print y */
    fprintf(outfile, "y[%u] = %1.5f + %1.5f*j\n", k[a2], yReal, yImag);
    /* print a */
    fprintf(outfile, "a[%u] = %1.5f + %1.5f*j\n", k[a2], a[k[a2]], a[k[a2]+M]);
}

void debug128(int a1, double * pW) {
    int col;
    if (a1 == 0) {
        /* print calculated singular values to file */
        fprintf(outfile, "singularvalue=[");
        for(col=1;col<=2*M;col++) {
            fprintf(outfile, "%1.5f, ", pW[col]);
        }
        fprintf(outfile, "];\n\n");
    }
    else {
        /* print inverted singular values to file */
        fprintf(outfile, "singinv=[");
        for(col=1;col<=2*M;col++) {
            fprintf(outfile, "%1.5f, ", pW[col]);
        }
        fprintf(outfile, "];\n\n");
    }
}

void debug256(int a1) {
    /* print differences between estimated and transmitted bits */
    fprintf(errfile, "estimated = %u, source = %u, index = %u\n", pEstimatedBits[a1], pSourceBits[a1], a1);
}

void debug1024(int a1, int column) 
{
    int row,col,i=0;
    
        /* print mlreceiver data to mimosim.out */
    if(a1==1) {
        fprintf(outfile, "G = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=N2;col++) {
                fprintf(outfile, "%1.5f, ", G[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
    else if(a1==2) {
        fprintf(outfile, "G_i = [");
        for(row=1;row<=N2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=M2;col++) {
                fprintf(outfile, "%1.5f, ", G_i[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
    else if(a1==3) {
        fprintf(outfile, "Gl = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=N2;col++) {
                fprintf(outfile, "%1.5f, ", G[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
    else if(a1==4) {
        fprintf(outfile, "G3_t = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=M2;col++) {
                fprintf(outfile, "%1.5f, ", G3_t[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
        fprintf(outfile, "Q_t = [");
        for(row=1;row<=N2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=M2;col++) {
                fprintf(outfile, "%1.5f, ", Q_t[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
    else if(a1==5) {
        fprintf(outfile, "H3 = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=M2;col++) {
                fprintf(outfile, "%1.5f, ", H3[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
    else if(a1==6) {
        fprintf(outfile, "frame = %d\ncolumn = %d\n", iFrameCounter+1, column-Lt);
        for(i=1;i<=M;i++)
            fprintf(outfile, "mod=(%1.5f + j*%1.5f)\n", pModulated[i][column], pModulated[i+M][column]);
        fprintf(outfile, "rec:\n");
        for(i=1;i<=N;i++)
            fprintf(outfile, "rec=(%1.5f + j*%1.5f)\n", pReceived[i][column], pReceived[i+N][column]);
        for(i=1;i<=N2;i++)
            fprintf(outfile, "translate=%1.5f, ", translate[i]);
        fprintf(outfile, "\n");
        for(i=1;i<=N;i++)
            fprintf(outfile, "rectrans=(%1.5f + j*%1.5f)\n", pReceived[i][column]+translate[i], pReceived[i+N][column]+translate[i+N]);
    }
    else if(a1==7) {
        for(i=1;i<=M;i++)
            fprintf(outfile, "x3=(%1.5f + j*%1.5f)\n", x3[i], x3[i+M]);
    }
    else if(a1==8) {
        for(i=1;i<=M;i++)
            fprintf(outfile, "u3=(%1.5f + j*%1.5f)\n", u3[i], u3[i+M]);
        for(i=1;i<=N;i++)
            fprintf(outfile, "foo=(%1.5f + j*%1.5f)\n", foo[i], foo[i+M]);
    }
    else if(a1==9)
        fprintf(outfile, "xest[%d]=%1.5f\n", i, xest[i]);
    else if(a1==10)
        fprintf(outfile, "ERROR# %d  BLOCKERROR# %d FRAME# %d\n", biterrors, blockerrors, iFrameCounter+1);
}

void debug2048(int a1) 
{
    int row,col;
    
    if(a1==1) {
            /* fprintf bstar M2xN2, mu M2xN2, B 1xM2, bt 1xN2 */
        fprintf(outfile, "\nbstar = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=N2;col++) {
                fprintf(outfile, "%1.5f, ", bstar[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
        fprintf(outfile, "\nmu = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=N2;col++) {
                fprintf(outfile, "%1.5f, ", mu[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
        fprintf(outfile, "\nB = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            fprintf(outfile, "%1.5f, ", B[row]);
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
        fprintf(outfile, "\nbt = [");
        for(row=1;row<=N2;row++) {
            fprintf(outfile, "[");
            fprintf(outfile, "%1.5f, ", bt[row]);
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
}

void debug4096(int i, int Mr)
{
    int row, col;
    
    if(i==0) {
            /* fprintf pZ N2xMr, pTT MrxMr */
        fprintf(outfile, "\npZ = [");
        for(row=1;row<=N2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=Mr;col++) {
                fprintf(outfile, "%1.5f, ", pZ[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
        fprintf(outfile, "\npTT = [");
        for(row=1;row<=Mr;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=Mr;col++) {
                fprintf(outfile, "%1.5f, ", pTT[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
    }
    if(i==1) {
            /* fprintf pTTinv MrxMr, and pT MrxN2 */
        fprintf(outfile, "\npTTinv = [");
        for(row=1;row<=Mr;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=Mr;col++) {
                fprintf(outfile, "%1.5f, ", pTTinv[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");
        fprintf(outfile, "\npT = [");
        for(row=1;row<=Mr;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=N2;col++) {
                fprintf(outfile, "%1.5f, ", pT[row][col]);
            }
            fprintf(outfile, "];\n");
        }
        fprintf(outfile, "];\n");        
    }
    if(i==2) {
            /* fprintf pHWork N2xM2 */
        fprintf(outfile, "\npHWork = [");
        for(row=1;row<=N2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=M2;col++) {
                fprintf(outfile, "%1.5f, ", pHWork[row][col]);
            }
            fprintf(outfile, "];\n");
        }
    }
    if(i==3) {
            /* fprintf psedoinverse M2xN2 */
        fprintf(outfile, "\nMPPI = [");
        for(row=1;row<=M2;row++) {
            fprintf(outfile, "[");
            for(col=1;col<=N2;col++) {
                fprintf(outfile, "%1.5f, ", storedpinvs[Mr][row][col]);
            }
            fprintf(outfile, "];\n");
        }
    }
}

void debugPrintMatrixR(int M, int N, double **H, char name[]) {
    int i,j;

    printf(name);
    for(i=1;i<=M;i++) {
	printf("[");
	for(j=1;j<=N;j++) {
	    printf("%1.5f ", H[i][j]);
	}
	printf("];\n");
    }
    printf("];\n");
}

void debugPrintMatrixC(int M, int N, double **H, char name[]) {
    int i,j;

    printf(name);
    for(i=1;i<=M;i++) {
	printf("[");
	for(j=1;j<=N;j++) {
	    printf("%1.5f + j*%1.5f, ", H[i][j], H[i+M][j]);
	}
	printf("];\n");
    }
    printf("];\n");
}

void debugPrintVectorI(int M, int *H, char name[]) {
    int i;

    printf(name);
    for(i=1;i<=M;i++) {
	printf("%d, ", H[i]);
    }
    printf("];\n");
}

void debugPrintVectorR(int M, double *H, char name[]) {
    int i;

    printf(name);
    for(i=1;i<=M;i++) {
	printf("%1.5f, ", H[i]);
    }
    printf("];\n");
}

void debugPrintVectorC(int M, double *H, char name[]) {
    int i;

    printf(name);
    for(i=1;i<=M;i++) {
	printf("%1.5f + j*%1.5f, ", H[i],H[i+M]);
    }
    printf("];\n");
}
