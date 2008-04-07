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

/* $Id: vb-wubben.c,v 1.4 2004/08/16 17:31:59 miguel Exp $ */

/* Wubben's sorted QR decomposition
 * H is of size NxM, where N >= M
 */
void ctQR_sorted( int N, int M, double **H, double **Q, double **R, unsigned int *perm ) {

    int    i, j, l, sc, temp;
    double smallest, mag;

    matrix_zero( R, M, M );

    xx {xh[mem_mz]+=2*M*M;} 

    matrix_copy( Q, H, N, M );

    xx {xh[mem_mc]+=4*N*M;} 

    for( i=1; i<=M; i++ )
	perm[i] = i;

    xx {xh[mem]+=M;}

    for( i=1; i<=M; i++ ) {
	
	smallest = 1e10;
	sc = 0;
	for( j=i; j<=M; j++ ) {
	    mag = col_magnitude( Q, N, j );
	    if( mag < smallest ) {
		smallest = mag;
		sc = j;
	    }
	    
	    xx {xh[mem]+=2*N;xh[mul]+=2*N;xh[add]+=2*N+1;}
	
	}

	col_exchange( Q, N, i, sc );
	col_exchange( R, M, i, sc );
	temp = perm[i];
	perm[i] = perm[sc];
	perm[sc] = temp;

	xx {xh[mem]+=8*N+8*M+4;}

	R[i][i] = sqrt( col_magnitude( Q, N, i ) );

	xx {xh[mem]+=2*N+1;xh[srt]++;xh[mul]+=2*N;xh[add]+=2*N;}

	col_divide( Q, N, i, R[i][i], 0 );

	xx {xh[mem]+=4*N;xh[mul]+=4*N+2;xh[add]+=2*N+1;xh[div]+=2*N;}

	for( j=i+1; j<=M; j++ ) {
	    R[i][j] = col_col_mult_real( Q, N,  i, j );
	    R[i+M][j] = col_col_mult_imag( Q, N, i, j);

	    xx {xh[mem]+=8*N+2;xh[mul]+=4*N;xh[add]+=4*N;}
	    
	    for( l=1; l<=N; l++ ) {
		Q[l][j] -= R[i][j] * Q[l][i] - R[i+M][j] * Q[l+N][i];
		Q[l+N][j] -= R[i+M][j] * Q[l][i] + R[i][j] * Q[l+N][i];
	    }
	    
	    xx {xh[mem]+=6*N;xh[add]+=4*N;xh[mul]+=4*N;}
	}
    }
}

int vb_wubben( void ) {

    /* function prototypes */
    void debug1024( int, int );
    void ctQR_sorted( int, int, double **, double **, double **, unsigned int * );
    void demodulate( double, double, int );

    /* local variables */
    int    column, i, j;
    int    infopointer;
    int    BLERflag, symflag, symindex;
    double tempr, tempi;

    /* infopointer is used to rebuild the data bit stream, once symbols
     * have been estimated. */
    infopointer = 1;

    /* NEEDFIX: 4-QAM is not tested - need several updates */

    /*
     * The first section of the receiver is done only once per
     * block. 
     *
     */
    /* debugPrintMatrixC( N, M, pH, "H = ["); */
    ctQR_sorted( N, M, pH, Q, G3, permut ); /* vector tempV keeps permutation order */
    /* debugPrintMatrixC( N, M, Q, "Q = ["); */
    /* debugPrintMatrixC( M, M, G3, "G3 = ["); */
    /* 
     * Here ends the first section. We now have all the data we need.
     * The next section is repeated for each received vector.
     * 
     */

    for ( column=(Lt+1); column<=(Lt+L); column++ ) {

	if (debug && ((cDebug & 1024) == 1024)) debug1024(6,column);
	
	/* find y = Q^H * x   (Wubben nomenclature)
	 *      x3 = Q^H * pReceived
	 */

	for( i=1; i<=M; i++ ) {
	    tempr = 0.0;
	    tempi = 0.0;
	    for( j=1; j<=N; j++ ) {
		tempr += Q[j][i] * pReceived[j][column] + \
		         Q[j+N][i] * pReceived[j+N][column];
		tempi += Q[j][i] * pReceived[j+N][column] - \
		         Q[j+N][i] * pReceived[j][column];
	    }
	    x3[i] = tempr;
	    x3[i+M] = tempi;
	}
	xx {xv[mem_mm]+=4*M*N+2*M;xv[add_mm]+=4*M*N;xv[mul_mm]+=4*M*N;}

	/* solve y = Rc slicing in each step */
	/* debugPrintVectorC( M, x3, "x3 = ["); */
	csolve_slice( G3, u3_, x3, M, M );
	/* debugPrintVectorC( M, u3_, "u3_ = ["); */
	/* debugPrintVectorI( M, permut, "permut = ["); */

	/* reorder estimated symbols */
	for( i=1; i<=M; i++ ) {
	    u3[permut[i]] = u3_[i];
	    u3[permut[i]+M] = u3_[i+M];
	}

	xx {xv[mem]+=4*M;}
	
	/* debugPrintVectorC( M, u3, "u3p = ["); */

	/* demodulate u3 */
	for( i=1; i<=M; i++ ) {
	    demodulate( u3[i], u3[i+M], infopointer );
	    infopointer ++;
	}

	xx {xv[mem]+=2*M;xv[add]+=5*M;} /* demodulate ops */

	infopointer = ( infopointer-M ) + ( M*iBitsPerSymbol );
    }

    /* count number of errors in this frame */

    BLERflag = 0;
    symindex = 0;
    symflag = 0;

    for( i=1; i<=iNumberInfoBitsFrame; i++ ) {
	/* print differences between estimated and transmitted bits */
	if (debug && ((cDebug & 256) == 256)) debug256(i);

	symindex ++;
	if( pEstimatedBits[i] != pSourceBits[i] ) {
	    biterrors ++;
	    if ( BLERflag == 0 ) {
		blockerrors ++;
		BLERflag = 1;
		if( ScreenFeedback ) {
		    printf( "Block error number: %d found!!\n", blockerrors );
		    printf( "  Frames simulated so far: %d\n", iFrameCounter );
		    printf( "  Current BLER = %1.3e\n", \
			    ( (float) blockerrors / (float) iFrameCounter) );
		}
	    }

	    if (debug && ((cDebug & 1024) == 1024)) debug1024(10,0);

	    if ( symflag == 0 ) {
		symerrors ++;
		symflag = 1;
	    }
	}
    }

    return 0;
}
