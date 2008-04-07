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

/* $Id: mlreceiver.c,v 1.3 2004/08/05 22:14:36 miguel Exp $ */

/* MIMO Simulation
   Miguel Bazdresch


   This file contains a ML MIMO receiver function based on Agrell.
   This function takes the output of the channel, namely pReceived,
   and applies the ML algorithm to it.
   The outputs of this function are pEstimatedBits and BER, BLER and SER.
   The main algorithm is as follows:
   
   i.     repeat for number of columns in block
     i.ii   repeat for number of info bits
       i.ii.i execute ML on the column
   ii.    count errors
   */

int mlreceiver(void)
{
    /* function prototypes */
    void debug1024( int a1, int column );
    void inverse_tqr( double **A, double **Ai, int r, int c );
    void lll( void );
    void rmult_mat_mat_tran( double **A, double **B, double **C, \
	    int a, int b, int c );
    void inverse_trans( double **A, double **Ai, int col );
    void tQR_trans( int m, int n, double **A, double **B, double **C );
    void inverse_tqr_trans( double **A, double **B, int m, int n );
    void decode( void );
    void demodulate( double a, double b, int c );

    /* local variables */
    int    row, col, column, i, j;
    int    infopointer;
    int    BLERflag, symflag, symindex;
    double temp, const_length;

    /* infopointer is used to rebuild the data bit stream, once symbols
     * have been estimated. */
    infopointer = 1;
    
    /* calculate constellation length -
     * it is needed to calculate the translate vector */ 
    if( cConstellationType == 0 ) const_length = e3;       /* 16-QAM */
    else if( cConstellationType == 1 ) const_length = e1;  /* 4-QAM  */

    /* NEEDFIX: 4-QAM is not tested - need demod update */
    
    /*
     * The first section of the receiver is done only once per
     * block. In this section G is reduced and H3 is calculated.
     * These are used to decode each received vector.
     *
     */

    /* Construct G from H: G = 2*e1*H'
     * H is N2xM2, G is M2xN2 */
    temp = 2*e1;
    for( row=1; row<=M2; row++ )
	for( col=1; col<=N2; col++ )
	    G[row][col] = pH[col][row] * temp;
    
    xx {xm[mem]+=2*M2*N2;xm[mul]+=M2*N2;}
    if(debug && ((cDebug & 1024) == 1024)) debug1024(1,0);

    /* create translate vector:
     * this vector is used to move the received point to the lattice
     * used for ML estimation.
     * translate is 1xN2, H is M2xN2 */
    
    for( i=1; i<=N2; i++ ) {
	
	double temp;
	
	temp = 0.0;
	for( j=1; j<=M2; j++ )
	    temp += pH[i][j];
	
	translate[i] = e3 * temp;
    }
    
    xx {xm[mem_mm]+=N2+N2*M2;xm[mul_mm]+=N2;xm[add_mm]+=M2*N2;}
    
    /* If configured, find G = LLL(G).
     * G is M2xN2; the LLL is stored in the same matrix G */
    
    if( DoLLL ) {

	/* Need inverse of G, see comments below
	 * Note that all operations regarding calculation of W are
	 * not counted towards the complexity totals, because W could
	 * in principle be found by lll() without extra cost */
	
	/* transpose G to tempM since inverse_tqr destroys it,
	 * and we want to find the inverse of G', not G */
	for( i=1; i<=M2; i++ )
	    for( j=1; j<=N2; j++)
		tempM[j][i] = G[i][j];
	
	/* R1 is the inverse of G' */
	inverse_tqr( tempM, R1, N2, M2 );
	
	if (debug && ((cDebug & 1024) == 1024)) debug1024(2,0);
	
	lll();
	
	if (debug && ((cDebug & 1024) == 1024)) debug1024(3,0);

	/* find W = LLL(G) * inv(G)
	 * G=M2xN2, inv(G)=N2xM2, LLL(G)=M2xN2, W=M2xM2
	 * I will re-use matrix bstar to store W */
	
        rmult_mat_mat_tran( bstar, G, R1, M2, N2, N2 );
	
    }
    
    /* Find G3*Q=G. If (A,B)=thinqr(G') then Q=A', G3=B'
     * G=M2xN2; G3=M2xM2; Q=M2xN2 */

    for( row=1; row<=M2; row++ )
	for( col=1; col<=row; col++ )
	    G3_t[row][col] = 0;

    xx {xm[mem_mz]+=(int)M2*M2/2;}
    
    tQR_trans( N2, M2, G, Q_t, G3_t );

    if (debug && ((cDebug & 1024) == 1024)) debug1024(4,0);
    
    /* Now calculate H3=inv(G3)
     * G3 = M2xM2;  H3 = M2xM2 */
    
#ifdef NR_LICENSED
    inverse_trans( G3_t, H3, M2 );
#else
    inverse_tqr_trans( G3_t, H3, M2, M2 );
#endif
    
    /* 
     * Here ends the first section. We now have G and H3. The next
     * section is repeated for each received vector.
     * 
     */
    
    for ( column=(Lt+1); column<=(Lt+L); column++ ) {
	
	if (debug && ((cDebug & 1024) == 1024)) debug1024(6,column);
	
	/* x3 = current column times Q_t
	 * x3 is 1xM2 */
	
	for( j=1; j<=N2; j++ )
	    pReceived[j][column] += translate[j];

	xx {xm[mem]+=2*N2;xm[add]+=N2;}
	
	for( i=1; i<=M2; i++ ) {
	    temp = 0.0;
	    for( j=1; j<=N2; j++ ) {
		temp += pReceived[j][column] * Q_t[j][i];
	    }
	    x3[i] = temp;
	}

	xx {xm[mem_mm]+=2*M2*N2+M2;xm[add_mm]+=M2*N2;xm[mul_mm]+=M2*N2;}
	if (debug && ((cDebug & 1024) == 1024)) debug1024(7,0);

	/* apply decode() to x3 and H3 */
	decode();

        /* If LLL reduction was performed on G, then we need to adjust u3
	 * by multiplying it by W (which is stored in bstar: u3 <- u3*W
	 * u3 is 1xM2, W is M2xM2 */
	if(DoLLL) {
	   for( i=1; i<=M2; i++ ) {
	       tempV[i] = u3[1]*bstar[1][i];
	       for( j=2; j<=M2; j++)
		   tempV[i] += u3[j] * bstar[j][i];
	   }
	   for( i=1; i<=M2; i++ )
	       u3[i] = rint( tempV[i] );
	   
	   xx {xl[mem]+=3*M2*M2;xl[add]+=M2*(M2-1)+2*M2;}
	   
	}

	/* slice u3 */
	for( i=1; i<=M2; i++ ) {
	    if( u3[i] > 3 ) u3[i] = 3;
	    else if( u3[i] < 0 ) u3[i] = 0;
	}
	
        xx {xm[mem]+=2*M2;xm[add]+=M2;}
	
	/* once u3 is sliced, demodulate it */
	for( i=1; i<=M; i++ ) {
	    demodulate( u3[i], u3[i+M], infopointer );
	    infopointer ++;
	}
	
	xx {xm[mem]+=2*M;xm[add]+=5*M;} /* demodulate ops */
	
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
