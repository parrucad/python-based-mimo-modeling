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

/* $Id: vb-reduced.c,v 1.10 2004/08/16 17:29:58 miguel Exp $ */

int vb_reduced( void ) {

    /* function prototypes */
    void debug1024( int a1, int column );
    void tQR_trans( int m, int n, double **A, double **B, double **C );
    void demodulate( double a, double b, int c );

    /* local variables */
    int    column, i, j, l;
    int    infopointer, index;
    int    BLERflag, symflag, symindex;
    double tempr, tempi, *ski, *skipm;
    double aki, akipm;

    /* infopointer is used to rebuild the data bit stream, once symbols
     * have been estimated. */
    infopointer = 1;
    
    /* NEEDFIX: 4-QAM is not tested - need several updates */
    
    /*
     * The first section of the receiver is done only once per
     * block. In this section G is reduced and G3 is calculated.
     * These are used to decode each received vector.
     *
     */

    /* Find G3*Q=G. If (A,B)=thinqr(G') then Q=A', G3=B'
     * G=MxN; G3=MxM; Q=MxN */
    ctQR( N, M, pH, Q1, G3 );

    /* copy G3 to G3_org */
    matrix_copy( G3_org, G3, M, M );
    xx {xv[mem_mc]+=4*M*M;}

    /* debugPrintMatrixC(M,N,Q1,"Q = ["); */

    /* initialize rc */
    for( i=1; i<=M; i++ )
	rc[i] = i;

    xx {xv[mem]+=M;}

    /* clear k[i] */
    ivector_zero( k, M );
    
    xx {xv[mem_mz]+=M;}

    /* find M pseudoinverses */
    
    /* clear storedpinvs
     * don't count it -- it's only needed when cycling receivers */
    for( i=1; i<=M; i++ )
	for( j=1; j<=M; j++ )
	    for( l=1; l<=M; l++ ) {
		storedpinvs[i][j][l] = 0;
		storedpinvs[i][j+M][l] = 0;
	    }

    for( i=1; i<=M; i++ ) {

	/* printf("i = %d\n",i); */

	/* debugPrintMatrixC(M,M,G3,"hw = ["); */
	
	pinv_vblast( G3, i );

	/* debugPrintMatrixC(M-i+1,M,storedpinvs[i],"mp = ["); */
	
	k[i] = argmin( storedpinvs[i], M-i+1, M );

	/* printf("ki = %d\n", k[i]); */
	
	/* find original ordering of columns */
	if( i ==  1 )
	    origor[1] = k[i];
	else
	    origor[i] = rc[k[i]];

	/* debugPrintVectorI(M,origor,"origor = ["); */
	
	xx{xv[mem]+=3;}

    }

    /* 
     * Here ends the first section. We now have all the data we need.
     * The next section is repeated for each received vector.
     * 
     */

    for ( column=(Lt+1); column<=(Lt+L); column++ ) {
	
	if (debug && ((cDebug & 1024) == 1024)) debug1024(6,column);
	
	/* x3 = Q' * current column
	 * x3 is 1xM2 */

	/* debug - print received vector */
	/* printf("r = [ "); */
	/* for( i=1; i<=N; i++ ) */
	    /* printf("%1.5f + j*%1.5f, ",pReceived[i][column],\ */
		    /* pReceived[i+N][column]); */
	/* printf(" ];\n"); */
	
	/* debugPrintMatrixC(M,N,Q1,"Q = ["); */

	for( i=1; i<=M; i++ ) {
	    tempr = 0.0;
	    tempi = 0.0;
	    for( j=1; j<=N; j++ ) {
		tempr += Q1[j][i] * pReceived[j][column] + \
		         Q1[j+N][i] * pReceived[j+N][column];
		tempi += Q1[j][i] * pReceived[j+N][column] - \
		         Q1[j+N][i] * pReceived[j][column];
	    }
	    x3[i] = tempr;
	    x3[i+M] = tempi;
	}

	/* debugPrintVectorC(M,x3,"x3 = ["); */
	
	xx {xv[mem_mm]+=4*M*N+2*M;xv[add_mm]+=4*M*N;xv[mul_mm]+=4*M*N;}

	/* apply vblast() to x3 and G3 to estimate u3 */
	for( i=1; i<=M; i++ ) {

	    /* printf("i = %d\n\n",i); */

	    /* repeat for every symbol in the column */
	    ski = storedpinvs[i][k[i]];
	    skipm = storedpinvs[i][k[i]+M-i+1];
	    
	    /* debugPrintVectorR(M,storedpinvs[i][k[i]],"ski = ["); */
	    /* debugPrintVectorR(M,storedpinvs[i][k[i]+M-i+1],"skip =
	     * ["); */

	    /* debugPrintVectorR(M2,x3,"x3 = ["); */
	    
	    /* step 9e */
	    yReal = 0.0;
	    yImag = 0.0;
	    for( index=1; index<=M; index++ ) {
		
		double skii,skipmi,ri,ripn;
		
		skii = ski[index];
		skipmi = skipm[index];
		ri = x3[index];
		ripn = x3[index+M];
		yReal += skii*ri - skipmi*ripn;
		yImag += skipmi*ri + skii*ripn;
	    }
	    
	    xx {xv[mul_mm]+=4*M;xv[add_mm]+=4*M;xv[mem_mm]+=4*M;}

	    /* step 9f */
	    aki = u3[origor[i]] = slice( yReal );
	    akipm = u3[origor[i]+M] = slice( yImag );

	    /* printf("y1 = %1.5f, y3 = %1.5f\n",yReal,yImag); */
	    /* printf("a1 = %1.5f, a2 = %1.5f\n",aki,akipm); */

	    /* debugPrintVectorC(M,u3,"u3 = ["); */
	    
            xx {xv[mem]+=2;xv[add]+=3;}
	    
	    /* printf("g3c = [ "); */
	    /* for( index=1; index<=M; index++ ) */
		/* printf("%1.5f + j*%1.5f, ",\ */
			/* G3_org[index][origor[i]],\ */
			/* G3_org[index+M][origor[i]]); */
	    /* printf(" ];\n"); */

	    /* step 9g */
	    for( index=1; index<=M; index++ ) {
		
		double G3i, G3ipn;

		G3i = G3_org[index][origor[i]];
		G3ipn = G3_org[index+M][origor[i]];
		x3[index] -= ( aki*G3i - akipm*G3ipn );
		x3[index+M] -= ( akipm*G3i + aki*G3ipn );
	    }

	    /* debugPrintVectorR(M2,x3,"x3 = ["); */
	    
	    xx {xv[mul_mm]+=4*M;xv[add_mm]+=4*M;xv[mem_mm]+=6*M;}
	    
	}

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
