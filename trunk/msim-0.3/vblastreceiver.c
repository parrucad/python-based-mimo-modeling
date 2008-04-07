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

/* $Id: vblastreceiver.c,v 1.6 2004/08/16 21:54:28 miguel Exp $ */

/* MIMO simulation
   Miguel Bazdresch

   VBLAST Receiver function
   This function takes the output of the channel, namely pReceivedReal
   and pReceivedImag, and applies the VBLAST algorithm to it.
   Hseedsave is used to reproduce the H matrices when perfect
   knowledge of them is assumed.  The outputs of this function are
   pEstimatedBits and BER, BLER and SER.  The pointers pHReal and
   pHImag are reused from function channel(). */

int vblastreceiver(void) 
{
    /* function prototypes */
    double slice( double a );
    int    argmin( double **G, int a, int b );
    void   demodulate( double x, double y, int infopointer );
    void   debug64( int a1, int a2 );
    void   debug256( int a1 );

    /* local variables */
    int    BLERflag, symflag, symindex;
    int    infopointer, cols;
    int    column=0, i, index;
    double *ski, *skipm, **r, aki, akipm;

    /* infopointer helps reconstruct the frame */
    infopointer = 1;

    r=pReceived;

    /* initialize rc */
    for( i=1; i<=M; i++ )
	rc[i] = i;

    xx {xv[mem]+=M;}

    /* clear k */
    ivector_zero( k, M );
    
    xx {xv[mem_mz]+=M;}
    
    /* copy pH to pHWork */
    if( iReceiverType == 1 ) cols=M2;
    else cols = M;
    matrix_copy( pHWork, pH, N, cols );

    xx {xv[mem_mc]+=4*N*cols;}
    
    for( i=1; i<=M; i++ ) {
	
	/* printf("i = %d\n",i); */
	
	/* debugPrintMatrixC(N,M,pHWork,"hw = ["); */
    
	/* step 9b */
	pinv_vblast( pHWork, i );
	
	/* debugPrintMatrixC(M-i+1,N,storedpinvs[i],"mp = ["); */
	
	/* step 9c
	 * depending on decision order, create sequence k[i] */
	switch( exp_vb_ordering ) {
	    case 0:
		k[i] = argmin( storedpinvs[i], M-i+1, N );
		/* printf("ki = %d\n", k[i]); */
		break;
	    case 1:
		/* fix this: column is undefined before i start
		 * processing the rx vectors */
		k[i] = snrorder( storedpinvs[i], column );
		break;
	    case 2:
		k[i] = randorder( i );
		break;
	    case 3:
		k[i] = fixedorder( i );
		break;
	    default:
		return 1;
	}

	/* find original ordering of columns */
	if( i ==  1 )
	    origor[1] = k[i];
	else
	    origor[i] = rc[k[i]];
	/* debugPrintVectorI(M,rc,"rc = ["); */
	/* debugPrintVectorI(M,origor,"origor = ["); */

	xx{xv[mem]+=3;}

    }
    
    /* Apply VBLAST to info columns of pReceived */
    for ( column=(Lt+1); column<=Lt+L; column++ ) {
	
	/* debug - print received vector */
	/* printf("r = [ "); */
	/* for( i=1; i<=N; i++ ) */
	    /* printf("%1.5f + j*%1.5f, ",r[i][column],r[i+N][column]);
	     * */
	/* printf(" ];\n"); */
	    
	for ( i=1; i<=M; i++ ) {
	    /* repeat for every symbol in the column */
	    ski = storedpinvs[i][k[i]];
	    skipm = storedpinvs[i][k[i]+M-i+1];

	    /* debugPrintVectorR(N,storedpinvs[i][k[i]],"ski = ["); */
	    /* debugPrintVectorR(N,storedpinvs[i][k[i]+M-i+1],\ */
		    /* "skip = [");  */
	    
	    /* step 9e */
	    yReal = 0.0;
	    yImag = 0.0;

	    for( index=1; index<=N; index++ ) {
		
		double skii,skipmi,ri,ripn;
		
		skii = ski[index];
		skipmi = skipm[index];
		ri = r[index][column];
		ripn = r[index+N][column];
		yReal += skii*ri - skipmi*ripn;
		yImag += skipmi*ri + skii*ripn;
	    }
	    
	    xx {xv[mul_mm]+=4*N;xv[add_mm]+=4*N;xv[mem_mm]+=4*N;}
	    
	    /* step 9f */
	    aki = a[origor[i]] = slice( yReal );
	    akipm = a[origor[i]+M] = slice( yImag );
	    
	    /* printf("y1 = %1.5f, y3 = %1.5f\n",yReal,yImag); */
	    /* printf("a1 = %1.5f, a2 = %1.5f\n",aki,akipm); */
	    /* debugPrintVectorI(M,origor,"origor = ["); */
	    /* debugPrintVectorC(M,a,"a = ["); */
	    
	    xx {xv[mem]+=2;xv[add]+=3;}
	    
	    /* step 9g */
	    for( index=1; index<=N; index++ ) {
		
		double pHi, pHipn;
		
		pHi = pH[index][origor[i]];
		pHipn = pH[index+N][origor[i]];
		r[index][column] -= ( aki*pHi - akipm*pHipn );
		r[index+N][column] -= (akipm*pHi + aki*pHipn);
	    }
	    
	    xx {xv[mul_mm]+=4*N;xv[add_mm]+=4*N;xv[mem_mm]+=6*N;}
	    if (debug && ((cDebug & 64) == 64)) debug64(column, i);
	    
	    /* debugPrintVectorC(M,a,"a = ["); */
	    
	} /* end of loop - column of received matrix completely processed */

	/* After estimating a whole column of the received matrix, convert
	   the symbols back to bits. */
	for( index=1; index<=M; index++ ) {
	    demodulate( a[index], a[index+M], infopointer );
	    infopointer++;
	}
	
	xx {xv[mem]+=2*M;xv[add]+=5*M;}
	
	infopointer = (infopointer-M)+(M*iBitsPerSymbol);
    }
    
    /* count number of errors in this frame */
    BLERflag = 0;
    symindex = 0;
    symflag = 0;
    
    for( index=1; index<=iNumberCodedBitsFrame; index++ ) {
	symindex++;
	if( pEstimatedBits[index] != pSourceBits[index]) {
	    biterrors++;
	    if ( BLERflag == 0) {
		blockerrors++;
		BLERflag = 1;
		
		if( ScreenFeedback) {
		    printf( "Block error number: %d found!!\n", blockerrors); 
		    printf( "  Frames simulated so far: %d\n", iFrameCounter); 
		    printf( "   Current BLER = %1.5f\n", \
			    ((float)blockerrors/(float)iFrameCounter)); 
		}
	    }
	    
	    if (symflag == 0) {
		symerrors++;
		symflag = 1;
	    }
	    
	    /* print differences between estimated and transmitted bits */
	    if ((cDebug & 256) == 256) debug256(index);
	}
	
	if (symindex == iBitsPerSymbol) { symindex = 0; symflag = 0;}
    }

    return 0;
}
