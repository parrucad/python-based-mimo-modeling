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

/* $Id: support.c,v 1.11 2004/08/16 21:54:05 miguel Exp $ */

/* MIMO simulation */

/* Miguel Bazdresch */

/* This function compresses H to remove column ki
 * and updates rc */
void strip_column( double **H, int rows, int ki, int i ) {

    int j, colstomove;

    colstomove = M - i - ki + 1;

    for( j=1; j<=colstomove; j++ ) {
	matrix_col_copy( H, H, rows, ki+j-1, ki+j );
	rc[ki+j-1] = rc[ki+j];
    }

    xx {xp[mem_mc]+=colstomove*(rows*4+2);}
}

/* This function compresses H to remove column ki */
void remove_column( double **H, int rows, int ki, int i ) {

    int j, colstomove;

    colstomove = M - i - ki + 1;

    for( j=1; j<=colstomove; j++ ) {
	matrix_col_copy( H, H, rows, ki+j-1, ki+j );
    }

    xx {xp[mem_mc]+=colstomove*(rows*4);}
}

/* pinv_vblast
 * returns pinv of H in storedpinvs[i]
 * H is MxN, storedpinvs[i] is NxM
 * method of calculation depends on iReceiverType
 * strips columns known to be zero */

void pinv_vblast( double **H, int i ) {

    /* function prototypes */
    void mult_mat_mat_tran(double **dest, double **M1, double **M2, int m, int r, int n);
    void mult_mat_tran_mat(double **dest, double **M1, double **M2, int m, int r, int n);
    void inverse(double **A, double **Ai, int col);
    void tQR(int M, int N, double **QR, double **Q, double **R);
    void solve(double **R, double *xls, double *b, int N);
    /* debug prototypes */
    void debug4096(int i, int Mr);
    /* local variables */
    int col, row, rc2;
#ifdef NR_LICENSED
    double wmax, wmin;
#endif

    /* strip input matrix of zeroed columns */
    rc2 = M - i + 1;

    if( i != 1 ) {
	if( iReceiverType == 5 )
	    strip_column( H, M, k[i-1], i-1 );
	else
	    strip_column( H, N, k[i-1], i-1 );
    }

    /* debugPrintMatrixC(N,rc2,H,"hr = ["); */

    /* do pinv */
    switch( iReceiverType ) {
	case 1:  /* svd */
#ifdef NR_LICENSED
	    nr_svdcmp( pZ, N2, red_cols, pW, pV );
	    /* identify and zero w[j] that are too small, and invert the others */
	    wmax = 0.0;
	    for( col=1; col<=red_cols; col++ )
		if( pW[col] > wmax )
		    wmax = pW[col];
	    /* change 2.2e-16 according to architecture -- see
	     * NR book for details */
	    wmin = N2*wmax*2.2e-16;
	    xx {xs[mem]+=red_cols;xs[add]+=red_cols;xs[mul]++;}
	    for( col=1;col<=red_cols;col++ ) {
		if( pW[col] < wmin )
		    pW[col] = 0.0;
		else {
		    pW[col] = 1.0 / pW[col];
		    xx {xs[div]++;}
		}
	    }
	    xx {xs[mem]+=red_cols*2;xs[add]+=red_cols;}
	    /* Multiply V by 1/W. Result is stored in V */
	    for( col=1;col<=red_cols;col++ ) {
		double temp=pW[col];
		for( row=1;row<=red_cols;row++ )
		    pV[row][col] *= temp;
	    }
	    xx {int rc=red_cols;xs[mem_mm]+=(rc+2*rc*rc);xs[mul_mm]+=rc*rc;}
	    /* Multiply T=VU'
	     * V is red_colsxred_cols; U' is red_cols x N2
	     * T is red_cols x N2 */
	    rmult_mat_mat_tran( pT, pV, pZ, red_cols, red_cols, N2); 
	    xx {int rc=red_cols;xs[mem_mm]+=2*rc*N2*(rc+1);
		xs[mul_mm]+=rc*rc*N2;xs[add_mm]+=rc*rc*N2;}
		/* print calculated singular values to file */
		if(debug && ((cDebug & 128) == 128)) debug128(0, pW);
#else
		printf("svd receiver not implemented, we're working on it...\n");
		exit(1);
#endif
		break;
	case 2:  /* formula */
#ifdef NR_LICENSED
		/* it calculates the pinv of H using the formula
		 *  A^+=(A^T*A)^(-1)*A^T */
		mult_mat_tran_mat( pTT, pZ, pZ, N, rc2, rc2 );
		
		for( i=1; i<=rc2; i++ )
		    for( j=1; j<=rc2; j++) {
			pTT[i][j+rc2] =- pTT[i+rc2][j];
			pTT[i+rc2][j+rc2] = pTT[i][j];
		    }
		xx {xf[mem_mc]+=4*rc2*rc2;xf[add]+=rc2*rc2;}

		inverse(pTT, pTTinv, red_cols);
		mult_mat_mat_tran(pT,pTTinv,pZ,rc2,rc2,N);
		/* mat. mult. counters */
		xx {xf[mem_mm]+=8*(N*rc2*(rc2+1));}
		xx {xf[mul_mm]+=8*rc2*rc2*N;xf[add_mm]+=8*rc2*rc2*N;}
#else
		printf("formula receiver not implemented, we're working on it...\n");
		exit(1);
#endif
		break;
	case 3: /* thin QR */
		matrix_copy( pZ, H, N, rc2 );

		ctQR( N, rc2, pZ, Q, R );
		/* debugPrintMatrixC(N,rc2,Q,"Q = ["); */
		/* debugPrintMatrixC(rc2,rc2,R,"R = ["); */
		/* now find pinv(pZ) */
		for( row=1; row<=N; row++ ) {
		    for( col=1; col<=rc2; col++ ) {
			b[col]     =  Q[row][col];   /* copy row i of Q to b */
			b[col+rc2] = -Q[row+N][col]; /* transpose of Q */
		    }
		    
		    csolve1( R, xls, b, rc2, rc2 );
		    
		    for( col=1; col<=rc2*2; col++ )
			/* copy xls to corresponding column of psinv */
			storedpinvs[i][col][row] = xls[col]; 
		}

		xx {xh[mem_mc]+=12*N*rc2;xh[add]+=N*rc2;}

		/* debugPrintMatrixC(rc2,N,storedpinvs[i],"mp = ["); */
		break;
		
	case 4: /* QR with updates */
		/* printf("i = %d ; rc2 = %d\n", i, rc2); */
		if( i == 1 ) {  /* do full QR */
		    matrix_copy( pZ, H, N, rc2 );
		    QR( Q, R, pZ, N, M );
		    /* debugPrintMatrixC( N, N, Q, "Q = ["); */
		    /* debugPrintMatrixC( N, rc2, R, "R = ["); */
		}
		else {

		    /* remove column k[i-1] from R */
		    remove_column( R, N, k[i-1], i-1 );
		    /* debugPrintMatrixC( N, N, Q, "Qn = ["); */
		    /* debugPrintMatrixC( N, rc2, R, "Rn = ["); */
		    updateqr( Q, R, N, rc2, k[i-1]);
		    /* debugPrintMatrixC( N, N, Q, "Qu = ["); */
		    /* debugPrintMatrixC( N, rc2, R, "Ru = ["); */
		}

		for( row=1; row <= N; row++ ) {
		    for( col=1; col<=N; col++ ) {
			b[col]   =  Q[row][col];   /* copy row i of Q to b */
			b[col+N] = -Q[row+N][col]; /* transpose of Q */
		    }

		    xx {xq[mem_mc]+=4*N;xq[add]+=N;}

		    csolve1( R, xls, b, N, rc2 );

		    for( col=1; col<=rc2*2; col++ )
			/* copy xls to corresponding column of psinv */
			storedpinvs[i][col][row] = xls[col];
		}
		/* debugPrintMatrixC(rc2,N,storedpinvs[i],"mp = ["); */
		break;

	case 5: /* LS-BLAST */
		/* pZ is N2*red_cols, Q is N2*red_cols, R is red_cols*red_cols */
		matrix_copy( pZ, H, M, rc2 );

		ctQR( M, rc2, pZ, Q, R );
		/* debugPrintMatrixC(N,rc2,Q,"Q = ["); */
		/* debugPrintMatrixC(rc2,rc2,R,"R = ["); */
		/* now find pinv(pZ) */
		for( row=1; row<=M; row++ ) {
		    for( col=1; col<=rc2; col++ ) {
			b[col] = Q[row][col]; /* copy row i of Q to b */
			b[col+rc2] =- Q[row+M][col]; /* transpose of Q */
		    }

		    csolve1( R, xls, b, rc2, rc2 );

		    for( col=1; col<=rc2*2; col++ )
			/* copy xls to corresponding column of psinv */
			storedpinvs[i][col][row] = xls[col]; 
		}

		xx {xh[mem_mc]+=12*M*rc2;xh[add]+=M*rc2;}

		break;

	default:
		printf( "Error: receiver %d not implemented, exiting...\n", iReceiverType);
		exit(1);
    }
}

int fixedorder(int i) {
    return i;
}

int randorder(int i) {

    double random,bin;
    static int order[30];
    int temp,flag,k,j;

    if(i==1) {
	bin=(double)(RAND_MAX)/(double)M;
	for(j=1;j<=M;j++) order[j]=0;
	for(j=1;j<=M;j++) {
	    /* repopulate order[] */
	    do {
		flag=0;
		random=(double)(rand());
		if (random==0) random=1;

		temp=(int)(ceil(random/bin));
		for(k=1;k<j;k++)
		    if(temp==order[k]) flag=1;
	    } while(flag==1);
	    order[j]=temp;
	}
    }
    return order[i];
}

int snrorder(double **G, int rcol) {

    /* local variables */
    int alreadyselected, i, row, col, ans=0;
    double largest, smallest, noiseR, noiseI, noiseMag; 

    /* find the ro for each symbol */

    smallest = 99999;
    largest = 0;
    for(row=1;row<=M;row++) {
	/* first check that this row has not been already selected. */
	alreadyselected = 0;
	for(i=1;i<=M;i++) {
	    if (k[i] == row) {
		alreadyselected = 1;
		break;
	    }
	}
	/* if this row has not been selected, then process it. */
	if (alreadyselected == 0) {
	    noiseR=0;
	    noiseI=0;
	    for(col=1;col<=N;col++) {
		noiseR+=G[row][col]*noisematrix[col][rcol]-G[row+M][col]*noisematrix[col+N][rcol];
		noiseI+=G[row][col]*noisematrix[col+N][rcol]+G[row+M][col]*noisematrix[col][rcol];
	    }
	    noiseMag=noiseR*noiseR+noiseI*noiseI;
	    if (noiseMag < smallest) {
		smallest = noiseMag;
		ans = row;
	    }
	}
    }
    return ans;
}

/* argmin */
/* This function examines the rows of G, and returns the row number with
   the minimum squared magnitude. It only examines the rows that have not
   been examined before - that is, rows for which k[row] != 0.
   G is size M2xN2. */
int argmin( double **G, int Rows, int Cols ) {

    /* local variables */
    int    row, col, ans = 0;
    double smallest, sum, *gr, *grpm;

    smallest = 1e10;
    
    for( row=1; row<=Rows; row++ ) {
	sum = 0.0;
	gr=G[row];
	grpm=G[row+Rows];
	for( col=1; col<=Cols; col++ ) {
	    double a,b;
	    a = gr[col];
	    b = grpm[col];
	    sum += a*a+b*b;
	}

	xx {xv[mem]+=2*Cols;xv[add]+=2*Cols;xv[mul]+=2*Cols;}

	if( sum < smallest ) {
	    smallest = sum;
	    ans = row;
	}

	xx {xv[add]++;}
    }
    return ans;
}

int argminR( double **G, int Rows, int Cols ) {

    /* local variables */
    unsigned int i, row, col, ans = 0, alreadyselected;
    double smallest, sum;

    smallest = 1e10;
    
    for( row=1; row<=Rows; row++ ) {
	/* first check that this row has not been already selected. */
	alreadyselected = 0;
	for( i=1; i<=Rows; i++ ) {
	    xx {xv[mem]++;}
	    
	    if( k[i] == row ) {
		alreadyselected = 1;
		break;
	    }
	}
	
	/* if this row has not been selected, then process it. */
	if( alreadyselected == 0 ) {
	    sum = 0.0;
	    for( col=1; col<=Cols; col++ ) {
		sum += G[row][col]*G[row][col];
	    }
	
	    xx {xv[mem]+=2*Cols;xv[add]+=Cols;xv[mul]+=Cols;}
	
	    if( sum < smallest ) {
		smallest = sum;
		ans = row;
	    }
	
	    xx {xv[add]++;}
	}
    }
    return ans;
}

/* slice */
/* this functions quantizes a number to its closest coordinate in the
 * selected constellation. */

double slice(double a) 
{
    /* which constellation to use depends on global variable cConstellationType */
    switch( cConstellationType ) {
	case 0:  /* 16-QAM */
	    if( a >= t3 ) {
		return e3;
	    }
	    if( a >= t2 ) {
		return e1;
	    }
	    if( a >= t1 ) {
		return -e1;
	    }
	    return -e3;
	    break;
	case 1: /* 4-QAM */
	    if( a>= 0 )
		return e1;
	    else
		return -e1;
	    break;
	default:
	    exit( 1 );
    }
    return 0;
}

/* demodulate */
/* This function returns as int the bits that correspond to symbol (x, y) */
void demodulate(double x, double y, int infopointer) {

    int a = 0, b = 0, c = 0, d = 0;
    double coord1, coord2, coord3, coord4;

    /* 
     * The expected coordinates of the data points depend on 
     * whether the receiver is BLAST or ML
     *
     */
    if( iReceiverType == 0 ) {
	coord1 = 0;
	coord2 = 1;
	coord3 = 2;
	coord4 = 3;
    } else {
	coord1 = -e3;
	coord2 = -e1;
	coord3 = e1;
	coord4 = e3;
    }
    
    /* The symbol to return depends on the constellation used. */
    switch (cConstellationType) {
	case 0: /* 16-QAM */
	    if (x == coord1) {
		if (y == coord1) {
		    a = 1;
		    b = 1;
		    c = 1;
		    d = 1;
		}
		if (y == coord2) {
		    a = 0;
		    b = 1;
		    c = 1;
		    d = 1;
		}
		if (y == coord3) {
		    a = 0;
		    b = 1;
		    c = 0;
		    d = 1;
		}
		if (y == coord4) {
		    a = 1;
		    b = 1;
		    c = 0;
		    d = 1;
		}
	    }
	    if (x == coord2) {
		if (y == coord1) {
		    a = 1;
		    b = 0;
		    c = 1;
		    d = 1;
		}
		if (y == coord2) {
		    a = 0;
		    b = 0;
		    c = 1;
		    d = 1;
		}
		if (y == coord3) {
		    a = 0;
		    b = 0;
		    c = 0;
		    d = 1;
		}
		if (y == coord4) {
		    a = 1;
		    b = 0;
		    c = 0;
		    d = 1;
		}
	    }
	    if (x == coord3) {
		if (y == coord1) {
		    a = 1;
		    b = 0;
		    c = 1;
		    d = 0;
		}
		if (y == coord2) {
		    a = 0;
		    b = 0;
		    c = 1;
		    d = 0;
		}
		if (y == coord3) {
		    a = 0;
		    b = 0;
		    c = 0;
		    d = 0;
		}
		if (y == coord4) {
		    a = 1;
		    b = 0;
		    c = 0;
		    d = 0;
		}
	    }
	    if (x == coord4) {
		if (y == coord1) {
		    a = 1;
		    b = 1;
		    c = 1;
		    d = 0;
		}
		if (y == coord2) {
		    a = 0;
		    b = 1;
		    c = 1;
		    d = 0;
		}
		if (y == coord3) {
		    a = 0;
		    b = 1;
		    c = 0;
		    d = 0;
		}
		if (y == coord4) {
		    a = 1;
		    b = 1;
		    c = 0;
		    d = 0;
		}
	    }
	    pEstimatedBits[infopointer] = a;
	    pEstimatedBits[infopointer+M] = b;
	    pEstimatedBits[infopointer+2*M] = c;
	    pEstimatedBits[infopointer+3*M] = d;
	    break;
	case 1: /* 4-QAM */
	    if( x == e1 ) {
		if( y == e1 ) {
		    a = 1;
		    b = 1;
		} else {
		    a = 1;
		    b = 0;
		}
	    } else {
	       	if( y == e1 ) {
		    a = 0;
		    b = 1;
		} else {
		    a = 0;
		    b = 0;
		}
	    }
	    pEstimatedBits[infopointer] = a;
	    pEstimatedBits[infopointer+M] = b;
	    break;
	default:
	    exit( 1 );
    }
}

void winddown(void) 
{
    /* This function calculates final error rate and writes
       results in file errfile */

    unsigned long long int iTotalBits;
    double tb;

    iTotalBits = iFrameCounter*iNumberInfoBitsFrame;
    tb=(double)iTotalBits;
    /* Calculate power if cDebug = 8 */
    if (debug && ((cDebug & 8) == 8)) {
	EstNoisePwr = 2*EstNoisePwr/(iFrameCounter*(Lt+L)*N2);
	EstRxPwr = 2*EstRxPwr/(iFrameCounter*(Lt+L)*N2);
	EstSentPwr = 2*EstSentPwr/(iFrameCounter*(Lt+L)*M2);
	fprintf(powerfile, "Estimated Noise Power (per antenna) = %1.5f\n", EstNoisePwr);
	fprintf(powerfile, "Estimated Sent Power (all antennas)= %1.5f\n", EstSentPwr);
	fprintf(powerfile, "Estimated Rx Power (per antenna) = %1.5f\n", EstRxPwr);
	if (fNoisePower != 0.0) {
	    fprintf(powerfile, "SNR (Rx/Noise) = %1.5fdB\n", 10*log10(EstRxPwr/EstNoisePwr));
	}
	else {
	    fprintf(powerfile, "No noise in channel\n");
	}
    }
    /* find end time */
    end_time = time(&timex);
    run_time = end_time - initial_time + suspend_time;
    fprintf(errfile, "#SNR:%1.1f,errors:%u,blocks:%u\n", fNoisePower,biterrors,blockerrors);
    if (biterrors > 0) {
	BER = biterrors/(double)(iNumberInfoBitsFrame*iFrameCounter);
	BLER = blockerrors/(double)(iFrameCounter);
	SER = symerrors/(double)(iNumberInfoBitsFrame*M*iFrameCounter/iBitsPerSymbol);
    }
    if (biterrors == 0) {
	BER = 0.0;
	BLER = 0.0;
	SER = 0.0;
    }
    printf("Results: Run time = %lu seconds.\n", run_time);
    printf("Results: L = %d\n",L);
    printf("Results: Simulated %d frames.\n", iFrameCounter);
    printf("Results: Simulated %d bits.\n", iFrameCounter*iNumberInfoBitsFrame);
    printf("Results: Bit errors: %d\n", biterrors);
    printf("Results: Block errors: %d\n", blockerrors);
    printf("Results: BER = %1.10f\n", BER);
    printf("Results: BLER = %1.10f\n", BLER);
    printf("Results: SER = %1.10f\n", SER);
    if (firstwd) fprintf(errfile, "#dB\tBER\tBLER\n%1.1f\t%1.3g\t%1.3g",fNoisePower,BER,BLER);
    fprintf(outfile, "Simulated %d frames.\n", iFrameCounter);
    fprintf(outfile, "Simulated %Lu bits.\n", iTotalBits);
    fprintf(outfile, "L = %d\n", L);
    fprintf(outfile, "  SNR = %1.1f\n", fNoisePower);
    fprintf(outfile, "  BER = %1.10f\n", BER);
    fprintf(outfile, "  BLER = %1.10f\n", BLER);
    fprintf(outfile, "  SER = %1.10f\n\n", SER);
    /* complexity totals and reports */
    xx {
	double tadd,taddmm,Tadd,tmul,tmulmm,Tmul,tmem,tmemmm,tmemmc;
	double Tdiv, Tsrt,tmemmz,Tmem,TOTAL;
	tadd=xp[add]+xv[add]+xm[add]+xd[add]+xl[add]+xh[add]+xq[add]+xs[add]+xf[add];
	taddmm=xp[add_mm]+xv[add_mm]+xm[add_mm]+xd[add_mm]+xl[add_mm]+xh[add_mm]+xq[add_mm]+xs[add_mm]+xf[add_mm];
	Tadd=tadd+taddmm;
	tmul=xp[mul]+xv[mul]+xm[mul]+xd[mul]+xl[mul]+xh[mul]+xq[mul]+xs[mul]+xf[mul];
	tmulmm=xp[mul_mm]+xv[mul_mm]+xm[mul_mm]+xd[mul_mm]+xl[mul_mm]+xh[mul_mm]+xq[mul_mm]+xs[mul_mm]+xf[mul_mm];
	Tmul=tmul+tmulmm;
	tmem=xp[mem]+xv[mem]+xm[mem]+xd[mem]+xl[mem]+xh[mem]+xq[mem]+xs[mem]+xf[mem];
	tmemmm=xp[mem_mm]+xv[mem_mm]+xm[mem_mm]+xd[mem_mm]+xl[mem_mm]+xh[mem_mm]+xq[mem_mm]+xs[mem_mm]+xf[mem_mm];
	tmemmc=xp[mem_mc]+xv[mem_mc]+xm[mem_mc]+xd[mem_mc]+xl[mem_mc]+xh[mem_mc]+xq[mem_mc]+xs[mem_mc]+xf[mem_mc];
	tmemmz=xp[mem_mz]+xv[mem_mz]+xm[mem_mz]+xd[mem_mz]+xl[mem_mz]+xh[mem_mz]+xq[mem_mz]+xs[mem_mz]+xf[mem_mz];
	Tmem=tmem+tmemmm+tmemmc+tmemmz;
	Tdiv=xp[div]+xv[div]+xm[div]+xd[div]+xl[div]+xh[div]+xq[div]+xs[div]+xf[div];
	Tsrt=xp[srt]+xv[srt]+xm[srt]+xd[srt]+xl[srt]+xh[srt]+xq[srt]+xs[srt]+xf[srt];
	TOTAL=Tadd+Tmul+Tmem+Tdiv+Tsrt;
	/* printf("xharit = %-16.0f, xhmem =
	 * %-16.0f\n",(double)(xh[add]+xh[add_mm]+ \ */
		/* xh[mul]+xh[mul_mm]+xh[div]+xh[srt]),
		 * (double)(xh[mem]+xh[mem_mm]+ \ */
		/* xh[mem_mc]+xh[mem_mz])); */
	/* printf("xvarit = %-16.0f, xvmem =
	 * %-16.0f\n",(double)(xv[add]+xv[add_mm]+ \ */
		/* xv[mul]+xv[mul_mm]+xv[div]+xv[srt]),
		 * (double)(xv[mem]+xv[mem_mm]+ \ */
		/* xv[mem_mc]+xv[mem_mz])); */
	/* this line only once */
	if( firstwd ) {
	    if( IndexComFile )
		fprintf(comfile,"\n\n#M=%d; N=%d\n", M, N);
	    else
		fprintf(comfile,"#M=%d; N=%d\n", M, N);
	    fprintf(comfile,"#Receiver = %d\n",iReceiverType);
	    fprintf(comfile,"#L = %d\n",L);
	    if(iReceiverType==0) fprintf(comfile,"#DoLLL = %d\n",DoLLL);
	    fprintf(comfile,"#1\t\t2\t\t3\t\t4\t\t5\t\t6\t\t7\t\t8\t\t9\t\t10\t\t");
	    fprintf(comfile,"11\t\t12\t\t13\t\t14\t\t15\t\t16\t\t");
	    fprintf(comfile,"17\t\t18\t\t19\t\t20\t\t21\t\t22\t\t23\t\t24\t\t25\t\t2");
	    fprintf(comfile,"6\t\t27\t\t28\t\t29\t\t30\t\t");
	    fprintf(comfile,"31\t\t32\t\t33\t\t34\t\t35\t\t36\t\t37\t\t38\n");
	    fprintf(comfile,"#dB\t\tBER\t\tBLER\t\t");
	    fprintf(comfile,"add\t\tadd%%\t\tadd_mm\t\tadd_mm%%\t\taddtot");
	    fprintf(comfile,"\t\tmul\t\tmul%%\t\tmul_mm\t\tmul_mm%%\t\tmultot\t\tmem");
	    fprintf(comfile,"\t\tmem%%\t\tmem_mm\t\tmem_mm%%\t\tmem_mc\t\tmem_mc%%\t\t");
	    fprintf(comfile,"mem_mz\t\tmem_mz%%\t\tmemtot\t\tdiv\t\tsrt\t\t");
	    fprintf(comfile,"%%add\t\t%%mul\t\t%%mem\t\t%%div\t\t%%srt\t\tadd/b\t\t");
	    fprintf(comfile,"mul/bit\t\tmem/b\t\tdiv/bit\t\tsrt/bit\t\ttot/bit");
	    fprintf(comfile,"\t\tL\t\tM\t\tN\n");
	}
	fprintf(comfile,"%-16.1f%-16.4e%-16.4e",fNoisePower,BER,BLER);
	fprintf(comfile,"%-16.0f%-16.3f",tadd,tadd/Tadd);
	fprintf(comfile,"%-16.0f%-16.3f%-16.0f",taddmm,taddmm/Tadd,Tadd);
	fprintf(comfile,"%-16.0f%-16.3f%-16.0f",tmul,tmul/Tmul,tmulmm);
	fprintf(comfile,"%-16.3f%-16.0f%-16.0f",tmulmm/Tmul,Tmul,tmem);
	fprintf(comfile,"%-16.3f%-16.0f%-16.3f",tmem/Tmem,tmemmm,tmemmm/Tmem);
	fprintf(comfile,"%-16.0f%-16.3f%-16.0f",tmemmc,tmemmc/Tmem,tmemmz);
	fprintf(comfile,"%-16.3f%-16.0f%-16.0f",tmemmz/Tmem,Tmem,Tdiv);
	fprintf(comfile,"%-16.0f%-16.3f%-16.3f",Tsrt,Tadd/TOTAL,Tmul/TOTAL);
	fprintf(comfile,"%-16.3f%-16.3f%-16.3f",Tmem/TOTAL,Tdiv/TOTAL,Tsrt/TOTAL);
	fprintf(comfile,"%-16.3f%-16.3f%-16.3f",Tadd/tb,Tmul/tb,Tmem/tb);
	fprintf(comfile,"%-16.3f%-16.3f%-16.3f",Tdiv/tb,Tsrt/tb,TOTAL/tb);
	fprintf(comfile,"%-16.0f%-16.0f%-16.0f\n",(float)L,(float)M,(float)N);
    }
}
