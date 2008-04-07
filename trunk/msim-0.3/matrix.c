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

/* $Id: matrix.c,v 1.9 2004/08/16 21:52:48 miguel Exp $ */

/* NOTE:
 * all operations are on complex, double matrices unless
 * otherwise noted.
 * the dimensions refer to the complex matrix, not the matrix
 * as stored in memory
 * */

/* -- copying operations -- */

/* dest = source */
inline void matrix_copy( double **dest, double **source, \
                 int rows, int cols )
{
    int i, j;

    for( i=1;i<=2*rows;i++ )
        for( j=1;j<=cols;j++ )
            dest[i][j] = source[i][j];
}

/* dest[:,i] = source[:,j] */
inline void matrix_col_copy( double **dest, double **source, \
	int rows, int i, int j ) {

    int index;

    for( index=1; index<=2*rows; index++ )
	dest[index][i] = source[index][j];
}

/* -- zeroing operations -- */

/* zero col of a matrix */
inline void matrix_zerocol(double **dest, int col_to_zero, int rows)
{
    int i;

    for(i=1;i<=2*rows;i++)
        dest[i][col_to_zero] = 0;
}

/* zero: zeroes matrix entries less than ZTH */
#define ZTH 1e-6
inline void matrix_zero_th( double **Z, int row, int col ) {

    int i, j;

    for( i=1; i<=row; i++ )
	for( j=1; j<=col; j++ )
	    if( fabs(Z[i][j]) < ZTH )
		Z[i][j] = 0;
}
#undef ZTH

/* zero complex matrix R */
void matrix_zero( double **R, int N, int M ) {

    int i, j;
    
    for( i=1; i<=2*N; i++ )
	for( j=1; j<=M; j++ )
	    R[i][j] = 0;
}

/* dest = 0 ; dest: integer real vector*/
inline void ivector_zero(int *dest, int rows)
{
    int i;

    for(i=1;i<=rows;i++)
        dest[i]=0;
}

/* -- multiplications -- */

/* dest = M1 * M2, real matrices */
void rmult_mat_mat(double **dest, double **M1, double **M2, int m, int p, int n)
{
    /* m = rows(M1), r = cols(M1), n=rows(M2) */
    
    int i, j, k;

    for(i=1;i<=m;i++) {
	double *desti=dest[i],*M1i=M1[i];
	for(j=1;j<=n;j++) {
	    desti[j]=M1i[1]*M2[1][j];
	}
    }
    for(i=1;i<=m;i++) {
	double *desti=dest[i],*M1i=M1[i];
        for(j=1;j<=n;j++) {
	    double destij=desti[j];
            for(k=2;k<=p;k++) {
                destij+=M1i[k]*M2[k][j];
	    }
	    desti[j]=destij;
	}
    }
}

/* dest = M1 * M2', real matrices */
void rmult_mat_mat_tran(double **dest, double **M1, double **M2, int m, int r, int n)
{
    /* m = rows(M1), r = cols(M1), n=rows(M2) */
    
    int i, j, k;

    for(i=1;i<=m;i++) {
	double *desti=dest[i],*M1i=M1[i];
	for(j=1;j<=n;j++) {
	    desti[j]=M1i[1]*M2[j][1];
	}
    }
    for(i=1;i<=m;i++) {
	double *desti=dest[i],*M1i=M1[i];
        for(j=1;j<=n;j++) {
	    double *M2j=M2[j],destij=desti[j];
            for(k=2;k<=r;k++) {
                destij+=M1i[k]*M2j[k];
	    }
	    desti[j]=destij;
	}
    }
}

/* dest = M1*M2' ; dest, M1, M2: complex matrices */
void mult_mat_mat_tran(double **dest, double **M1, double **M2, int m, int r, int n)
{
    /* m = rows(M1), r = cols(M1), n=rows(M2) */
    
    int i, j, k;

    for(i=1;i<=m;i++) {
	double *M1i=M1[i],*M1ipm=M1[i+m];
	double *desti=dest[i],*destipm=dest[i+m];;
	for(j=1;j<=n;j++) {
	    desti[j]=M1i[1]*M2[j][1]+M1ipm[1]*M2[j+n][1];
	    destipm[j]=-M1i[1]*M2[j+n][1]+M1ipm[1]*M2[j][1];
	}
    }
    for(i=1;i<=m;i++) {
	double *M1i=M1[i],*M1ipm=M1[i+m];
	double *desti=dest[i],*destipm=dest[i+m];;
	for(j=1;j<=n;j++) {
	    double *M2j=M2[j],*M2jpn=M2[j+n];
	    double destij=desti[j],destipmj=destipm[j];
	    for(k=2;k<=r;k++) {
		/* dest[i][j]+=M1[i][k]*M2[j][k]; */
		double M1ik=M1i[k],M1ipmk=M1ipm[k];
		double M2jk=M2j[k],M2jpnk=M2jpn[k];
		destij+=M1ik*M2jk+M1ipmk*M2jpnk;
		destipmj+=-M1ik*M2jpnk+M1ipmk*M2jk;
	    }
	    desti[j]=destij;
	    destipm[j]=destipmj;
	}
    }
}

/* dest = M1'*M2; dest, M1, M2: complex matrices */
void mult_mat_tran_mat(double **dest, double **M1, double **M2, int m, int r, int n)
{
    /* m = rows(M1), r = cols(M1), n = cols(M2) */

    int i, j, k;
    double *desti,*destipr;

    for(i=1;i<=r;i++) {
	desti=dest[i];destipr=dest[i+r];
	for(j=1;j<=n;j++) {
	    desti[j]=M1[1][i]*M2[1][j]+M1[1+m][i]*M2[1+m][j];
	    destipr[j]=M1[1][i]*M2[1+m][j]-M1[1+m][i]*M2[1][j];
	}
    }
    for(i=1;i<=r;i++) {
	desti=dest[i];destipr=dest[i+r];
	for(j=1;j<=n;j++) {
	    double destij=desti[j],destiprj=destipr[j];
	    for(k=2;k<=m;k++) {
		/* dest[i][j]+=M1[k][i]*M2[k][j]; */
		destij+=M1[k][i]*M2[k][j]+M1[k+m][i]*M2[k+m][j];
		destiprj+=M1[k][i]*M2[k+m][j]-M1[k+m][i]*M2[k][j];
	    }
	    desti[j]=destij;
	    destipr[j]=destiprj;
	}
    }
}

/* C = A.*B (cross product) ; C, A, B: complex vectors */
inline void vector_mult(double *C, double *A, double *B, int rows )
{
    int i;

    for( i=1; i<=rows; i++ )
        cmul( &C[i], &C[i+rows], A[i], A[i+rows], B[i], B[i+rows] );
}

/* -- operations on rows,cols of a matrix -- */

/* find the magnitude of column j of matrix Q; Q has N rows */
inline double col_magnitude( double **Q, int N, int j ) {

    double answer;
    int    i;

    answer = 0;
    for( i=1; i<=2*N; i++ )
	answer += Q[i][j] * Q[i][j];

    return answer;
}

/* exchange columns i and k of complex matrix A, which has row rows */
inline void col_exchange( double **A, int row, int i, int k ) {

    int    index;
    double temp;

    for( index=1; index<=2*row; index++ ) {
	temp = A[index][i];
	A[index][i] = A[index][k];
	A[index][k] = temp;
    }
}

/* divide column j of A by div; A has N rows */
void col_divide( double **A, int N, int j, double dreal, double dimag ) {

    int    index;
    double cd, temp_r, temp_i;

    cd = dreal*dreal + dimag*dimag;

    for( index=1; index<=N; index++ ) {
	temp_r = A[index][j] * dreal + A[index+N][j] * dimag;
	temp_i = A[index+N][j] * dreal - A[index][j] * dimag;
	    
	A[index][j] = temp_r / cd;
	A[index+N][j] = temp_i / cd;
    }

}

/* multiply hermitian of column i of A by column j of A; A has N rows
 * and return real part of multiplication */
inline double col_col_mult_real( double **A, int N, int i, int j) {
    
    int    row;
    double answer;

    answer = 0;
    for( row=1; row<=N; row++ )
	answer += A[row][i] * A[row][j] + A[row+N][i] * A[row+N][j];

    return answer;

}

/* multiply hermitian of column i of A by column j of A; A has N rows
 * and return imag part of multiplication */
inline double col_col_mult_imag( double **A, int N, int i, int j) {
    
    int    row;
    double answer;

    answer = 0;
    for( row=1; row<=N; row++ )
	answer += A[row][i] * A[row+N][j] - A[row+N][i] * A[row][j];

    return answer;

}

/* -- matrix decompositions and other operations -- */

/* tQR
 * computes the thin QR decomp. of real matrix QR
 * NOTE: This function destroys QR!!
 */ 
void tQR(int m, int n, double **QR, double **Qd, double **Rd)
{
    int k, j, i;
    double temp;

    /* Dim: Q is mxn
     *      R is nxn */
    
    /* main tQR loop */
    for( k=1; k<=n; k++ ) {
	
	/* norm(G[1:m,k]) */
        temp = 0.0;
        for( i=1; i<=m; i++ )
            temp += QR[i][k] * QR[i][k];
        temp = pow( temp, 0.5 );
        Rd[k][k] = temp;
	
	xx {xh[mem]+=(m+1);xh[srt]++;xh[add]+=m;xh[mul]+=m;}
	
	/* rest of algorithm */
        for( i=1; i<=m; i++ )
            Qd[i][k] = QR[i][k] / temp;
	
	xx {xh[mem]+=2*m;xh[div]+=m;}
	
        for( j=k+1; j<=n; j++ ) {
            temp = 0.0;
            for( i=1; i<=m; i++ )
                temp += Qd[i][k] * QR[i][j];
            Rd[k][j] = temp;
	    
	    xx {xh[mem]+=(2*m+1);xh[mul]+=m;xh[add]+=m;}
	    
            for( i=1; i<=m; i++ )
                QR[i][j] -= Qd[i][k] * temp;
	    
	    xx {xh[mem]+=3*m;xh[add]+=m;xh[mul]+=m;}
        }
    }
}

/* tQR_trans
 * computes the thin QR decomp. of real matrix QR' */
void tQR_trans(int M, int N, double **QR, double **Q, double **R)
{
    int k, j, i;
    double temp;

    /* clear Q and R
     * Dim: QR' is NxM
     *      Q is NxM
     *      R is MxM
     */
    
    /* main tQR loop */
    for( k=1; k<=N; k++ ) {
	
	/* norm(G[1:m,k]) */
        temp = 0.0;
        for( i=1; i<=M; i++)
            temp += QR[k][i] * QR[k][i];
        temp = pow( temp,0.5 );
        R[k][k] = temp;
	
	xx {xh[mem]+=(2*M+1);xh[srt]++;xh[add]+=M;xh[mul]+=M;}
	
	/* rest of algorithm */
        for( i=1; i<=M; i++)
            Q[i][k] = QR[k][i] / temp;
	
	xx {xh[mem]+=2*M;xh[div]+=M;}
	
        for( j=k+1; j<=N; j++ ) {
            temp = 0.0;
            for( i=1; i<=M; i++ )
                temp += Q[i][k] * QR[j][i];
            R[k][j] = temp;
	    
	    xx {xh[mem]+=(2*M+1);xh[mul]+=M;xh[add]+=M;}
	    
            for( i=1; i<=M; i++)
                QR[j][i] -= Q[i][k] * temp;
	    
	    xx {xh[mem]+=2*M;xh[add]+=M;xh[mul]+=M;}
        }
    }
}

/* complex thinqr
 * computes the thin QR decomp. of complex matrix QR
 * destroys original matrix QR */
void ctQR(int rows, int cols, double **QR, double **Q, double **R)
{
    int k, j, i;
    double tempr, tempi;
    long long int *yy;

    if( iReceiverType == 0 ) yy=xp;
    else yy=xh;

    /* Dim: Q is rowsxcols (complex) rows>cols
     *      R is colsxcols (complex) */
    
    for( i=1; i<=cols; i++ )
	for( j=1; j<=i; j++ ) {
	    R[i][j] = 0.0;
	    R[i+cols][j] = 0.0;
	    
	    xx {yy[mem_mz]+=2;}
	}
    
    /* main tQR loop */
    for( k=1; k<=cols; k++ ) {

	/* find norm(G[1:m,k]) */
	tempr = 0.0;
	for( i=1; i<=rows; i++) {
	    tempr += QR[i][k]*QR[i][k] + QR[i+rows][k]*QR[i+rows][k];
	}
        tempr = sqrt( tempr );
        R[k][k] = tempr;
	
	xx {yy[mem]+=(1+2*rows);yy[srt]++;yy[add]+=2*rows;yy[mul]+=2*rows;}
	
        for( i=1; i<=rows; i++ ) {
            Q[i][k] = QR[i][k] / tempr;
	    Q[i+rows][k] = QR[i+rows][k] / tempr;
	}
	
	xx {yy[mem]+=4*rows;yy[div]+=2*rows;}
	
        for( j=(k+1); j<=cols; j++) {
	    
	    double Qik,QRij,Qipmk,QRipmj;
	    
            tempr = 0.0;
	    tempi = 0.0;
	    
            for( i=1; i<=rows; i++) {
		Qik = Q[i][k];
		QRij = QR[i][j];
		Qipmk = Q[i+rows][k];
		QRipmj = QR[i+rows][j];
		tempr += ( Qik*QRij + Qipmk*QRipmj );
	        tempi += ( Qik*QRipmj - Qipmk*QRij );
	    }
            R[k][j] = tempr;
	    R[k+cols][j] = tempi;
            for( i=1; i<=rows; i++) {
		Qik = Q[i][k];
		Qipmk = Q[i+rows][k];
                QR[i][j] -= ( Qik*tempr - Qipmk*tempi );
		QR[i+rows][j] -= ( Qik*tempi + Qipmk*tempr );
	    }
	    
	    xx {yy[mem]+=(10*rows+2);yy[add]+=8*rows;yy[mul]+=8*rows;}
	    
	}
    }
}

/* solve: finds xls such that UT * xls = b
 * UT is a real upper-triangular matrix */
void solve(double **UT, double *xls, double *b, int n)
{
    int i,j;
    double temp;
    
    if( fabs( UT[n][n] ) > 1e-6 )
	xls[n] = b[n] / UT[n][n];
    else
	xls[n] = 0;
    
    for( i=n-1; i>=1; i-- ) {
        temp = 0.0;
        for( j=i+1;j<=n;j++ )
            temp += UT[i][j] * xls[j];
	if( fabs( UT[i][i] ) > 1e-6 )
	    xls[i] = ( 1 / UT[i][i] ) * ( b[i] - temp );
	else
	    xls[i]=0;
    }
}

/* finds the inverse of non-singular complex upper-triangular
 * matrix A and stores it in T 
 * A is NxN */
void inverse_ut_complex(int N, double **A, double **T) {

    int col, i, j;
    double r, im;
    
    for( col=N; col>=1; col-- ) {
	
	/* T[col][col] = 1/A[N][N]; */
	crep( A[col][col], A[col+N][col], &T[col][col], &T[col+N][col] );
	
	xx {xh[mem]+=4;xh[mul]+=8;xh[add]+=3;xh[div]++;}

	for( i=col-1; i>=1; i-- ) {
	    T[i][col] = 0;
	    T[i+N][col] = 0;
	    for( j=i+1; j<=col; j++ ) {
		/* T[i][col] -= A[i][j] * T[j][col]; */

		cmul( &r, &im, A[i][j], A[i+N][j], T[j][col], T[j+N][col] );
		T[i][col] -= r;
		T[i+N][col] -= im;

		xx {xh[mem]+=8;xh[mul]+=4;xh[add]+=4;}
	    }
	    
	    /* T[i][col] /= A[i][i]; */
	    cdiv( T[i][col], T[i+N][col], A[i][i], A[i+N][i], \
		  &T[i][col], &T[i+N][col] );
		
	    xx {xh[mem]+=4;xh[mul]+=8;xh[add]+=3;xh[div]++;}

	}
    }
}

/* solve_UT: finds xls such that UT*xls=b
 * UT is a real upper-triangular matrix */
void solve_UT(double **UT, double *xls, double *b, int n)
{
    int i,j;
    double temp;
    
    if( fabs( UT[n][n] ) > 1e-6 )
	xls[n] = b[n] / UT[n][n];
    else
	xls[n] = 0;
    
    for( i=n-1; i>=1; i-- ) {
        temp = 0.0;
        for( j=i+1;j<=n;j++ )
            temp += UT[j][i] * xls[j];
	if( fabs( UT[i][i] ) > 1e-6 )
	    xls[i] = ( 1 / UT[i][i] ) * ( b[i] - temp );
	else
	    xls[i]=0;
    }
}

/* inverse_UT: finds the inverse of an upper-triangular real matrix
 * of size rowsxrows */
void inverse_UT( double **A, double **Ai, int rows)
{
    void solve_UT(double **UT, double *xls, double *b, int n);

    int i, j;

    for( i=1; i<=rows; i++ ) {
	for( j=1; j<=rows; j++ )
	    b[j] = 0;
	b[i] = 1;
	
	solve_UT( A, xls, b, rows );

	for( j=1; j<=rows; j++ )
	    Ai[j][i] = xls[j];
    }
}
	
/* inverse_LT_trans: finds the inverse of the transpose of a lower-triangular
 * real matrix A and stores it in Ai. Matrix A is not affected by this
 * function. Matrix A is rowsxrows */
void inverse_LT_trans( double **A, double **Ai, int rows )
{
    void solve_LT_trans(double **LT, double *xls, double *b, int n);
    
    int i, j;

    for( i=1; i<=rows; i++ ) {
	for( j=1; j<=rows; j++ )
	    b[j] = 0;
	b[i] = 1;

	solve_LT_trans( A, xls, b, rows );

	for( j=1; j<=rows; j++ )
	    Ai[j][i] = xls[j];
    }
}

/* solve_LT_trans: finds xls such that UT*xls=b
 * LT is a real lower-triangular matrix */
void solve_LT_trans(double **LT, double *xls, double *b, int n)
{
    int i,j;
    double temp;

    if( fabs( LT[1][1] ) > 1e-6 )
	xls[1] = b[1] / LT[1][1];
    else
	xls[1] = 0;

    for( i=2; i<=n; i++ ) {
	temp = 0.0;
	for( j=1; j<i; j++ )
	    temp += LT[j][i] * xls[j];
	if( fabs( LT[i][i] ) > 1e-6 )
	    xls[i] = ( b[i] - temp ) / LT[i][i];
	else
	    xls[i] = 0;
    }
}

/* solve_withcounters: finds xls such that UT*xls=b
 * includes complexity counters
 * UT is a real upper-triangular matrix */
void solve_withcounters(double **UT, double *xls, double *b, int n)
{
    int i,j;
    double temp;
    
    if( fabs( UT[n][n] ) > 1e-6 ) {
	xls[n] = b[n] / UT[n][n];

	xx {xh[mem]+=3;xh[div]++;}
    }
    else {
	xls[n] = 0;
	
        xx {xh[mem]++;}
    }

    xx {xh[mem]++;xh[add]++;} 
    
    for( i=n-1; i>=1; i-- ) {
        temp = 0.0;
        for( j=i+1; j<=n; j++ )
            temp += UT[i][j] * xls[j];

	xx {xh[mem]+=2*n;xh[mul]+=n;xh[add]+=n;}
	
	if( fabs( UT[i][i] ) > 1e-6 ) {
	    xls[i] = ( 1 / UT[i][i] ) * ( b[i] - temp );

	    xx {xh[mem]+=3;xh[add]++;xh[div]++;}
	}
	else {
	    xls[i] = 0;

	    xx {xh[mem]++;}
	}

	xx {xh[mem]++;xh[add]++;}
    }
}

/* finds xls such that R*xls=b for complex numbers */
void csolve1(double **R, double *xls, double *b, int Rrows, int Rcols)
{
    int    i, j, m, n;
    double tempr, tempi, cd2, Rnn, Rnpmn, bn, bnpm;
    unsigned long long *yy;

    if( (iReceiverType == 3) || (iReceiverType == 5) ) yy=xh;
    else if ( iReceiverType == 0 ) yy=xp;
    else yy=xq;
    
    m=Rrows;
    n=Rcols;
    
    Rnn = R[n][n];
    Rnpmn = R[n+m][n];
    if( fabs(Rnn) < 1e-6 ) { R[n][n] = 0; Rnn = 0; xx {yy[mem]++;} }
    if( fabs(Rnpmn) < 1e-6 ) { R[n+m][n] = 0; Rnpmn = 0; xx {yy[mem]++;} }
    
    xx {yy[mem]+=2;yy[add]+=4;}

    bn = b[n];
    bnpm = b[n+m];
    cd2 = Rnn*Rnn + Rnpmn*Rnpmn;
    
    xx {yy[mem]+=2;yy[mul]+=2;yy[add]++;}

    if( cd2 != 0 ) {
	/* xls[N]=b[N]/R[N][N]; */
	xls[n] = (Rnn*bn + Rnpmn*bnpm) / cd2;
	xls[2*n] = (Rnn*bnpm - Rnpmn*bn) / cd2;
	
	xx {yy[mem]+=2;yy[mul]+=4;yy[add]+=2;yy[div]+=2;}
	
    } else {
	xls[n] = 0;
	xls[2*n] = 0;
	
	xx {yy[mem]+=2;}
	
    }
    
    for( i=n-1; i>=1; i--) {
	
	double Rii,Ripmi;
	
	Rii = R[i][i];
	Ripmi = R[i+m][i];
	if( fabs(Rii) < 1e-6 ) { Rii = R[i][i] = 0; xx {yy[mem]++;} }
	if( fabs(Ripmi) < 1e-6) { Ripmi = R[i+m][i] = 0; xx {yy[mem]++;} }
	
	xx {yy[mem]+=2;yy[add]+=4;}
	
	if( Rii == 0 && Ripmi == 0 ) {
	    xls[i] = 0;
	    xls[i+n] = 0;
	    
	    xx {yy[mem]+=2;}
	    
	} else {
	    tempr = 0.0;
	    tempi = 0.0;
	    for( j=i+1; j<=n; j++ ) {
		tempr += R[i][j]*xls[j] - R[i+m][j]*xls[j+n];
		tempi += R[i][j]*xls[j+n] + R[i+m][j]*xls[j];
		
		xx {yy[mem]+=4;yy[mul]+=4;yy[add]+=4;}
		
	    }
	    /* xls[i]=(1/R[i][i])*(b[i]-temp); */
	    tempr = b[i] - tempr;
	    tempi = b[i+m] - tempi;
	    cd2 = Rii*Rii + Ripmi*Ripmi;
	    xls[i] = (tempr*Rii + tempi*Ripmi) / cd2;
	    xls[i+n] = (tempi*Rii - tempr*Ripmi) / cd2;
	    
	    xx {yy[mem]+=4;yy[mul]+=6;yy[div]+=2;yy[add]+=5;}
	    
	}
    }
}

/* finds xls such that R*xls=b for complex numbers, but does slicing in
 * each intermediate step */
void csolve_slice(double **R, double *xls, double *b, int Rrows, int Rcols)
{
    int i,j,m,n;
    double tempr, tempi, cd2, Rnn, Rnpmn, bn, bnpm;
    unsigned long long *yy;

    yy=xv;
    
    m=Rrows;
    n=Rcols;
    
    Rnn = R[n][n];
    Rnpmn = R[n+m][n];
    if( fabs(Rnn) < 1e-6 ) { R[n][n] = 0; Rnn = 0; xx {yy[mem]++;} }
    if( fabs(Rnpmn) < 1e-6 ) { R[n+m][n] = 0; Rnpmn = 0; xx {yy[mem]++;} }
    
    xx {yy[mem]+=2;yy[add]+=4;}

    bn = b[n];
    bnpm = b[n+m];
    cd2 = Rnn*Rnn + Rnpmn*Rnpmn;
    
    xx {yy[mem]+=2;yy[mul]+=2;yy[add]++;}

    if( cd2 != 0 ) {
	/* xls[N]=b[N]/R[N][N]; */
	xls[n] = slice( (Rnn*bn + Rnpmn*bnpm) / cd2 );
	xls[2*n] = slice( (Rnn*bnpm - Rnpmn*bn) / cd2 );
	
	xx {yy[mem]+=2;yy[mul]+=4;yy[add]+=2;yy[div]+=2;}
	
    } else {
	xls[n] = 0;
	xls[2*n] = 0;
	
	xx {yy[mem]+=2;}
	
    }
    
    for( i=n-1; i>=1; i--) {
	
	double Rii,Ripmi;
	
	Rii = R[i][i];
	Ripmi = R[i+m][i];
	if( fabs(Rii) < 1e-6 ) { Rii = R[i][i] = 0; xx {yy[mem]++;} }
	if( fabs(Ripmi) < 1e-6) { Ripmi = R[i+m][i] = 0; xx {yy[mem]++;} }
	
	xx {yy[mem]+=2;yy[add]+=4;}
	
	if( Rii == 0 && Ripmi == 0 ) {
	    xls[i] = 0;
	    xls[i+n] = 0;
	    
	    xx {yy[mem]+=2;}
	    
	} else {
	    tempr = 0.0;
	    tempi = 0.0;
	    for( j=i+1; j<=n; j++ ) {
		tempr += R[i][j]*xls[j] - R[i+m][j]*xls[j+n];
		tempi += R[i][j]*xls[j+n] + R[i+m][j]*xls[j];
		
		xx {yy[mem]+=4;yy[mul]+=4;yy[add]+=4;}
		
	    }
	    tempr = b[i] - tempr;
	    tempi = b[i+m] - tempi;
	    cd2 = Rii*Rii + Ripmi*Ripmi;
	    xls[i] = slice( (tempr*Rii + tempi*Ripmi) / cd2 );
	    xls[i+n] = slice( (tempi*Rii - tempr*Ripmi) / cd2 );
	    
	    xx {yy[mem]+=4;yy[mul]+=6;yy[div]+=2;yy[add]+=8;}
	    
	}
    }
}

/* find givens transformation of complex fr, fi, gr, gi */
void cgivens(double fr, double fi, double gr, double gi, \
	double *c, double *sr, double *si)
{
    double f2,g2,fg2,d1,temp;

    f2 = fr*fr + fi*fi;
    g2 = gr*gr + gi*gi;
    fg2 = f2 + g2;
    d1 = 1 / sqrt( f2*fg2 );
    *c = f2 * d1;
    fg2 = fg2 * d1;
    *sr = fr * d1;
    *si = fi * d1;
    /* s=conj(g)*s */
    temp = (*sr)*gr + (*si)*gi;
    *si = -(*sr)*gi + (*si)*gr;
    *sr = temp;
    
    xx {xq[mem]+=6;xq[mul]+=13;xq[div]++;xq[srt]++;xq[add]+=6;}
    
}

/* QR decomposition */
void QR( double **Qi, double **Ri, double **Z, int rows, int cols )
{
    int    i, j, l;
    double c, si, sr, tau1R, tau1I, tau2R, tau2I;

    /* Q is rows*rows, R is rows*cols */
    /* Qi=I */
    for( i=1; i<=rows; i++ ) {
	/* real part of Qi is equal to I */
	Qi[i][i] = 1;
	for( j=1; j<i; j++ )
	    Qi[i][j] = 0;
	for( j=i+1; j<=rows; j++ )
	    Qi[i][j] = 0;
	/* imaginary part of Qi is equal to 0 */
	for( j=1; j<=rows; j++ )
	    Qi[i+rows][j] = 0;
    }
    
    xx {xq[mem_mz]+=rows*rows;}
    
    for( j=1; j<=cols; j++ ) {
	for(i=rows;i>=( j+1); i--) {
	    cgivens( Z[i-1][j], Z[i-1+rows][j], \
		     Z[i][j], Z[i+rows][j], &c, &sr, &si );
	    
	    xx {xq[mem]+=4;}
	    
	    for( l=1; l<=cols; l++ ) {
		tau1R = Z[i-1][l];
		tau1I = Z[i-1+rows][l];
		tau2R = Z[i][l];
		tau2I = Z[i+rows][l];
		
		/* real part of Z */
		Z[i-1][l] =  c*tau1R + sr*tau2R - si*tau2I;
	        Z[i][l]   = -sr*tau1R - si*tau1I + c*tau2R;
		
		/* imag part of Z */
		Z[i-1+rows][l] =  c*tau1I + sr*tau2I + si*tau2R;
		Z[i+rows][l]   = -sr*tau1I + si*tau1R + c*tau2I;
	    }
	    
	    for( l=1; l<=rows; l++ ) {
		tau1R = Qi[l][i-1];
		tau1I = Qi[l+rows][i-1];
		tau2R = Qi[l][i];
		tau2I = Qi[l+rows][i];
		
		/* real part of Qi */
		Qi[l][i-1] =  c*tau1R + sr*tau2R + si*tau2I;
		Qi[l][i]   = -sr*tau1R + si*tau1I + c*tau2R;
		
		/* imag part of Qi */
		Qi[l+rows][i-1] =  c*tau1I + sr*tau2I - si*tau2R;
		Qi[l+rows][i]   = -sr*tau1I - si*tau1R + c*tau2I;
	    }
	    
	    xx {xq[mem]+=8*(cols+rows);xq[mul]+=12*(cols+rows);xq[add]+=10*(cols+rows);}
	}
    }
    /* copy Z to R */
    for(i=1;i<=rows;i++)
	for(j=1;j<=cols;j++) {
	    Ri[i][j]=Z[i][j];
	    Ri[i+rows][j]=Z[i+rows][j];
	}
    xx {xq[mem_mc]+=rows*cols*4;}
}

void updateqr(double **Q, double **R, int m, int n, int rc)
{
    int    index, j;
    double c, si, sr; /* givens rotation */
    double tau1R, tau1I, tau2R, tau2I;

    for( index=rc; index<=n; index++ ) {
	cgivens( R[index][index], R[index+m][index], \
		 R[index+1][index], R[index+1+m][index], &c, &sr, &si);
	
	xx {xq[mem]+=4;}
	
	/* apply the rotation to R */
	/* only rows index and index+1 are affected */
	for( j=1; j<=n; j++ ) {
	    tau1R = R[index][j];
	    tau1I = R[index+m][j];
	    tau2R = R[index+1][j];
	    tau2I = R[index+1+m][j];
	    
	    /* real part of R */
	    R[index][j]=c*tau1R+sr*tau2R-si*tau2I;
	    R[index+1][j]=-sr*tau1R-si*tau1I+c*tau2R;
	    
	    /* imag part of R */
	    R[index+m][j]=c*tau1I+sr*tau2I+si*tau2R;
	    R[index+1+m][j]=-sr*tau1I+si*tau1R+c*tau2I;

	}
	/* apply the rotation to Q */
	/* only columns index and index+1 are affected */
	for( j=1; j<=m; j++ ) {
	    tau1R = Q[j][index];
	    tau1I = Q[j+m][index];
	    tau2R = Q[j][index+1];
	    tau2I = Q[j+m][index+1];
	    
	    /* real part of Q */
	    Q[j][index]   =  c*tau1R + sr*tau2R + si*tau2I;
	    Q[j][index+1] = -sr*tau1R + si*tau1I + c*tau2R;
	    
	    /* imag part of Q */
	    Q[j+m][index]   =  c*tau1I + sr*tau2I - si*tau2R;
	    Q[j+m][index+1] = -sr*tau1I - si*tau1R + c*tau2I;
	}
	
	xx {xq[mem]+=8*(m+n);xq[mul]+=12*(n+m);xq[add]+=10*(n+m);}
    }
}

/*
 * inverse_tqr calculates the (pseudo)inverse of a matrix A and stores
 * it in Ai. A is mxn, Ai is nxm. It uses the thin QR decomposition
 * to find the inverse.
 * Note: This function destroys the original matrix.
 * 
 */
void inverse_tqr(double **A, double **Ai, int m, int n) {

    int row, col;
    
    /* function prototpypes */
    void tQR( int m, int n, double **A, double **B, double **C );
    void solve( double **A, double *xls, double *b, int n );
    
    for( row=1; row<=n; row++ )
	for( col=1; col<=row; col++ )
	    R[row][col] = 0;
    
    /* Q = mxn, R = nxn */
    tQR( m, n, A, Q, R ); 
	
    for( row=1; row<=m; row++ ) {
	
	/* for( col=1; col<=n; col++ ) */
	    /* b[col] = Q[row][col]; */

	solve( R, xls, Q[row], n );

	for( col=1; col<=n; col++ )
	    Ai[col][row] = xls[col];

    }
}

/*
 * inverse_tqr_trans calculates the (pseudo)inverse of matrix A' and stores
 * it in Ai. A is mxn (m<n). It uses the thin QR decomposition
 * to find the inverse.
 * Note that this function destroys the original matrix.
 * 
 */
void inverse_tqr_trans(double **A, double **Ai, int m, int n) {

    void tQR_trans( int m, int n, double **A, double **Q, double **R );
    void solve_withcounters( double **A, double *xls, double *b, int n );
    
    int row,col;

    /* In thin QR, if A = mxn, then Q = mxn and R = nxn */
    
    for( row=1; row<=m; row++ )
	for( col=1; col<=row; col++ )
	    R[row][col] = 0;

    xx {xh[mem_mz]+=rint(m*m/2);}

    tQR_trans( m, n, A, Q, R );

    for( row=1; row<=m; row++ ) {

	/* for( col=1; col<=n; col++ ) */
	    /* b[col] = Q[row][col]; */

	solve_withcounters( R, xls, Q[row], n );

	for( col=1; col<=n; col++ )
	    Ai[col][row] = xls[col];

	xx {xh[mem_mc]+=n;}

    }
}
