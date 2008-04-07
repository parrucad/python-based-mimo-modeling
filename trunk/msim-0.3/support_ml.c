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

/* $Id */


/* support_ml.c
   Contains all matrix functions needed to implement ML
   Matlab code in comments
   */

/* This function used in LLL */
void red(int k, int l)
{
    int col, i;
    double q;

    /* % begin red */
    /* if abs(mu(k,k-1)) > 0.5, */
    if(debug && ((cDebug & 2048) == 2048))
	fprintf(outfile, "mu[k][l]=%1.5f  inside red\n", mu[k][l]);
    xx {xl[mem]++;xl[add]++;}
    if (fabs(mu[k][l]) > 0.5) {
	/* q = floor(0.5+mu(k,k-1)); */
	q = floor(0.5+mu[k][l]);
	xx {xl[mem]++;xl[add]+=2;}
	if(debug && ((cDebug & 2048) == 2048))
	    fprintf(outfile, "q = %1.5f\n", q);
	/* b(k,:)=b(k,:)-q*b(k-1,:); */
	for(col=1;col<=N2;col++) {
	    G[k][col]=G[k][col]-q*G[l][col];
	}
	xx {xl[mem]+=3*N2;xl[add]+=N2;xl[mul]+=N2;}
	/* mu(k,k-1)=mu(k,k-1)-q; */
	mu[k][l] = mu[k][l]-q; xx{xl[mem]+=2;xl[add]++;}
	if(debug && ((cDebug & 2048) == 2048))
	    fprintf(outfile, "mu[k][l] = %1.5f\n", mu[k][l]);
	/* for i = 1:k-2, */
	/* mu(k,i) = mu(k,i) - q*mu(k-1,i); */
	for(i=1;i<=l-1;i++) {
	    mu[k][i]=mu[k][i]-q*mu[l][i];
	    xx {xl[mem]+=3;xl[add]++;xl[mul]++;}
	    if(debug && ((cDebug & 2048) == 2048))
		fprintf(outfile, "mu[k][i] = %1.5f\n", mu[k][i]);
	}
    }
}

/* This function used in LLL */
void swap(int k, int kmax)
{
    int col, j, i;
    double t, temp, mut, Bt;

    /* temp = b(k,:); */
    /* b(k,:) = b(k-1,:); */
    /* b(k-1,:) = temp; */
    for(col=1;col<=N2;col++) {
	temp=G[k][col];
	G[k][col]=G[k-1][col];
	G[k-1][col]=temp;
    }
    xx {xl[mem]+=4*N2;}
    /* if k > 2, */
    if (k>2) {
	/* for j = 1:k-2, */
	for(j=1;j<=k-2;j++) {
	    /* temp = mu(k,j); */
	    /* mu(k,j) = mu(k-1,j); */
	    /* mu(k-1,j) = temp; */
	    temp=mu[k][j];
	    mu[k][j]=mu[k-1][j];
	    mu[k-1][j]=temp;
	    xx {xl[mem]+=4;}
	}
    }
    /* mut = mu(k,k-1); */
    mut = mu[k][k-1];
    /* Bt = B(k) + mut^2*B(k-1); */
    Bt = B[k] + mut*mut*B[k-1];
    /* mu(k,k-1) = mut*B(k-1)/Bt; */
    mu[k][k-1] = mut*B[k-1]/Bt;
    xx {xl[mem]+=4;xl[add]++;xl[mul]+=3;xl[div]++;}
    /* bt = bstar(k-1,:); */
    for(col=1;col<=N2;col++)
	bt[col]=bstar[k-1][col];
    xx {xl[mem_mc]+=2*N2;}
    /* bstar(k-1,:) = bstar(k,:) + mut*bt; */
    for(col=1;col<=N2;col++)
	bstar[k-1][col]=bstar[k][col]+mut*bt[col];
    xx {xl[mem]+=3*N2;xl[mul]+=N2;xl[add]+=N2;}
    /* bstar(k,:) = mu(k,k-1)*bstar(k,:) + (B(k)/Bt)*bt; */
    for(col=1;col<=N2;col++)
	bstar[k][col]=mu[k][k-1]*bstar[k][col]+(B[k]/Bt)*bt[col];
    xx {xl[mem]+=N2*5;xl[mul]+=2*N2;xl[add]+=N2;xl[div]+=N2;}
    /* B(k) = B(k-1)*B(k)/Bt; */
    B[k]=B[k-1]*B[k]/Bt;
    /* B(k-1) = Bt; */
    B[k-1]=Bt;
    xx {xl[mem]+=4;xl[mul]++;xl[div]++;}
    if(debug && ((cDebug & 2048) == 2048)) {
	fprintf(outfile, "k = %d  inside swap\n", k);
	fprintf(outfile, "B[k-1] = %1.5f, B[k] = %1.5f, Bt = %1.5f\n", B[k-1], B[k], Bt);
    }
    /* for i = k+1:kmax, */
    for(i=k+1;i<=kmax;i++) {
	/* t = mu(i,k); */
	t = mu[i][k];
	/* mu(i,k) = mu(i,k-1) - mut*t; */
	mu[i][k] = mu[i][k-1] - mut*t;
	/* mu(i,k-1) = t + mu(k,k-1)*mu(i,k); */
	mu[i][k-1] = t + mu[k][k-1]*mu[i][k];
	xx {xl[mem]+=5;xl[mul]+=2;xl[add]+=2;}
    }
}

void lll(void)
    /* This function calculates the LLL reduction of G */
{

    /* prototypes */
    void red(int k, int l);
    void swap(int k, int kmax);
    void debug2048(int a1);


    /*  function L = lll(b)
	%% lll
	%% This function returns the LLL reduction of input matrix b,
	%% according to Cohen p.86 and errata.
	%% The rows of B are the basis vectors to be reduced.


	disp('starting lll');

	G = b;

	% step 1 */

    /* local variables */
    int    col, k, kmax, l, j, itemp;
    int    marker=1;
    double temp;

    k = 2;
    kmax = 1;

    /* bstar(1,:) = b(1,:); */
    /* B(1)=b(1,:)*b(1,:)'; */
    B[1] = 0;
    for( col=1; col<=N2; col++ ) {
	temp = G[1][col];
	bstar[1][col] = temp;
	B[1] += temp * temp;
    }
    
    xx {xl[mem]+=(1+3*N2);xl[add]+=N2;xl[mul]+=N2;}
    
    /* n = M2 */
    /* while k <= n, */
    while( k <= M2 ) {
	/* % step 2 */
	/* if k > kmax, */
	if( k > kmax ) { /* step 2 */
	    /* kmax = k; */
	    kmax = k;
	    /* bstar(k,:) = b(k,:); */
	    for( col=1; col<=N2; col++ ) {
		bstar[k][col] = G[k][col];
	    }
	    
	    xx {xl[mem_mc]+=N2;}
	    
	    /* for j = 1:(k-1), */
	    for( j=1; j<=k-1; j++ ) {
		/* mu(k,j) = b(k,:)*bstar(j,:)'/B(j); */
		temp = 0;
		for( col=1; col<=N2; col++ )
		    temp += G[k][col] * bstar[j][col];
		
		xx {xl[mem_mm]+=2*N2;xl[add_mm]+=N2;xl[mul_mm]+=N2;}
		
		mu[k][j] = temp / B[j]; xx{ xl[mem]+=2;xl[div]++;}
		
		/* bstar(k,:)=bstar(k,:)-mu(k,j)*bstar(j,:); */
		for( col=1; col<=N2; col++ )
		    bstar[k][col] = bstar[k][col] - mu[k][j]*bstar[j][col];
		
		xx {xl[mem]+=4*N2;xl[mul]+=N2;xl[add]+=N2;}
	    }
	    /* B(k) = bstar(k)*bstar(k)'; */
	    temp = 0;
	    for( col=1; col<=N2; col++ )
		temp += bstar[k][col] * bstar[k][col];
	    B[k] = temp;
	    
	    xx {xl[mem]+=(1+2*N2);xl[mul]+=N2;xl[add]+=N2;}
	    if (debug && ((cDebug & 2048) == 2048)) {
		fprintf(outfile, "k = %d  in step 2\n", k);
		debug2048(1);
	    }
	    
	    /* if B(k) == 0, */
	    /* disp('error. vectors do not form a basis'); */
	    if (B[k] == 0) {
		printf("ERROR: vectors do not form a basis\n");
		exit(1);
	    }
	}
	/* % step 3 */
	/* while 1, */
	while( 1 ) {
	    if(debug && ((cDebug & 2048) == 2048)) {
		fprintf(outfile, "marker = %d\n", marker);
		marker++;
		fprintf(outfile, "red 1 start at beg of step 3\n");
		debug2048(1);
	    }
	    
	    red(k,k-1);
	    
	    if (debug && ((cDebug & 2048) == 2048)) {
		fprintf(outfile, "red 1 end at beg of step 3\n");
		debug2048(1);
	    }
	    
	    xx {xl[mem]+=3;xl[mul]+=2;xl[add]+=2;}
	    
	    /* if B(k) < (0.75 - mu(k,k-1)^2)*B(k-1), */
	    if( B[k] < ( (0.75 - mu[k][k-1] * mu[k][k-1] ) * B[k-1] ) ) {
		if (debug && ((cDebug & 2048) == 2048)) {
		    fprintf(outfile, "swap start\n");
		    debug2048(1);
		}
		
		swap(k, kmax);
		
		if (debug && ((cDebug & 2048)) == 2048) {
		    fprintf(outfile, "swap end\n");
		    debug2048(1);
		}
		/* k = max([2, k-1]); */
		itemp = 2;
		if( (k-1) > 2 )
		    itemp = k-1;
		k = itemp;
	    }
	    else {
		/* for l = k-2:-1:1, */
		for(l=k-2;l>=1;l--) {
		    if ((debug && (cDebug & 2048)) == 2048) {
			fprintf(outfile, "red 2 start  at end of step 3\n");
			debug2048(1);
		    }
		    
		    red(k,l);
		    
		    if ((debug && (cDebug & 2048)) == 2048) {
			fprintf(outfile, "red 2 end  at end of step 3\n");
			debug2048(1);
		    }
		}
		
		k = k + 1;
		break;
	    }
	}
    }
}

/* thinqr */
void thinqr(void)
{
    /* thinqr destroys the matrix it's operating on, so
     * a copy is created in Gtemp */

    /* Matlab code for thinqr function */
    /* function [q,r] = thinqr(a) */
    /* m = size(a,1); % number of rows of a
       n = size(a,2); % number of columns of a

       for k = 1:n,
       r(k,k) = norm(a(1:m,k));
       q(1:m,k) = a(1:m,k)/r(k,k);
       for j = (k+1):n,
       r(k,j) = q(1:m,k)'*a(1:m,j);
       a(1:m,j) = a(1:m,j)-q(1:m,k)*r(k,j);
       end
       end 
       %% Note that NORM(V) = norm(V,2) where V is a vector, and:
       %% NORM(V,P) = sum(abs(V).^P)^(1/P)
       */

    int k, j, i, m, n;

    /* This function operates on G' (N2xM2, mxn) */

    m=N2;
    n=M2;

    /* copy G into Gtemp, both are M2xN2 */
    for(i=1;i<=M2;i++)
	for(j=1;j<=N2;j++)
	    G_copy[i][j] = G[i][j];
    /* clear G3_t and Q_t. G3_t is M2xM2 (nxn), Q_t is N2xM2 (mxn) */
    for(i=1;i<=n;i++)
	for(j=1;j<=n;j++) {
	    G3_t[i][j] = 0.0;
	    Q_t[i][j] = 0.0;
	}
    for(i=n+1;i<=m;i++)
	for(j=1;j<=n;j++)
	    Q_t[i][j] = 0.0;
    /* main tqr loop */
    for(k=1;k<=n;k++) {
	/* norm(G[1:m,k]) */
	for(i=1;i<=m;i++)
	    G3_t[k][k] += SQR(G_copy[k][i]);
	G3_t[k][k] = pow(G3_t[k][k],0.5);
	/* rest of algorithm */
	for(i=1;i<=m;i++)
	    Q_t[i][k] = G_copy[k][i]/G3_t[k][k];
	for(j=(k+1);j<=n;j++) {
	    G3_t[k][j] = 0.0;
	    for(i=1;i<=m;i++)
		G3_t[k][j] += Q_t[i][k]*G_copy[j][i];
	    for(i=1;i<=m;i++)
		G_copy[j][i] -= Q_t[i][k]*G3_t[k][j];
	}
    }
}

/* decode */
/* Agrell algorithm */
void decode(void) 
{
    /* x3 is 1xM2
     * H3 is M2xM2 */

    /* function prototypes */
    double round(double);

    int n=M2;
    int k, i, j;
    double y, aux, bestdist, newdist, temp;
    /* double *dist; 1xM2
       double **e; M2xM2
       double *u; 1xM2
       double *step; 1xM2
       double *u3; 1xM2 */
    /* Original matlab code:
     * n=size(H3,1);
     bestdist=inf;
     k=n;
     dist(k)=0;
     aux=x3*H3;
     for i=1:n,
     e(k,i)=aux(i);
     end;
     u(k)=round(e(k,k));
     y=(e(k,k)-u(k))/H3(k,k);
     if (y>0), step(k)=1;
     else step(k) = -1;
     end;
     while 1,
     newdist=dist(k)+y^2;
     if newdist < bestdist,
     if k ~= 1,
     for i = 1:(k-1),
     e(k-1,i)=e(k,i)-y*H3(k,i);
     end;
     k=k-1;
     dist(k)=newdist;
     u(k)=round(e(k,k));
     y=(e(k,k)-u(k))/H3(k,k);
     if (y>0), step(k)=1;
     else step(k)=-1;
     end;
     else
     u3=u;
     bestdist=newdist;
     k=k+1;
     u(k)=u(k)+step(k);
     y=(e(k,k)-u(k))/H3(k,k);
     if step(k)>0, aux=1;
     else aux=-1;
     end
     step(k)=-step(k)-aux;
     end
     else
     if k==n,
     break;
     else
     k=k+1;
     u(k)=u(k)+step(k);
     y=(e(k,k)-u(k))/H3(k,k);
     if step(k)>0, aux=1;
     else aux=-1;
     end
     step(k)=-step(k)-aux;
     end
     end
     end */
    /* initial values */
    for( i=1; i<=M2; i++ ) {
	dist[i] = 0.0;
	u[i] = 0.0;
	u3[i] = 0.0;
	step[i] = 0.0;
	for( j=1; j<=M2; j++ )
	    e[i][j] = 0.0;
    }
    
    xx {xd[mem]+=4*M2;xd[mem_mz]+=M2*M2;}
    
    /* start */
    bestdist = 1e100;
    k = n;
    
    /* dist[k] = 0.0; already done in initial values above */
    for( i=1; i<=M2; i++ ) {
	temp = e[k][i];
	for( j=1; j<=M2; j++ )
	    temp += x3[j] * H3[j][i];
	e[k][i] = temp;
    }
    
    u[k] = round( e[k][k] );
    y = (e[k][k] - u[k]) / H3[k][k];
    
    xx {xd[mem]+=3;xd[mem_mm]+=2*M2*M2+2*M2;xd[add_mm]+=M2*M2; \
	xd[div]++;xd[mul_mm]+=M2*M2;xd[add]+=2;}
    
    if( y > 0 ) step[k] = 1;
    else step[k] = -1;
    
    xx {xd[mem]++;xd[add]++;}
    
    while( 1 ) {
	newdist = dist[k] + y*y;
	
	xx {xd[add]+=2;xd[mul]++;}
	
	if( newdist < bestdist ) {
	    if( k != 1 ) {
		for( i=1; i<=k-1; i++ ) {
		    e[k-1][i] = e[k][i] - y*H3[k][i];
		    
		    xx {xd[mem]+=3;xd[add]++;xd[mul]++;}
		}
		
		k = k-1;
		dist[k] = newdist;
		u[k] = round( e[k][k] );
		y = (e[k][k] - u[k]) / H3[k][k];
		
		xx {xd[mem]+=4;xd[add]+=2;xd[div]++;}
		
		if( y > 0 ) step[k] = 1;
		else step[k] = -1;
		
		xx {xd[mem]++;xd[add]++;}
	    }
	    else {
		for( i=1; i<=M2; i++ )
		    u3[i] = rint( u[i] );
		bestdist = newdist;
		k = k + 1;
		u[k] = u[k] + step[k];
		y = ( e[k][k] - u[k] ) / H3[k][k];
		
		xx {xd[mem_mc]+=2*M2;xd[mem]+=5;xd[div]++;xd[add]+=2;}
		
		if( step[k] > 0 ) aux = 1;
		else aux = -1;
		step[k] = -step[k] - aux;
		
		xx {xd[mem]+=2;xd[add]+=3;}
	    }
	}
	else {
	    if( k==n ) break;
	    else {
		k = k+1;
		u[k] = u[k] + step[k];
		y = ( e[k][k] - u[k] ) / H3[k][k];
		
		xx {xd[mem]+=5;xd[add]+=2;xd[div]++;}
		
		if( step[k] > 0 ) aux = 1;
		else aux = -1;
		step[k] = -step[k] - aux;
		
		xx {xd[mem]+=2;xd[add]+=3;}
	    }
	}
    }
}
