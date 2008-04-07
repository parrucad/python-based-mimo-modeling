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

/* $Id: memoryalloc.c,v 1.7 2004/08/16 17:25:26 miguel Exp $ */

/* Miguel Bazdresch */

/* MIMO Simulation */

/* Memory allocation functions */
/* These two functions allocate and deallocate all needed */
/* memory. */

int allocate_mem( int L_max )
{
    unsigned int *nr_ivector(int nl, int nh);
    unsigned int **nr_imatrix(int nrl, int nrh, int ncl, int nch);
    double       **nr_matrix(int nrl, int nrh, int ncl, int nch);
    double       *nr_vector(int nl, int nh);

    int          i;
    
    /* source memory */
    if ((pSourceBits = nr_ivector(1, iNumberInfoBitsFrame)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 1;
    }
    /* training sequence memory */
    if ((pTrainSequence = nr_ivector(1, iNumberTrainBitsFrame)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 2;
    }
    /* transmit matrix */
    if ((pTransmitMatrix = nr_imatrix(1, M, 1, (Lt+L_max)*iBitsPerSymbol)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 3;
    }
    /* modulator memory */
    if ((pModulated=nr_matrix(1, M2, 1, Lt+L_max)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 4;
    }
    /* received memory */
    if ((pReceived=nr_matrix(1, N2, 1, Lt+L_max)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 6;
    }
    if ((noisematrix=nr_matrix(1, N2, 1, Lt+L_max)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 6;
    }
    if ((testvector=nr_vector(1, N2)) == NULL) {
	printf("  Memory: Error allocating memory!\n");
	return 6;
    }
    /* channel matrix */
    if ((pH=nr_matrix(1, N2, 1, M2)) == NULL) {
        printf("  Memory: Error allocating memory!\n");
        return 8;
    }
    /* estimated bits */
    if ((pEstimatedBits=nr_ivector(1, iNumberInfoBitsFrame)) == NULL) {
        printf("  Memory: Error allocating memory!");
        return 12;
    }
   /* VBLAST receiver memory */
    if( iReceiverType >= 1 && iReceiverType <= 6 ) {
	/* VBLAST */
	/* used by svdcmp */
	if ((rv1=nr_vector(1, M2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 26;
        }
        if ((pV=nr_matrix(1, M2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 16;
        }
        if ((pW=nr_vector(1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 17;
        }
        if ((r=nr_vector(1, N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 10;
        }
        if ((k=nr_ivector(1, M)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 13;
        }
        if ((pT=nr_matrix(1, M2, 1, N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 14;
        }
        if ((pTT=nr_matrix(1, M2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 15;
        }
        if ((pTTinv=nr_matrix(1, M2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 15;
        }
        if ((pHWork=nr_matrix(1, N2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 21;
        }
        if ((pHTemp=nr_matrix(1, N2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 22;
        }
        if ((pZ=nr_matrix(1, N2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 22;
        }
	if ((Q=nr_matrix(1, N2, 1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((R=nr_matrix(1, N2, 1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((a=nr_vector(1, M2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
	origor = nr_ivector( 1, M );
	if ( origor == NULL ) {
	    printf( "  Memory: Error allocating memory!" );
	    return 23;
	}
	rc = nr_ivector( 1, M );
	if ( rc == NULL ) {
	    printf( "  Memory: Error allocating memory!" );
	    return 23;
	}
	permut = nr_ivector( 1, M );
	if ( permut == NULL ) {
	    printf( "  Memory: Error allocating memory!" );
	    return 23;
	}
	u3_= nr_vector( 1, M2 );
	if (u3_ == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
	if ((b=nr_vector(1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
	if ((xls=nr_vector(1, M2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
        if ((wg=nr_vector(1, N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 25;
        }
        if ((col_inv=nr_vector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 46;
        }
        if ((storedpinvs=(double ***) malloc((size_t) \
			 ((M+1)*sizeof(double ***)))) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 30;
	}
	for(i=1;i<=M;i++) {
            if ((storedpinvs[i]=nr_matrix(1, M2, 1, N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 30;
            }
        }
#ifdef NR_LICENSED
        if ((ludcmp_vv=nr_vector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 37;
        }
        if ((ludcmp_indx=nr_ivector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 38;
        }
#endif
	/* vblast-reduced stuff */
	if ((Q1=nr_matrix(1, N2, 1, M2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
        if ((u3=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 43;
        }
        if ((x3=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 44;
        }
	if ((xls=nr_vector(1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
	if ((b=nr_vector(1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
        if ((G3=nr_matrix(1,M2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 32;
        }
        if ((G3_org=nr_matrix(1,M2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 32;
        }

    }
    else if( iReceiverType == 0 ) {
	/* ML */
	if ((b=nr_vector(1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
	if ((xls=nr_vector(1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 23;
	}
	if ((tempM=nr_matrix(1, N2, 1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((tempV=nr_vector(1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((Q=nr_matrix(1, N2, 1, M2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((R=nr_matrix(1, M2, 1, M2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((Q1=nr_matrix(1, M2, 1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
	if ((R1=nr_matrix(1, N2, 1, N2)) == NULL) {
	    printf("  Memory: Error allocating memory!");
	    return 22;
	}
        if ((G=nr_matrix(1,M2,1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 31;
        }
	if ((G_t=nr_matrix(1,N2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 31;
        }
        if ((G3_t=nr_matrix(1,M2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 32;
        }
        if ((Q_t=nr_matrix(1,N2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 32;
        }
        if ((H3=nr_matrix(1,M2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 32;
        }
        if ((bstar=nr_matrix(1,M2,1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 33;
        }
        if ((mu=nr_matrix(1,M2,1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 34;
        }
        if ((B=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 35;
        }
        if ((bt=nr_vector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 36;
        }
        if ((dist=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 39;
        }
        if ((e=nr_matrix(1,M2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 40;
        }
        if ((u=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 41;
        }
        if ((step=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 42;
        }
        if ((u3=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 43;
        }
        if ((xest=nr_vector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 44;
        }
        if ((x3=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 44;
        }
        if ((G_copy=nr_matrix(1,M2,1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 45;
        }
        if ((G_i=nr_matrix(1,N2,1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 46;
        }
        if ((foo=nr_vector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 45;
        }
        if ((col_inv=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 46;
        }
	if ((translate=nr_vector(1,N2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 47;
        }
#ifdef NR_LICENSED
        if ((ludcmp_vv=nr_vector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 37;
        }
        if ((ludcmp_indx=nr_ivector(1,M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 38;
        }
#endif
        if ((storedpinvs=(double ***) malloc((size_t) \
			 ((M2+1)*sizeof(double ***)))) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 30;
	}
	for(i=1;i<=M2;i++) {
            if ((storedpinvs[i]=nr_matrix(1, M2, 1, M2)) == NULL) {
            printf("  Memory: Error allocating memory!");
            return 39;
            }
        }

    }

    return 0;
}
