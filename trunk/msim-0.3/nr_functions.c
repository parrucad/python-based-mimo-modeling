/*
 *  This file contains the NR functions used in MIMO Simulation.
 * This code has been released into the public domain by its authors,
 * and has been taken from the book Numerical Recipes in C.
 * Besides, some of this code has been modified by Miguel Bazdresch to
 * suit the needs of the MSIM simulator.
 */

/* $Id: nr_functions.c,v 1.2 2004/07/31 00:03:44 miguel Exp $ */

/* cvector (from nrutil.c) Modified by MB*/

#define NR_END 1
#define FREE_ARG char*

unsigned char *nr_cvector(int nl, int nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
    unsigned char *v;

    v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
    if (!v) return NULL;
    return v-nl+NR_END;
}

/* vector */
double *nr_vector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) return NULL;
    return v-nl+NR_END;
}

/* ivector */
unsigned int *nr_ivector(int nl, int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    unsigned int *v;

    v=(unsigned int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned int)));
    if (!v) return NULL;
    return v-nl+NR_END;
}

/* imatrix, adapted from imatrix by MB (from nrutil.c) */

unsigned int **nr_imatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    unsigned int **m;

    /* allocate pointers to rows */
    m=(unsigned int **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned int*)));
    if (!m) return NULL;
    m += NR_END;
    m -= nrl;


    /* allocate rows and set pointers to them */
    m[nrl]=(unsigned int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned int)));
    if (!m[nrl]) return NULL;
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

/* nr_matrix, from nrutil.c */

double **nr_matrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;

    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    if (!m) return NULL;
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    if (!m[nrl]) return NULL;
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

/* free_vector */
void nr_free_vector(double *v, int nl)
/* free a double vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

/* free_ivector */
void nr_free_ivector(unsigned int *v, int nl)
/* free an int vector allocated with ivector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

/* free_matrix */

void nr_free_matrix(double **m, int nrl, int ncl)
/* free a double matrix allocated by matrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

/* free_imatrix */
void nr_free_imatrix(unsigned int **m, int nrl, int ncl)
/* free a double matrix allocated by matrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
