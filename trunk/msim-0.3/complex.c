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

/* $Id: complex.c,v 1.2 2004/08/16 17:19:12 miguel Exp $ */

/* complex division
 * assumes small numbers different from zero
 * x = (a+bi)/(c+di); real(x) stored in *real; imag(d) in *imag */
void cdiv( double a, double b, double c, double d, double *real, double *imag ) {

    double cd;

    cd = 1 / (c*c + d*d) ;

    *real = ( a*c + b*d ) * cd;
    *imag = ( b*c - a*d ) * cd;
}

/* complex reciprocal
 * assumes small numbers different from zero
 * x = 1/(c+di); real(x) stored in *real; imag(d) in *imag */
void crep( double c, double d, double *real, double *imag ) {

    double cd;
    
    cd = 1 / (c*c + d*d) ;
    
    *real =  c * cd;
    *imag = -d * cd;
}

/* complex multiplication
 * assumes small numbers different from zero
 * (real+i*imag) = (a+bi)*(c+di) */
void cmul( double *real, double *imag, double a, double b, double c, double d ) {

    *real = a*c - b*d;
    *imag = a*d + b*c;
}
