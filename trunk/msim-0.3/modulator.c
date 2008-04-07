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

/* $Id: modulator.c,v 1.3 2004/08/16 17:26:49 miguel Exp $ */

/* Miguel Bazdresch

   MIMO Simulation

   Modulator  function

   This function modulates the frames produced by the coder.
   The function returns a 0 if everything went alright, and
   returns different than zero if there was an error.
   The pointer to the bits that are to be transmitted is
   pTransmitMatrix.

   The type of modulation is determined by variable
   cConstellationType, and the energy of the constellation is
   determined by fConstellationEnergy. */

int modulator (void) {

    /* function prototypes */
    void debug8(int a1, double d1, double d2);
    void debug32(void);

    /* local variables */
    unsigned int row, col, i;
    unsigned int mod_column;
    unsigned int current_symbol[5];

    /* case statement to handle the different modulations */
    switch(cConstellationType) {
        case 0: /* 16-QAM */
	    /* process four bits in pTransmitMatrix at a time.
	       Each transmitter has its own modulator, so modulation is done
	       one row at a time. */
            for (row=1;row<=M;row++) {
                col = 0;
                mod_column = 1;
                while (mod_column <= (Lt+L)) {
                    for(i=1;i<=4;i++){
                        current_symbol[i] = pTransmitMatrix[row][col+i];
                    }
                    if (current_symbol[1]==1) {
                        if (current_symbol[2]==1) {
                            if ((current_symbol[3] == 1 && current_symbol[4] == 1)) {
                                /* 15 */
                                pModulated[row][mod_column] = -e3;
                                pModulated[row+M][mod_column] = -e3;
                            }
                            if ((current_symbol[3] == 1 && current_symbol[4] == 0)) {
                                /* 14 */
                                pModulated[row][mod_column] = e3;
                                pModulated[row+M][mod_column] = -e3;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 1)) {
                                /* 13 */
                                pModulated[row][mod_column] = -e3;
                                pModulated[row+M][mod_column] = e3;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 0)) {
                                /* 12 */
                                pModulated[row][mod_column] = e3;
                                pModulated[row+M][mod_column] = e3;
                            }
                        }
                        else {
                            if ((current_symbol[3] == 1 && current_symbol[4] == 1)) {
                                /* 11 */
                                pModulated[row][mod_column] = -e1;
                                pModulated[row+M][mod_column] = -e3;
                            }
                            if ((current_symbol[3] == 1 && current_symbol[4] == 0)) {
                                /* 10 */
                                pModulated[row][mod_column] = e1;
                                pModulated[row+M][mod_column] = -e3;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 1)) {
                                /* 9 */
                                pModulated[row][mod_column] = -e1;
                                pModulated[row+M][mod_column] = e3;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 0)) {
                                /* 8 */
                                pModulated[row][mod_column] = e1;
                                pModulated[row+M][mod_column] = e3;
                            }
                        }
                    }
                    else {
                        if (current_symbol[2]==1) {
                            if ((current_symbol[3] == 1 && current_symbol[4] == 1)) {
                                /* 7 */
                                pModulated[row][mod_column] = -e3;
                                pModulated[row+M][mod_column] = -e1;
                            }
                            if ((current_symbol[3] == 1 && current_symbol[4] == 0)) {
                                /* 6 */
                                pModulated[row][mod_column] = e3;
                                pModulated[row+M][mod_column] = -e1;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 1)) {
                                /* 5 */
                                pModulated[row][mod_column] = -e3;
                                pModulated[row+M][mod_column] = e1;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 0)) {
                                /* 4 */
                                pModulated[row][mod_column] = e3;
                                pModulated[row+M][mod_column] = e1;
                            }
                        }
                        else {
                            if ((current_symbol[3] == 1 && current_symbol[4] == 1)) {
                                /* 3 */
                                pModulated[row][mod_column] = -e1;
                                pModulated[row+M][mod_column] = -e1;
                            }
                            if ((current_symbol[3] == 1 && current_symbol[4] == 0)) {
                                /* 2 */
                                pModulated[row][mod_column] = e1;
                                pModulated[row+M][mod_column] = -e1;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 1)) {
                                /* 1 */
                                pModulated[row][mod_column] = -e1;
                                pModulated[row+M][mod_column] = e1;
                            }
                            if ((current_symbol[3] == 0 && current_symbol[4] == 0)) {
                                /* 0 */
                                pModulated[row][mod_column] = e1;
                                pModulated[row+M][mod_column] = e1;
                            }
                        }
                    }
                    col += 4;
                    mod_column++;
                }
            }
            break;
	case 1: /* 4-QAM */
	    /* process two bits in pTransmitMatrix at a time.
	     * Each transmitter has its own modulator, so modulation is done
	     * one row at a time.
	     */
	    {
		for (row=1;row<=M;row++) {
		    col = 0;
		    mod_column = 1;
		    while (mod_column <= (Lt+L)) {
			for(i=1;i<=2;i++){
			    current_symbol[i] = pTransmitMatrix[row][col+i];
			}
			if( current_symbol[1] == 1 && current_symbol[2] == 1 ) {
			    pModulated[row][mod_column] = e1;
			    pModulated[row+M][mod_column] = e1;
			}
			if( current_symbol[1] == 1 && current_symbol[2] == 0 ) {
			    pModulated[row][mod_column] = e1;
			    pModulated[row+M][mod_column] = -e1;
			}
			if( current_symbol[1] == 0 && current_symbol[2] == 0 ) {
			    pModulated[row][mod_column] = -e1;
			    pModulated[row+M][mod_column] = -e1;
			}
			if( current_symbol[1] == 0 && current_symbol[2] == 1 ) {
			    pModulated[row][mod_column] = -e1;
			    pModulated[row+M][mod_column] = e1;
			}
			mod_column++;
			col+=2;
		    }
		}
	    }
	    break;
	default:
	    return 1;
    }
    
    if( debug ) {
	if ( debug && (cDebug & 32) == 32 )
	    debug32(); /* print modulated matrix if debug is enabled*/
	if ( debug && (cDebug & 8) == 8 )
	    debug8( 1, 0.0, 0.0 ); /* calculate power in modulated matrix */
    }
    
    return 0;
}
