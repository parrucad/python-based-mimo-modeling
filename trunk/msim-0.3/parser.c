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

/* $Id: parser.c,v 1.6 2004/08/12 20:42:27 miguel Exp $ */

/* MIMO Simulation */

/* Miguel Bazdresch */

/* Parser  function */
/* This function parses a file named mimosim.ini, and */
/* changes the value of the control variables accordingly. */
/* The function returns a 0 if everything went alright, and */
/* returns different than zero if there was an error. */

int parser (void) {

    /* variable declaration */
    FILE *comfile;
    char pCommand[30];
    /*    char cArg; */
    int iArg;
    long lArg;
    double fArg;

    /* open the file */
    if ((comfile = fopen("mimosim.ini", "r")) == NULL) {
        printf("  Parser: Error opening file mimosim.ini!\n");
        fclose(comfile);
        return 1;
    }
    else {
        while (feof(comfile)==0) {
            strcpy(pCommand, "\0");
            if (fscanf(comfile, "%s", pCommand)) {
                if(strcmp(pCommand, "\0"))
                    printf("  Parser: Command read is: %s\n", pCommand);
                if (!strcmp(pCommand, "Debug")) {
                    fscanf(comfile, "%d", &iArg);
                    cDebug = iArg;
                    printf("  Parser: set Debug variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "ReceiverType")) {
                    fscanf(comfile, "%d", &iArg);
                    iReceiverType = iArg;
                    printf("  Parser: set ReceiverType variable to: %u\n", iArg);
                }
                if ((!strcmp(pCommand, "CycleVblastTypes")) ||
		       	(!strcmp(pCommand, "cyclevblasttypes"))) {
                    fscanf(comfile, "%d", &iArg);
                    CycleVblastTypes = iArg;
                    printf("  Parser: set CycleVblastTypes variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "ConstellationSymbolEnergy")) {
                    fscanf(comfile, "%lf", &fArg);
                    fConstellationSymbolEnergy = fArg;
                    printf("  Parser: set ConstellationSymbolEnergy variable to: %f\n", fArg);
                }
                if (!strcmp(pCommand, "Lt")) {
                    fscanf(comfile, "%d", &iArg);
                    Lt = iArg;
                    printf("  Parser: set Lt variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand,  "L")) {
                    fscanf(comfile, "%d", &iArg);
                    L = iArg;
                    printf("  Parser: set L variable to: %u\n", iArg);
                }
		if ((!strcmp(pCommand, "CycleL")) || (!strcmp(pCommand, "cyclel"))) {
		    fscanf(comfile, "%d", &iArg);
		    CycleL = iArg;
		    printf("  Parser: set CycleL variable to: %u\n", iArg);
		}
                if (!strcmp(pCommand, "ConstellationType")) {
                    fscanf(comfile, "%d", &iArg);
                    cConstellationType = iArg;
                    printf("  Parser: set ConstellationType variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "SourceType")) {
                    fscanf(comfile, "%d", &iArg);
                    cSourceType = iArg;
                    printf("  Parser: set SourceType variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "SourceSeed")) {
                    fscanf(comfile, "%ld", &lArg);
                    iSourceSeed = lArg;
                    printf("  Parser: set SourceSeed variable to: %lu\n", lArg);
                }
                if (!strcmp(pCommand, "NoiseSeed")) {
                    fscanf(comfile, "%ld", &lArg);
                    iNoiseSeed = lArg;
                    printf("  Parser: set NoiseSeed variable to: %lu\n", lArg);
                }
                if (!strcmp(pCommand, "TransmitAntennas")) {
                    fscanf(comfile, "%d", &iArg);
                    cTransmitAntennas = iArg;
                    printf("  Parser: set TransmitAntennas variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "ReceiveAntennas")) {
                    fscanf(comfile, "%d", &iArg);
                    cReceiveAntennas = iArg;
                    printf("  Parser: set ReceiveAntennas variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "NoisePower")) {
                    fscanf(comfile, "%lf", &fArg);
                    fNoisePower = fArg;
                    printf("  Parser: set NoisePower variable to: %f\n", fArg);
                }
                if (!strcmp(pCommand, "MatrixHComponentGen")) {
                    fscanf(comfile, "%d", &iArg);
                    cMatrixHComponentGen = iArg;
                    printf("  Parser: set MatrixHComponentGen variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "HSeed")) {
                    fscanf(comfile, "%ld", &lArg);
                    iHSeed = lArg;
                    printf("  Parser: set HSeed variable to: %lu\n", lArg);
                }
                if (!strcmp(pCommand, "MethodChannelEstimation")) {
                    fscanf(comfile, "%d", &iArg);
                    cMethodChannelEstimation = iArg;
                    printf("  Parser: set MethodChannelEstimation variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "ChannelCode")) {
                    fscanf(comfile, "%d", &iArg);
                    cChannelCode = iArg;
                    printf("  Parser: set ChannelCode variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "TypeTrainingSeq")) {
                    fscanf(comfile, "%d", &iArg);
                    cTypeTrainingSeq = iArg;
                    printf("  Parser: set TypeTrainingSeq variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "NumberOfErrors")) {
                    fscanf(comfile, "%d", &iArg);
                    iNumberOfErrors = iArg;
                    printf("  Parser: set NumberOfErrors variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "RunFrames")) {
                    fscanf(comfile, "%d", &iArg);
                    iRunFrames = iArg;
                    printf("  Parser: set RunFrames variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "Resume")) {
                    fscanf(comfile, "%d", &iArg);
                    iResume = iArg;
                    printf("  Parser: set Resume variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "experimental")) {
                    fscanf(comfile, "%d", &iArg);
                    experimental = iArg;
                    printf("  Parser: set experimental variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "exp_vb_ordering")) {
                    fscanf(comfile, "%d", &iArg);
                    exp_vb_ordering = iArg;
                    printf("  Parser: set exp_vb_ordering variable to: %u\n", iArg);
                }
                if (!strcmp(pCommand, "exp_vb_cancellation")) {
                    fscanf(comfile, "%d", &iArg);
		    exp_vb_cancellation = iArg;
		    printf("  Parser: set exp_vb_cancellation variable to: %u\n", iArg);
		}
		if ((!strcmp(pCommand, "NIterations")) ||
		       	(!strcmp(pCommand, "niterations"))) {
		    fscanf(comfile, "%d", &iArg);
		    NIterations = iArg;
		    printf("  Parser: set NIterations variable to: %u\n", iArg);
		}
		if ((!strcmp(pCommand, "NStep")) || (!strcmp(pCommand, "nstep"))) {
		    fscanf(comfile, "%lf", &fArg);
		    NStep = fArg;
		    printf("  Parser: set NStep variable to: %u\n", iArg);
		}
		if (!strcmp(pCommand, "DoLLL")) {
		    fscanf(comfile, "%d", &iArg);
		    DoLLL = iArg;
		    printf("  Parser: set DoLLL variable to: %u\n", iArg);
		}
		if (!strcmp(pCommand, "CycleLLL")) {
		    fscanf(comfile, "%d", &iArg);
		    CycleLLL = iArg;
		    printf("  Parser: set CycleLLL variable to: %u\n", iArg);
		}
		if (!strcmp(pCommand, "ScreenFeedback")) {
		    fscanf(comfile, "%d", &iArg);
		    ScreenFeedback = iArg;
		    printf("  Parser: set ScreenFeedback variable to: %u\n", iArg);
		}
		if ( (!strcmp(pCommand, "LIterations")) ||
		     (!strcmp(pCommand, "literations")) ) {
		    fscanf(comfile, "%d", &iArg);
		    LIterations = iArg;
		    printf("  Parser: set LIterations variable to: %u\n", iArg);
		}
		if ( (!strcmp(pCommand, "LStep")) ||
		     (!strcmp(pCommand, "lstep")) ) {
		    fscanf(comfile, "%d", &iArg);
		    LStep = iArg;
		    printf("  Parser: set LStep variable to: %u\n", iArg);
		}
		if (!strcmp(pCommand, "IndexComFile")) {
		    fscanf(comfile, "%d", &iArg);
		    IndexComFile = iArg;
		    printf("  Parser: set IndexComFile variable to: %u\n", iArg);
		}
	    }
	    else {
		printf("Error reading from file.\n");
		fclose(comfile);
		return 1;
	    }
	}
	return 0;
    }
}
