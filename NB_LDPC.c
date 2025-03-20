/*!
 * \file NB_LDPC.c
 * \brief Non-binary LDPC reduced complexity decoder with horizontal scheduling
 * \author Sadam Hussein
 * \copyright BSD copyright
 * \date 03/03/2015
 * \details
 
   This implements an Extended Min-Sum Decoder for Non-Binary LDPC codes.
   The decoder uses horizontal scheduling with a layered architecture
   and syndrome-based decoding approach.
   
   The algorithm provides an efficient reduced-complexity implementation
   for decoding non-binary LDPC codes in Galois Fields.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "./include/NB_LDPC.h"

/// preprocessing directives
//#define CCSK // use of Code-shift keying modulation

/*!
 * \fn int main(int argc, char * argv[])
 * \brief Main program for simulating non-binary LDPC codes
 *
 * This program simulates the performance of non-binary LDPC codes using
 * a reduced complexity decoder. It takes various parameters as input and
 * calculates performance metrics like Frame Error Rate (FER) and Bit Error Rate (BER).
 *
 * Inputs:
 *    NbMonteCarlo     : Number of frames to simulate
 *    NbIterMax        : Maximum number of decoding iterations
 *    FileMatrix       : Filename containing the parity-check matrix
 *    EbN              : Signal-to-noise ratio in dB (Eb/No)
 *    n_vc             : Size of truncated messages from Variable to Check nodes
 *    n_cv             : Size of truncated messages from Check to Variable nodes
 *    Offset           : Offset correction factor (typically 0.4 -- 1)
 *    NbOper           : Maximum number of operations for sorting
 *
 * Output:
 *    Frame Error Rate and other performance metrics for the given parameters
 *
 * Input File: 'FileMatrix' should be an ASCII file with parity-check matrix in aList format.
 * Output File: Results are saved in the ./data directory
 *
 * Usage example (Linux):
 *    compile using make
 *    ./essai 2000 10 ./matrices/KN/N576_K480_GF64.txt 3.5 10 20 0.3 25
 *
 * Example output:
 *    <0> FER = 40/751 = 0.053262 BER = 520 / x = 0.001443 avr_it = 2.58
 */
int main(int argc, char * argv[])
{
    int         k, l, n, iteration, i, g;
    int         **infoBits, *infoSymbols, **codeBits, *codeSymbols;
    int         *decodedSymbols, *codeWord;
    int         frameCount, maxIterations;
    float       snrValue;
    int         maxOperations, numMonteCarloTrials;
    float       offsetFactor;
    char        *outputFileName, *matrixFileName, *simulationName;
    int         syndrom = 0, numBitErrors, numFrameErrors = 0, numUndetectedErrors = 0;
    int         totalBitErrors = 0;

    code_t code;
    table_t table;
    decoder_t decoder;

    int checkNode;

    int randomSeed = -1; // initialization of random generator
    srand(2);

    /*
     * Command line arguments
     */
    //  printf("argc: %d\n", argc);
    if (argc < 9)
    {
        printf("File:\n %s\n ", argv[0]);
        //printf usage_txt of type: static char usage_txt[]
        printf("%s", usage_txt);

        return (EXIT_FAILURE);
    }
    outputFileName = malloc(STR_MAXSIZE);
    matrixFileName = malloc(STR_MAXSIZE);
    simulationName = malloc(STR_MAXSIZE);

    numMonteCarloTrials = atoi(argv[1]);
    maxIterations       = atoi(argv[2]);
    strcpy(matrixFileName, argv[3]);
    snrValue            = atof(argv[4]);
    decoder.n_vc        = atoi(argv[5]);
    decoder.n_cv        = atoi(argv[6]);
    offsetFactor        = atof(argv[7]);
    maxOperations       = atoi(argv[8]);
    
    printf(" Monte-Carlo simulation of Non-Binary LDPC decoder \n\n");
    printf("Simulation parameters:\n");
    printf("\n\t NbMonteCarlo     : %d", numMonteCarloTrials);
    printf("\n\t NbIterMax        : %d", maxIterations);
    printf("\n\t FileMatrix       : %s", matrixFileName);
    printf("\n\t snrValue       : %g", snrValue);
    printf("\n\t n_vc            : %d", decoder.n_vc);
    printf("\n\t n_cv            : %d", decoder.n_cv);
    printf("\n\t Offset           : %g", offsetFactor);
    printf("\n\t NbOper           : %d\n", maxOperations);

    //check if FileMatrix exists
    FILE *file;
    if ((file = fopen(matrixFileName, "r")) == NULL)
    {
        printf("File %s not found\n", matrixFileName);
        return (EXIT_FAILURE);
    }
    fclose(file);


    printf("Load code  ... ");
    LoadCode(matrixFileName, &code);
    printf(" OK \n Load table ...");
    LoadTables(&table, code.GF, code.logGF);
    printf("OK \n Allocate decoder ... ");
    AllocateDecoder(&code, &decoder);
    printf("OK \n Gaussian Elimination ... ");
    GaussianElimination(&code, &table);
    printf(" OK \n");

    // Output results to a file
    FILE *resultFile;
    char note[20] = "GF64_CCSK";
    printf("\n\t Note             : %s\n", note);
    char fileName[60];
    time_t startTime;
    time_t endTime;
    double executionTime;
    char* timeString;

    #ifdef CCSK
    // CCSK: build CCSK table
    csk_t csk;
    csk.PNsize = code.GF;
    printf("\n\t PN is generated using an LFSR \n");
    allocate_csk(&csk, csk.PNsize);
    PNGenerator(&csk); // Generate a PN sequence for csk modulation
    build_natural_csk_mapping(code.GF, &csk, table.BINGF); // Fills the csk_arr with shifted versions of PN sequence
    csk.PNsize = 6;  // For "punctured" CCSK mapping
    #endif

    sprintf(fileName, "./data/N%d_GF%d_nm%d_%s.txt", code.N, code.GF, decoder.n_cv, note);

    startTime = time(NULL);
    timeString = ctime(&startTime);
    printf("Simulation started at time: %s \n", timeString);

    /*
     * Memory allocation
     */
    codeBits = (int **)calloc(code.N, sizeof(int *));
    for (n = 0; n < code.N; n++)  
        codeBits[n] = (int *)calloc(code.logGF, sizeof(int));
    
    infoBits = (int **)calloc(code.K, sizeof(int *));
    for (k = 0; k < code.K; k++) 
        infoBits[k] = (int *)calloc(code.logGF, sizeof(int));

    codeSymbols = (int *)calloc(code.N, sizeof(int));
    infoSymbols = (int *)calloc(code.K, sizeof(int));
    codeWord = (int *)calloc(code.N, sizeof(int));
    decodedSymbols = (int *)calloc(code.N, sizeof(int));

    // Check that check node degree (dc) is constant
    int dcMin = 100;
    int dcMax = 0;
    for (checkNode = 0; checkNode < code.M; checkNode++)
    {
        if (dcMax < code.rowDegree[checkNode])
        {
            dcMax = code.rowDegree[checkNode];
        }
        if (dcMin > code.rowDegree[checkNode])
        {
            dcMin = code.rowDegree[checkNode];
        }
    }
    if (dcMin != dcMax)
    {
        printf("d_c is not constant: dc_min= %d ; dc_max=%d !!!!!! \n", dcMin, dcMax);
    }

    int totalIterations;
    softdata_t messageVtoC_temp[dcMax][code.GF];
    softdata_t messageVtoC_temp2[dcMax][code.GF];
    totalIterations = 0;

    for (frameCount = 1; frameCount <= numMonteCarloTrials; frameCount++)
    {
        /* Generate uniformly distributed information bits */
        RandomBinaryGenerator(code.N, code.M, code.GF, code.logGF, infoBits, infoSymbols, table.BINGF, &randomSeed);

        /* Encode the information bits to a (non binary) codeword */
        Encoding(&code, &table, codeWord, codeBits, infoSymbols);

        /* Noisy channel (AWGN) */
        #ifdef CCSK
        ModelChannel_AWGN_BPSK_CSK(&csk, &code, &decoder, &table, codeWord, snrValue, &randomSeed);
        #endif
        #ifndef CCSK
        ModelChannel_AWGN_BPSK(&code, &decoder, &table, codeBits, snrValue, &randomSeed);
        #endif

        /***********************************************/
        /* Implementation of the horizontal scheduling */

        // Initialize VtoC with intrinsic values
        for (n = 0; n < code.N; n++)
        {
            for (k = 0; k < code.GF; k++)
            {
                decoder.VtoC[n][k] = decoder.intrinsic_LLR[n][k];
            }
        }

        /* Decoding iterations */
        for (iteration = 0; iteration < maxIterations - 1; iteration++)
        {
            for (checkNode = 0; checkNode < code.M; checkNode++) /* Loop for the M Check nodes */
            {
                for (i = 0; i < code.rowDegree[checkNode]; i++)
                {
                    for (k = 0; k < code.GF; k++)
                    {
                        messageVtoC_temp[i][k] = decoder.VtoC[code.mat[checkNode][i]][k];
                        messageVtoC_temp2[i][k] = messageVtoC_temp[i][k];
                    }
                }

                // Sorting Mvc values
                for (i = 0; i < code.rowDegree[checkNode]; i++)
                {
                    for (k = 0; k < decoder.n_vc; k++)
                    {
                        decoder.M_VtoC_LLR[i][k] = +1e5;
                        decoder.M_VtoC_GF[i][k] = 0;
                        for (g = 0; g < code.GF; g++)
                        {
                            if (messageVtoC_temp2[i][g] < decoder.M_VtoC_LLR[i][k])
                            {
                                decoder.M_VtoC_LLR[i][k] = messageVtoC_temp2[i][g];
                                decoder.M_VtoC_GF[i][k] = g;
                            }
                        }
                        messageVtoC_temp2[i][decoder.M_VtoC_GF[i][k]] = +1e5;
                    }

                    // Normalization
                    for (g = 1; g < decoder.n_vc; g++)
                    {
                        decoder.M_VtoC_LLR[i][g] = decoder.M_VtoC_LLR[i][g] - decoder.M_VtoC_LLR[i][0];
                    }
                    decoder.M_VtoC_LLR[i][0] = 0.0;
                }

                CheckPassLogEMS(checkNode, &decoder, &code, &table, maxOperations, offsetFactor);

                // Compute soft output
                for (i = 0; i < code.rowDegree[checkNode]; i++)
                {
                    for (k = 0; k < code.GF; k++)
                    {
                        decoder.APP[code.mat[checkNode][i]][k] = decoder.M_CtoV_LLR[i][k] + messageVtoC_temp[i][k];
                        decoder.VtoC[code.mat[checkNode][i]][k] = decoder.M_CtoV_LLR[i][k] + decoder.intrinsic_LLR[code.mat[checkNode][i]][k]; // Compute Mvc and save RAM
                    }
                }
            } /* End of the node update */

            Decision(decodedSymbols, decoder.APP, code.N, code.GF);
            syndrom = Syndrom(&code, decodedSymbols, &table);
            
            if (syndrom == 0)
                break;
        }

        totalIterations = totalIterations + iteration + 1;

        /* Compute the Bit Error Rate (BER) */
        numBitErrors = 0;
        for (k = 0; k < code.K; k++)
        {
            for (l = 0; l < code.logGF; l++)
                if (table.BINGF[decodedSymbols[k]][l] != codeBits[k][l])
                    numBitErrors++;
        }

        totalBitErrors = totalBitErrors + numBitErrors;
        if (numBitErrors != 0)
        {
            numFrameErrors++;
            if (syndrom == 0)
                numUndetectedErrors++;
        }
        
        if (frameCount % 10 == 0)
        {
            printf("\r<%d> FER= %d / %d = %f BER= %d / x = %f  avr_it=%.2f",
                numUndetectedErrors, numFrameErrors, frameCount, (double)(numFrameErrors) / frameCount, 
                totalBitErrors, (double)totalBitErrors / (frameCount * code.K * code.logGF), 
                (double)(totalIterations) / frameCount);
            fflush(stdout);
        }

        if (numFrameErrors == 40)
            break;
    }

    printf("\r<%d> FER= %d / %d = %f BER= %d / x = %f  avr_it=%.2f",
        numUndetectedErrors, numFrameErrors, frameCount, (double)(numFrameErrors) / frameCount,
        totalBitErrors, (double)totalBitErrors / (frameCount * code.K * code.logGF),
        (double)(totalIterations) / frameCount);

    printf(" \n results are printed in file %s \n", fileName);

    endTime = time(NULL);
    timeString = ctime(&endTime);
    executionTime = difftime(endTime, startTime);
    resultFile = fopen(fileName, "a");
    
    if ((resultFile) == NULL)
    {
        printf(" \n !! file not found \n ");
    }
    else
    {
        fprintf(resultFile, " SNR:%.2f: \t FER= %d / %d = %f ", snrValue, numFrameErrors, frameCount, (double)(numFrameErrors) / frameCount);
        fprintf(resultFile, " \t BER= %d / x = \t %f  avr_it= \t %.2f \t time: %s", totalBitErrors, 
                (double)totalBitErrors / (double)(frameCount * code.K * code.logGF), 
                (double)(totalIterations) / frameCount, timeString);
    }
    fclose(resultFile);

    printf("\n");
    printf("Simulation complete at time: %s", timeString);
    printf("execution time:%0.2f", executionTime);

    getchar();

    // Free allocated memory
    free(outputFileName);
    free(matrixFileName);
    free(simulationName);
    free(decodedSymbols);
    free(codeWord);
    free(infoSymbols);
    free(codeSymbols);

    for (n = 0; n < code.N; n++) free(codeBits[n]);
    free(codeBits);
    for (k = 0; k < code.K; k++) free(infoBits[k]);
    free(infoBits);

    FreeCode(&code);
    FreeDecoder(&decoder);
    FreeTable(&table);
    return (EXIT_SUCCESS);
}