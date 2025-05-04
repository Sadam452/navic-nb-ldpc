/*!
 * \file channel.c
 * \brief AWGN and Rayleigh channel
 */

#include <math.h>
#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
#include "./include/channel.h"



#define APSK
//#define QAM
//#define QAM_R //rotated QAM

//#define rayleigh_fading
//#define rayleigh_fading_SSD
//#define erasure

/*!
 * \fn void ModelChannel_AWGN_BPSK (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,int *init_rand)
 * \brief BPSK modulation + AWGN noise on the codeword.
 *                 The function computes the intrinsic_LLRs corresponding to the noisy observations.
 * Inputs
 * 	- NBIN : Binary copy of the codeword
 * 	- EbN  : Signal to noise ratio in terms of Eb/No (dB)
 * Outputs
 *      - decoder->intrinsic_LLR
 */
void ModelChannel_AWGN_BPSK (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,int *init_rand)
{
    const int N = code->N;
    const int logGF = code->logGF;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[4096];

    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) NoisyBin[n] = calloc(logGF,sizeof(float));

    /* Binary-input AWGN channel : */

    sigma = sqrt(1.0/(2*code->rate*pow(10,EbN/10.0)));  // considering EbNo
    //sigma = sqrt(1.0/(pow(10,EbN/10.0))); // considering SNR
    for (n=0; n<N; n++)
    {
        for (q=0; q<logGF; q++)
        {

            u=My_drand48(init_rand);
            while(u==0.0)
			{
                u=My_drand48(init_rand);
			}
            v=My_drand48(init_rand);
            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
            NoisyBin[n][q] = BPSK(NBIN[n][q]) + sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v);
        }

    }


    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {
        for (g=0; g<code->GF; g++)
        {
            TMP[g]=0.0;
            for (q=0; q<logGF; q++)
            {
                //TMP[g] = TMP[g] + SQR(NoisyBin[n][q]-BPSK(table->BINGF[g][q]))/(2.0*SQR(sigma));
                TMP[g] = TMP[g] - NoisyBin[n][q]*BPSK(table->BINGF[g][q]);
            }

        }

        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }


    }

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);
}