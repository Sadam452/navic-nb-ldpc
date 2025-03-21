# C code for NB-LDPC simulation

You can use this code to simulate NB-LDPC matrices using the Extented-Min Sum (EMS) algorithm.
 The Check Node (CN) is processed using Forward Backward(FB) algorithm. The FB algorithm splits CN in elementary CNs (ECN).

# usage

## input argument
 there are 8 arguments
 1.		NbMonteCarlo     : # simulated frames
 1.		NbIterMax        : # of maximum decoding iteration
 1.		FileMatrix       : File name of the parity-check matrix
 1.		snrValue         : snrValue
 1.		n_vc             : size of truncated messages from Variable to Check
 1.		n_cv	         : size of truncated messages from Check to Variable
 1.		Offset           : offset correction factor (0.4 -- 1)
 1.		NbOper           : Maximum number of operations for sorting
 
 ## output

Frame Error Rate (FER) and Bit Error Rate (BER) and average number of iterations for the given simulation parameters

## input and output files
 * Input File : 'FileMatrix' is an ASCII file with the parity-check matrix description in aList format.
 * Output File : a txt file giving in the ./data forder giving SNR, BER, average number of iterations and time of the end of simulation

## Simulation on Linux

> make
> ./essai 2000 10 ./matrices/KN/N576_K480_GF64.txt 3.5 20 20 0.3 25

## Simulation results

> Monte-Carlo simulation of Non-Binary LDPC decoder
>
> Simulation parameters:
>
>        NbMonteCarlo     : 2000
>         NbIterMax        : 10
 >        FileMatrix       : ./matrices/KN/N1200_K600_GF64_Navic.txt
>         Eb/No (dB)       : 3.5
>         n_vc            : 20
>         n_cv            : 20
>         Offset           : 0.3
>         NbOper           : 25
> Load code  ...  
> Normal alist format is used! 
> LDPC code parameters: 
>         N      :200 
>         K      :100 
>         M      :100 
>         CR     :0.5 
>         GF     :64 
>         logGF  :6
> OK 
> Load table ...OK 
> Allocate decoder ... OK 
> Gaussian Elimination ...  OK 
>
>         Note             : GF64_CCSK
> Simulation started at time: Fri Mar 21 01:53:13 2025
 
> <0> FER= 0 / 2001 = 0.000000 BER= 0 / x = 0.000000  avr_it=2.02 
> results are printed in file ./data/N200_GF64_nm20_GF64_CCSK.txt 

> Simulation complete at time: Fri Mar 21 01:53:21 2025
> execution time:8.00