/*!
 * \file bubble_decoder.c
 * \brief Implementation of check node processing using the bubble decoder algorithm
 * 
 * This file contains implementations for the bubble check node processing
 * algorithm, which is used in decoding non-binary LDPC codes.
 */

#include <stdlib.h>
#include <stdio.h>
#include "./include/syndrome_decoder.h"
#include "./include/bubble_decoder.h"

/**
 * \fn find_maximum
 * \brief Finds the index of the maximum value in an array
 * 
 * @param values Array of float values to search
 * @param length Length of the array
 * @return Index of the maximum value
 */
int find_maximum(float *values, int length)
{
    float max_value = values[0];
    int max_index = 0;
    
    for (int i = 1; i < length; i++) {
        if (values[i] > max_value) {
            max_value = values[i];
            max_index = i;
        }
    }
    
    return max_index;
}

/**
 * \fn find_minimum
 * \brief Finds the index of the minimum value in an array
 * 
 * @param values Array of float values to search
 * @param length Length of the array
 * @return Index of the minimum value
 */
int find_minimum(float *values, int length)
{
    float min_value = values[0];
    int min_index = 0;
    
    for (int i = 1; i < length; i++) {
        if (values[i] < min_value) {
            min_value = values[i];
            min_index = i;
        }
    }
    
    return min_index;
}

/**
 * \fn ElementaryStep
 * \brief Performs elementary bubble check node processing
 * 
 * This function implements the core bubble check node algorithm, combining
 * two input messages to produce an output message.
 * 
 * @param input1 First input LLR values
 * @param input2 Second input LLR values
 * @param indice_input1 First input GF indices
 * @param indice_input2 Second input GF indices
 * @param output Output LLR values
 * @param indice_out Output GF indices
 * @param add_gf GF addition table
 * @param gf_size Size of the Galois Field
 * @param max_size Maximum size of inputs/outputs
 * @param num_operations Number of operations to perform
 * @return 0 on success
 */
int ElementaryStep(float *input1, float *input2, int *indice_input1, int *indice_input2,
                   float *output, int *indice_out, int **add_gf, int gf_size, 
                   int max_size, int num_operations)
{
    const int NUM_BUBBLES = 8;
    const int HALF_BUBBLES = NUM_BUBBLES / 2;
    const float INFINITY_VALUE = 1e5;
    
    // Initialize temporary arrays for output
    float temp_output[max_size];
    int temp_indice_out[max_size];
    
    // Track which GF values have been added to output
    int *gf_values_seen = (int *)calloc(gf_size, sizeof(int));
    for (int i = 0; i < gf_size; i++) {
        gf_values_seen[i] = -1;
    }
    
    // Initialize output arrays
    for (int i = 0; i < max_size; i++) {
        temp_output[i] = INFINITY_VALUE;
        temp_indice_out[i] = -1;
    }
    
    // Pre-compute additions for initial bubbles
    float addition_table[max_size][max_size];
    
    // Horizontal pre-computation
    for (int j = 0; j < HALF_BUBBLES; j++) {
        for (int i = 0; i < max_size; i++) {
            addition_table[i][j] = input1[i] + input2[j]; // Min-sum algorithm
        }
    }
    
    // Vertical pre-computation
    for (int i = 0; i < HALF_BUBBLES; i++) {
        for (int j = HALF_BUBBLES; j < max_size; j++) {
            addition_table[i][j] = input1[i] + input2[j]; // Min-sum algorithm
        }
    }
    
    // Initialize competitor matrix
    float competitor_values[NUM_BUBBLES];
    int competitor_row[NUM_BUBBLES];
    int competitor_col[NUM_BUBBLES];
    
    // Vertical competitors
    for (int j = 0; j < HALF_BUBBLES; j++) {
        competitor_values[j] = addition_table[j][0];
        competitor_row[j] = j;
        competitor_col[j] = 0;
    }
    
    // Horizontal competitors
    for (int j = 0; j < HALF_BUBBLES; j++) {
        competitor_values[j + HALF_BUBBLES] = addition_table[HALF_BUBBLES][j];
        competitor_row[j + HALF_BUBBLES] = HALF_BUBBLES;
        competitor_col[j + HALF_BUBBLES] = j;
    }
    
    // Fill output arrays
    int output_index = 0;
    
    for (int op = 0; op < num_operations; op++) {
        int pos = find_minimum(competitor_values, NUM_BUBBLES);
        
        // Check for out of bounds
        if ((indice_input1[competitor_row[pos]] == -1) || (indice_input2[competitor_col[pos]] == -1)) {
            break;
        }
        
        // Calculate GF addition (XOR for binary case)
        int gf_result = indice_input1[competitor_row[pos]] ^ indice_input2[competitor_col[pos]];
        
        // Add to output if not already present
        if (gf_values_seen[gf_result] == -1) {
            temp_output[output_index] = competitor_values[pos];
            temp_indice_out[output_index] = gf_result;
            gf_values_seen[gf_result] = 1;
            output_index++;
        }
        
        // Stop if we've filled the output
        if (output_index == max_size) {
            break;
        }
        
        // Update competitor with next value
        if (pos >= HALF_BUBBLES) {
            competitor_row[pos]++;
        } else {
            competitor_col[pos]++;
        }
        
        // Check bounds
        if (competitor_row[pos] >= max_size || competitor_col[pos] >= max_size) {
            competitor_values[pos] = INFINITY_VALUE;
        } else {
            competitor_values[pos] = addition_table[competitor_row[pos]][competitor_col[pos]];
        }
    }
    
    // Copy temporary arrays to output
    for (int i = 0; i < max_size; i++) {
        output[i] = temp_output[i];
        indice_out[i] = temp_indice_out[i];
    }
    
    free(gf_values_seen);
    return 0;
}

/**
 * \fn ElementaryStep_nm
 * \brief Performs elementary bubble check node processing with different input/output sizes
 * 
 * This variant of ElementaryStep handles different sizes for input and output arrays.
 * 
 * @param input1 First input LLR values
 * @param input2 Second input LLR values
 * @param indice_input1 First input GF indices
 * @param indice_input2 Second input GF indices
 * @param output Output LLR values
 * @param indice_out Output GF indices
 * @param add_gf GF addition table
 * @param gf_size Size of the Galois Field
 * @param input1_size Size of the first input
 * @param input2_size Size of the second input
 * @param output_size Size of the output
 * @param num_operations Number of operations to perform
 * @return 0 on success
 */
int ElementaryStep_nm(float *input1, float *input2, int *indice_input1, int *indice_input2,
                      float *output, int *indice_out, int **add_gf, int gf_size,
                      int input1_size, int input2_size, int output_size, int num_operations)
{
    const int NUM_BUBBLES = 8;
    const int HALF_BUBBLES = NUM_BUBBLES / 2;
    const float INFINITY_VALUE = 1e5;
    
    // Initialize temporary arrays for output
    float temp_output[output_size];
    int temp_indice_out[output_size];
    
    // Track which GF values have been added to output
    int *gf_values_seen = (int *)calloc(gf_size, sizeof(int));
    for (int i = 0; i < gf_size; i++) {
        gf_values_seen[i] = -1;
    }
    
    // Initialize output arrays
    for (int i = 0; i < output_size; i++) {
        temp_output[i] = INFINITY_VALUE;
        temp_indice_out[i] = -1;
    }
    
    // Pre-compute additions for initial bubbles
    float addition_table[input1_size + 1][input2_size + 1];
    
    // Horizontal pre-computation
    for (int j = 0; j < HALF_BUBBLES; j++) {
        for (int i = 0; i < input1_size; i++) {
            addition_table[i][j] = input1[i] + input2[j]; // Min-sum algorithm
        }
        addition_table[input1_size][j] = INFINITY_VALUE;
    }
    
    // Vertical pre-computation
    for (int i = 0; i < HALF_BUBBLES; i++) {
        for (int j = HALF_BUBBLES; j < input2_size; j++) {
            addition_table[i][j] = input1[i] + input2[j]; // Min-sum algorithm
        }
        addition_table[i][input2_size] = INFINITY_VALUE;
    }
    
    // Initialize competitor matrix
    float competitor_values[NUM_BUBBLES];
    int competitor_row[NUM_BUBBLES];
    int competitor_col[NUM_BUBBLES];
    
    // Vertical competitors
    for (int j = 0; j < HALF_BUBBLES; j++) {
        competitor_values[j] = addition_table[j][0];
        competitor_row[j] = j;
        competitor_col[j] = 0;
    }
    
    // Horizontal competitors
    for (int j = 0; j < HALF_BUBBLES; j++) {
        competitor_values[j + HALF_BUBBLES] = addition_table[HALF_BUBBLES][j];
        competitor_row[j + HALF_BUBBLES] = HALF_BUBBLES;
        competitor_col[j + HALF_BUBBLES] = j;
    }
    
    // Fill output arrays
    int output_index = 0;
    
    for (int op = 0; op < num_operations; op++) {
        int pos = find_minimum(competitor_values, NUM_BUBBLES);
        
        // Check for out of bounds
        if ((indice_input1[competitor_row[pos]] == -1) || (indice_input2[competitor_col[pos]] == -1)) {
            break;
        }
        
        // Calculate GF addition (XOR for binary case)
        int gf_result = indice_input1[competitor_row[pos]] ^ indice_input2[competitor_col[pos]];
        
        // Add to output if not already present
        if (gf_values_seen[gf_result] == -1) {
            temp_output[output_index] = competitor_values[pos];
            temp_indice_out[output_index] = gf_result;
            gf_values_seen[gf_result] = 1;
            output_index++;
        }
        
        // Stop if we've filled the output
        if (output_index == output_size) {
            break;
        }
        
        // Update competitor with next value
        if (pos >= HALF_BUBBLES) {
            competitor_row[pos]++;
        } else {
            competitor_col[pos]++;
        }
        
        // Check bounds
        if (competitor_row[pos] >= input1_size || competitor_col[pos] >= input2_size) {
            competitor_values[pos] = INFINITY_VALUE;
        } else {
            competitor_values[pos] = addition_table[competitor_row[pos]][competitor_col[pos]];
        }
    }
    
    // Copy temporary arrays to output
    for (int i = 0; i < output_size; i++) {
        output[i] = temp_output[i];
        indice_out[i] = temp_indice_out[i];
    }
    
    free(gf_values_seen);
    return 0;
}

/**
 * \fn CheckPassLogEMS
 * \brief Process the Check to Variable messages in the decoding graph
 * 
 * This function implements the check node processing for general row degree.
 * It processes Variable to Check messages and produces Check to Variable messages
 * using the bubble decoder algorithm.
 * 
 * @param node Current check node index
 * @param decoder Decoder state
 * @param code Code parameters
 * @param table Lookup tables for GF operations
 * @param num_operations Number of operations to perform
 * @param offset Offset value for LLR calculation
 */
void CheckPassLogEMS(int node, decoder_t *decoder, code_t *code, table_t *table, 
                     int num_operations, float offset)
{
    const int max_size = decoder->n_vc;
    const float INFINITY_VALUE = 1e5;
    int row_degree = code->rowDegree[node];
    int num_steps = 2 * (row_degree - 2);
    int t, k, g, k1, Stp=0;
    
    // Allocate temporary buffers for forward/backward recursion
    float *forward = (float *)calloc(max_size, sizeof(float));
    int *forward_indices = (int *)calloc(max_size, sizeof(int));
    float *backward = (float *)calloc(max_size, sizeof(float));
    int *backward_indices = (int *)calloc(max_size, sizeof(int));
    
    float *forward1 = (float *)calloc(max_size, sizeof(float));
    int *forward1_indices = (int *)calloc(max_size, sizeof(int));
    float *backward1 = (float *)calloc(max_size, sizeof(float));
    int *backward1_indices = (int *)calloc(max_size, sizeof(int));
    
    // Allocate intermediate results storage
    float **intermediate_values = (float **)calloc(num_steps, sizeof(float *));
    int **intermediate_indices = (int **)calloc(num_steps, sizeof(int *));
    
    for (k = 0; k < num_steps; k++) {
        intermediate_values[k] = (float *)calloc(max_size, sizeof(float));
        intermediate_indices[k] = (int *)calloc(max_size, sizeof(int));
        
        for (k1 = 0; k1 < max_size; k1++) {
            intermediate_values[k][k1] = INFINITY_VALUE;
            intermediate_indices[k][k1] = -1;
        }
    }
    
    // Initialize temporary buffers
    for (k = 0; k < max_size; k++) {
        forward[k] = INFINITY_VALUE;
        backward[k] = INFINITY_VALUE;
        forward_indices[k] = -1;
        backward_indices[k] = -1;
    }
    
    // Apply GF rotation to variable-to-check messages
    for (t = 0; t < row_degree; t++) {
        for (k = 0; k < max_size; k++) {
            if (decoder->M_VtoC_GF[t][k] != -1) {
                int gf_value = decoder->M_VtoC_GF[t][k];
                // GF multiplication
                gf_value = table->MULDEC[gf_value][code->matValue[node][t]];
                decoder->M_VtoC_GF[t][k] = gf_value;
            }
        }
    }
    
    // Initialize forward/backward messages
    for (k = 0; k < max_size; k++) {
        forward[k] = decoder->M_VtoC_LLR[0][k];
        forward_indices[k] = decoder->M_VtoC_GF[0][k];
        backward[k] = decoder->M_VtoC_LLR[row_degree - 1][k];
        backward_indices[k] = decoder->M_VtoC_GF[row_degree - 1][k];
    }
    
    // Forward-backward recursion
    for (int kk = 1; kk < row_degree - 1; kk++) {
        // Copy current messages
        for (k = 0; k < max_size; k++) {
            forward1[k] = decoder->M_VtoC_LLR[kk][k];
            forward1_indices[k] = decoder->M_VtoC_GF[kk][k];
            backward1[k] = decoder->M_VtoC_LLR[row_degree - kk - 1][k];
            backward1_indices[k] = decoder->M_VtoC_GF[row_degree - kk - 1][k];
        }
        
        // Store intermediate values
        for (k = 0; k < max_size; k++) {
            intermediate_values[kk - 1][k] = forward[k];
            intermediate_indices[kk - 1][k] = forward_indices[k];
            intermediate_values[num_steps - kk][k] = backward[k];
            intermediate_indices[num_steps - kk][k] = backward_indices[k];
        }
        
        // Forward step
        ElementaryStep(forward, forward1, forward_indices, forward1_indices,
                       forward, forward_indices, table->ADDGF, code->GF, 
                       max_size, num_operations);
        
        // Backward step
        ElementaryStep(backward, backward1, backward_indices, backward1_indices,
                       backward, backward_indices, table->ADDGF, code->GF, 
                       max_size, num_operations);
    }
    
    // Update the first and last CtoV messages
    for (k = 0; k < max_size; k++) {
        // Last vector
        decoder->M_CtoV_LLR[row_degree - 1][k] = forward[k];
        decoder->M_CtoV_GF[row_degree - 1][k] = forward_indices[k];
        // First vector
        decoder->M_CtoV_LLR[0][k] = backward[k];
        decoder->M_CtoV_GF[0][k] = backward_indices[k];
    }
    
    // Update intermediate CtoV messages using merging
    for (k = 0; k < row_degree - 2; k++) {
        ElementaryStep(intermediate_values[k], intermediate_values[row_degree - 2 + k],
                      intermediate_indices[k], intermediate_indices[row_degree - 2 + k],
                      forward, forward_indices, table->ADDGF, code->GF, 
                      max_size, num_operations);
        
        for (g = 0; g < max_size; g++) {
            decoder->M_CtoV_LLR[k + 1][g] = forward[g];
            decoder->M_CtoV_GF[k + 1][g] = forward_indices[g];
        }
    }
    
    // Apply inverse GF rotation to check-to-variable messages
    float llr_tmp[code->GF];
    
    for (t = 0; t < row_degree; t++) {
        // Find the stopping point
        for (k = 0; k < max_size; k++) {
            if (decoder->M_CtoV_GF[t][k] == -1) {
                Stp = k;
                break;
            } else {
                Stp = max_size;
            }
        }
        
        // GF division
        for (k = 0; k < Stp; k++) {
            decoder->M_CtoV_GF[t][k] = table->DIVDEC[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }
        
        if (decoder->M_CtoV_LLR[t][Stp - 1] > 0) {
            // Apply offset to all LLR values
            for (k = 0; k < code->GF; k++) {
                llr_tmp[k] = decoder->M_CtoV_LLR[t][Stp - 1] + offset;
            }
            
            // Update specific LLR values
            for (k = 0; k < Stp; k++) {
                llr_tmp[decoder->M_CtoV_GF[t][k]] = decoder->M_CtoV_LLR[t][k];
            }
            
            // Update all CtoV messages
            for (k = 0; k < code->GF; k++) {
                decoder->M_CtoV_LLR[t][k] = llr_tmp[k];
            }
        } else {
            // Reset all LLR values
            for (k = 0; k < code->GF; k++) {
                decoder->M_CtoV_LLR[t][k] = 0;
            }
        }
        
        // Update GF indices
        for (k = 0; k < code->GF; k++) {
            decoder->M_CtoV_GF[t][k] = k;
        }
    }
    
    // Free allocated memory
    for (k = 0; k < num_steps; k++) {
        free(intermediate_indices[k]);
        free(intermediate_values[k]);
    }
    free(intermediate_indices);
    free(intermediate_values);
    free(forward);
    free(forward_indices);
    free(backward);
    free(backward_indices);
    free(forward1);
    free(forward1_indices);
    free(backward1);
    free(backward1_indices);
}

/**
 * \fn CheckPassLogEMS_dc3
 * \brief Process the Check to Variable messages for row degree 3
 * 
 * This function is an optimized version of CheckPassLogEMS specifically for
 * check nodes with row degree 3.
 * 
 * @param node Current check node index
 * @param decoder Decoder state
 * @param code Code parameters
 * @param table Lookup tables for GF operations
 * @param num_operations Number of operations to perform
 * @param offset Offset value for LLR calculation
 */
void CheckPassLogEMS_dc3(int node, decoder_t *decoder, code_t *code, table_t *table, 
                         int num_operations, float offset)
{
    const int n_vc = decoder->n_vc;
    const int n_cv = decoder->n_cv;
    int t, c, k, Stp=0;
    float llr_tmp[code->GF];
    
    // Apply GF rotation to variable-to-check messages
    for (t = 0; t < code->rowDegree[node]; t++) {
        for (k = 0; k < n_vc; k++) {
            if (decoder->M_VtoC_GF[t][k] != -1) {
                c = decoder->M_VtoC_GF[t][k];
                // GF multiplication
                c = table->MULDEC[c][code->matValue[node][t]];
                decoder->M_VtoC_GF[t][k] = c;
            }
        }
    }
    
    // Specialized elementary steps for row degree 3
    // Forward step
    ElementaryStep_nm(decoder->M_VtoC_LLR[0], decoder->M_VtoC_LLR[1],
                     decoder->M_VtoC_GF[0], decoder->M_VtoC_GF[1],
                     decoder->M_CtoV_LLR[2], decoder->M_CtoV_GF[2],
                     table->ADDGF, code->GF, n_vc, n_vc, n_cv, num_operations);
    
    // Backward step
    ElementaryStep_nm(decoder->M_VtoC_LLR[1], decoder->M_VtoC_LLR[2],
                     decoder->M_VtoC_GF[1], decoder->M_VtoC_GF[2],
                     decoder->M_CtoV_LLR[0], decoder->M_CtoV_GF[0],
                     table->ADDGF, code->GF, n_vc, n_vc, n_cv, num_operations);
    
    // Merge step
    ElementaryStep_nm(decoder->M_VtoC_LLR[0], decoder->M_VtoC_LLR[2],
                     decoder->M_VtoC_GF[0], decoder->M_VtoC_GF[2],
                     decoder->M_CtoV_LLR[1], decoder->M_CtoV_GF[1],
                     table->ADDGF, code->GF, n_vc, n_vc, n_cv, num_operations);
    
    // Apply inverse GF rotation to check-to-variable messages
    for (t = 0; t < code->rowDegree[node]; t++) {
        // Find the stopping point
        for (k = 0; k < n_cv; k++) {
            if (decoder->M_CtoV_GF[t][k] == -1) {
                Stp = k;
                break;
            } else {
                Stp = n_cv;
            }
        }
        
        // GF division
        for (k = 0; k < Stp; k++) {
            decoder->M_CtoV_GF[t][k] = table->DIVDEC[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }
        
        // Apply offset to all LLR values
        for (k = 0; k < code->GF; k++) {
            llr_tmp[k] = decoder->M_CtoV_LLR[t][Stp - 1] + offset;
        }
        
        // Update specific LLR values
        for (k = 0; k < Stp; k++) {
            llr_tmp[decoder->M_CtoV_GF[t][k]] = decoder->M_CtoV_LLR[t][k];
        }
        
        // Update all CtoV messages
        for (k = 0; k < code->GF; k++) {
            decoder->M_CtoV_GF[t][k] = k;
            decoder->M_CtoV_LLR[t][k] = llr_tmp[k];
        }
    }
}