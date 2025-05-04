#ifndef STRUCT_H_INCLUDED
#define STRUCT_H_INCLUDED

/*!
 * \file struct.h
 * \brief Header for data structures used in LDPC decoding
 */

//! \struct GF_element
//! \brief Represents an element in GF(q) with corresponding LLR
typedef struct {
    int GF;
    float LLR;
} GF_element;

//! \struct syndrome_type
//! \brief Represents a syndrome value with configuration info
typedef struct {
    int GF;
    float LLR;
    int config; // Configuration number
} syndrome_type;

//! \struct GF_syndrome_type
//! \brief Used in the decoding process to store minimum LLRs and their indices
typedef struct {
    float min1;
    float min2;
    float min3;
    float min4;
    int index1;
    int index2;
    int index3;
    float LLR;
} GF_syndrome_type;

/**********************************************************/
/**************** LDPC Code Parameters Description ********/
/**********************************************************/

/*!
 * \struct code_t
 * \brief LDPC code parameters with sparse matrix representation
 * 
 * This structure contains all necessary parameters for a given LDPC code.
 * These parameters are usually loaded from a text file that holds the 
 * parity-check matrix.
 */
typedef struct {
    int N;             // Number of columns in H
    int M;             // Number of rows in H
    int K;             // Number of information symbols: K = N - M
    int GF;            // Field order (e.g., GF=64 or GF=256)
    int logGF;         // logGF = log2(GF), i.e., the number of bits in each symbol
    int **mat;         // Parity-check matrix: 2D array listing VNs (columns) involved in each check node (row)
    int **matValue;    // Non-binary coefficients of the parity-check matrix (GF(q) values per constraint)
    int *rowDegree;    // Degree of each check node (row)
    int *columnDegree; // Degree of each variable node (column)
    int nbBranch;      // Number of edges in the Tanner graph
    float rate;        // Code rate
    int **matUT;       // Upper Triangular form used for encoding (after Gaussian elimination)
    int *Perm;         // Permutation array used in Gaussian elimination
} code_t;

/*!
 * \struct table_t
 * \brief Lookup tables for arithmetic in GF(q)
 */
typedef struct {
    int **BINGF;    // Binary representation of GF(q) symbols
    int **ADDGF;    // Addition table for GF(q)
    int **MULGF;    // Multiplication table for GF(q)
    int **DIVGF;    // Division table for GF(q)
    int *DECGF;     // Binary-to-decimal mapping of GF(q) symbols
    int *GFDEC;     // Decimal-to-GF(q) symbol mapping
    int **MULDEC;   // Multiplication in decimal domain
    int **DIVDEC;   // Division in decimal domain
} table_t;

/**********************************************************/
/****************** EMS Decoder Parameters ****************/
/**********************************************************/

typedef float softdata_t;

/*!
 * \struct decoder_t
 * \brief Contains all parameters and data for EMS decoding
 */
typedef struct {
    int N;
    int n_cv;          // Top-n_cv most reliable symbols passed from Check Node to Variable Node
    int n_vc;          // Top-n_vc most reliable symbols passed from Variable Node to Check Node
    int nbBranch;      // Number of edges in the Tanner graph
    softdata_t **CtoV; // LLRs from Check to Variable nodes (size: nbBranch × nbMax)
    softdata_t **VtoC; // LLRs from Variable to Check nodes (size: nbBranch × nbMax)
    softdata_t **intrinsic_LLR; // Channel LLRs (size: N × GF)
    int **intrinsic_GF;         // GF symbols associated with intrinsic LLRs (same size)
    softdata_t **APP;           // A Posteriori Probabilities (same size as intrinsic_LLR)
    softdata_t **M_VtoC_LLR;    // LLR inputs to Check node processor
    int **M_VtoC_GF;            // GF indices corresponding to M_VtoC_LLR
    softdata_t **M_CtoV_LLR;    // LLR outputs from Check node processor
    int **M_CtoV_GF;            // GF indices corresponding to M_CtoV_LLR
} decoder_t;

/*!
 * \brief Binary image of the GF(64) field
 * Primitive polynomial: P(x) = x^6 + x + 1
 */
static const int BinGF_64[64][6]=
{
    {0,0,0,0,0,0}, //0 => 0
    {1,0,0,0,0,0}, //32 => 1
    {0,1,0,0,0,0}, //16 => 2
    {0,0,1,0,0,0}, //8 => 3
    {0,0,0,1,0,0}, //4 => 4
    {0,0,0,0,1,0}, //2 => 5
    {0,0,0,0,0,1}, //1 => 6
    {1,1,0,0,0,0}, // 48 => 7
    {0,1,1,0,0,0}, //24 => 8
    {0,0,1,1,0,0},
    {0,0,0,1,1,0},
    {0,0,0,0,1,1},
    {1,1,0,0,0,1},
    {1,0,1,0,0,0},
    {0,1,0,1,0,0},
    {0,0,1,0,1,0},
    {0,0,0,1,0,1},
    {1,1,0,0,1,0},
    {0,1,1,0,0,1},
    {1,1,1,1,0,0},
    {0,1,1,1,1,0},
    {0,0,1,1,1,1},
    {1,1,0,1,1,1},
    {1,0,1,0,1,1},
    {1,0,0,1,0,1},
    {1,0,0,0,1,0},
    {0,1,0,0,0,1},
    {1,1,1,0,0,0},
    {0,1,1,1,0,0},
    {0,0,1,1,1,0},
    {0,0,0,1,1,1},
    {1,1,0,0,1,1},
    {1,0,1,0,0,1},
    {1,0,0,1,0,0},
    {0,1,0,0,1,0},
    {0,0,1,0,0,1},
    {1,1,0,1,0,0},
    {0,1,1,0,1,0},
    {0,0,1,1,0,1},
    {1,1,0,1,1,0},
    {0,1,1,0,1,1},
    {1,1,1,1,0,1},
    {1,0,1,1,1,0},
    {0,1,0,1,1,1},
    {1,1,1,0,1,1},
    {1,0,1,1,0,1},
    {1,0,0,1,1,0},
    {0,1,0,0,1,1},
    {1,1,1,0,0,1},
    {1,0,1,1,0,0},
    {0,1,0,1,1,0},
    {0,0,1,0,1,1},
    {1,1,0,1,0,1},
    {1,0,1,0,1,0},
    {0,1,0,1,0,1},
    {1,1,1,0,1,0},
    {0,1,1,1,0,1},
    {1,1,1,1,1,0},
    {0,1,1,1,1,1},
    {1,1,1,1,1,1},
    {1,0,1,1,1,1},
    {1,0,0,1,1,1},
    {1,0,0,0,1,1},
    {1,0,0,0,0,1}, //33 => 63
};
#endif