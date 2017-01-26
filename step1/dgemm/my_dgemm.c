/*
 * --------------------------------------------------------------------------
 * BLISLAB 
 * --------------------------------------------------------------------------
 * Copyright (C) 2016, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * bl_dgemm.c
 *
 *
 * Purpose:
 * this is the main file of blislab dgemm.
 *
 * Todo:
 *
 *
 * Modification:
 *
 * 
 * */
 

#include "bl_dgemm.h"


void bl_dgemm(
    int    m,
    int    n,
    int    k,
    double *A,
    int    lda,
    double *B,
    int    ldb,
    double *C,        // must be aligned
    int    ldc        // ldc must also be aligned
)
{
  int    i, j, p;
  double *cp;
  double *cp1;
  double *cp2;
  double *cp3;
  double mR;
  double nR;

  // Early return if possible
  if ( m == 0 || n == 0 || k == 0 ) {
    printf( "bl_dgemm(): early return\n" );
    return;
  }

  for ( j = 0; j < n; j +=4 ) {                    // 2-th loop

      cp = &C [ j*ldc];
      cp1 = &C [ j+1*ldc];
      cp2 = &C [ j+2*ldc];
      cp3 = &C [ j+3*ldc];

      for ( i = 0; i < m; i += 4 ) {                // 1-th loop
 
          
          for ( p = 0; p < k; p +=4 ) {            // 0-th loop

              //C[ j * ldc + i ] += A[ p * lda + i ] * B[ j * ldb + p ];

              *(cp+0) = A( i, p ) * B( p, j );
              *(cp+0) = A( i, p+1 ) * B( p+1, j ); 
              *(cp+0) = A( i, p+2 ) * B( p+2, j );
              *(cp+0) = A( i, p+3 ) * B( p+3, j );

              *(cp+1) = A( i+1, p ) * B( p, j );
              *(cp+1) = A( i+1, p+1 ) * B( p+1, j ); 
              *(cp+1) = A( i+1, p+2 ) * B( p+2, j );
              *(cp+1) = A( i+1, p+3 ) * B( p+3, j );

              *(cp+2) = A( i+2, p ) * B( p, j );
              *(cp+2) = A( i+2, p+1 ) * B( p+1, j ); 
              *(cp+2) = A( i+2, p+2 ) * B( p+2, j );
              *(cp+2) = A( i+2, p+3 ) * B( p+3, j );
 
              *(cp+3) = A( i+3, p ) * B( p, j );
              *(cp+3) = A( i+3, p+1 ) * B( p+1, j ); 
              *(cp+3) = A( i+3, p+2 ) * B( p+2, j );
              *(cp+3) = A( i+3, p+3 ) * B( p+3, j );
              cp += 4;   


              *(cp1+0) = A( i, p ) * B( p, j );
              *(cp1+0) = A( i, p+1 ) * B( p+1, j ); 
              *(cp1+0) = A( i, p+2 ) * B( p+2, j );
              *(cp1+0) = A( i, p+3 ) * B( p+3, j );

              *(cp1+1) = A( i+1, p ) * B( p, j );
              *(cp1+1) = A( i+1, p+1 ) * B( p+1, j ); 
              *(cp1+1) = A( i+1, p+2 ) * B( p+2, j );
              *(cp1+1) = A( i+1, p+3 ) * B( p+3, j );

              *(cp1+2) = A( i+2, p ) * B( p, j );
              *(cp1+2) = A( i+2, p+1 ) * B( p+1, j ); 
              *(cp1+2) = A( i+2, p+2 ) * B( p+2, j );
              *(cp1+2) = A( i+2, p+3 ) * B( p+3, j );
 
              *(cp1+3) = A( i+3, p ) * B( p, j );
              *(cp1+3) = A( i+3, p+1 ) * B( p+1, j ); 
              *(cp1+3) = A( i+3, p+2 ) * B( p+2, j );
              *(cp1+3) = A( i+3, p+3 ) * B( p+3, j );
              cp1 += 4;

              *(cp2+0) = A( i, p ) * B( p, j );
              *(cp2+0) = A( i, p+1 ) * B( p+1, j ); 
              *(cp2+0) = A( i, p+2 ) * B( p+2, j );
              *(cp+0) = A( i, p+3 ) * B( p+3, j );

              *(cp2+1) = A( i+1, p ) * B( p, j );
              *(cp2+1) = A( i+1, p+1 ) * B( p+1, j ); 
              *(cp2+1) = A( i+1, p+2 ) * B( p+2, j );
              *(cp2+1) = A( i+1, p+3 ) * B( p+3, j );

              *(cp2+2) = A( i+2, p ) * B( p, j );
              *(cp2+2) = A( i+2, p+1 ) * B( p+1, j ); 
              *(cp2+2) = A( i+2, p+2 ) * B( p+2, j );
              *(cp2+2) = A( i+2, p+3 ) * B( p+3, j );
 
              *(cp2+3) = A( i+3, p ) * B( p, j );
              *(cp2+3) = A( i+3, p+1 ) * B( p+1, j ); 
              *(cp2+3) = A( i+3, p+2 ) * B( p+2, j );
              *(cp2+3) = A( i+3, p+3 ) * B( p+3, j );
              cp2 += 4; 

              *(cp3+0) = A( i, p ) * B( p, j );
              *(cp3+0) = A( i, p+1 ) * B( p+1, j ); 
              *(cp3+0) = A( i, p+2 ) * B( p+2, j );
              *(cp3+0) = A( i, p+3 ) * B( p+3, j );

              *(cp3+1) = A( i+1, p ) * B( p, j );
              *(cp3+1) = A( i+1, p+1 ) * B( p+1, j ); 
              *(cp3+1) = A( i+1, p+2 ) * B( p+2, j );
              *(cp3+1) = A( i+1, p+3 ) * B( p+3, j );

              *(cp3+2) = A( i+2, p ) * B( p, j );
              *(cp3+2) = A( i+2, p+1 ) * B( p+1, j ); 
              *(cp3+2) = A( i+2, p+2 ) * B( p+2, j );
              *(cp3+2) = A( i+2, p+3 ) * B( p+3, j );
 
              *(cp3+3) = A( i+3, p ) * B( p, j );
              *(cp3+3) = A( i+3, p+1 ) * B( p+1, j ); 
              *(cp3+3) = A( i+3, p+2 ) * B( p+2, j );
              *(cp3+3) = A( i+3, p+3 ) * B( p+3, j );
              cp3 += 4;  
             //Each operand is a MACRO defined in bl_dgemm() function.
             // C( i, j ) += A( i,p+1) * B( p+1,j);
          }                                      // End 0-th loop
      }                                          // End 1-th loop
  }                                              // End 2-th loop

}



