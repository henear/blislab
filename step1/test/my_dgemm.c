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
  
  // Early return if possible
  if ( m == 0 || n == 0 || k == 0 ) {
    printf( "bl_dgemm(): early return\n" );
    return;
  }
  double *cp;
  double *cp1;
  double *cp2;
  double *cp3;
                     
       for ( j = 0; j < n; j +=1 ) {
        for ( i = 0; i < m; i +=4 ) {                    	
	  
	  cp = &C(i, j); 
	  cp1 = &C(i+1, j);
	  cp2 = &C(i+2, j);
          cp3 = &C(i+3, j);
	 
	 
          for ( p = 0; p < k; p +=8 ) {            
              //C[ j * ldc + i ] += A[ p * lda + i ] * B[ j * ldb + p ];
	      //Each operand is a MACRO defined in bl_dgemm() function.
              *cp+= A( i, p ) * B( p, j ); 	
	      *cp+= A( i, p+1)* B(p+1,j );
	      *cp+= A( i, p+2)* B(p+2,j );
	      *cp+= A( i, p+3)* B(p+3,j ); 
	      *cp+= A( i, p+4)* B(p+4,j );
	      *cp+= A( i, p+5)* B(p+5,j );
	      *cp+= A( i, p+6)* B(p+6,j ); 
	      *cp+= A( i, p+7)* B(p+7,j ); 

 	      *cp1+= A( i+1, p ) * B( p, j ); 	
	      *cp1+= A( i+1, p+1)* B(p+1,j );
	      *cp1+= A( i+1, p+2)* B(p+2,j );
	      *cp1+= A( i+1, p+3)* B(p+3,j ); 
	      *cp1+= A( i+1, p+4)* B(p+4,j );
	      *cp1+= A( i+1, p+5)* B(p+5,j );
	      *cp1+= A( i+1, p+6)* B(p+6,j ); 
	      *cp1+= A( i+1, p+7)* B(p+7,j );

              *cp2+= A( i+2, p ) * B( p, j ); 	
	      *cp2+= A( i+2, p+1)* B(p+1,j );
	      *cp2+= A( i+2, p+2)* B(p+2,j );
	      *cp2+= A( i+2, p+3)* B(p+3,j ); 
	      *cp2+= A( i+2, p+4)* B(p+4,j );
	      *cp2+= A( i+2, p+5)* B(p+5,j );
	      *cp2+= A( i+2, p+6)* B(p+6,j ); 
	      *cp2+= A( i+2, p+7)* B(p+7,j );
 
	      *cp3+= A( i+3, p ) * B( p, j ); 	
	      *cp3+= A( i+3, p+1)* B(p+1,j );
	      *cp3+= A( i+3, p+2)* B(p+2,j );
	      *cp3+= A( i+3, p+3)* B(p+3,j ); 
	      *cp3+= A( i+3, p+4)* B(p+4,j );
	      *cp3+= A( i+3, p+5)* B(p+5,j );
	      *cp3+= A( i+3, p+6)* B(p+6,j ); 
	      *cp3+= A( i+3, p+7)* B(p+7,j );  

 	      
        }
                                         
      }                                          
  }                                              

}


