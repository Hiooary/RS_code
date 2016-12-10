#ifndef _RS_H
#define _RS_H

#include <math.h>
#include <stdio.h>

#define mm 4 /* RS code over GF(2**8) - change to suit */
#define nn 15 /* nn=2**mm -1 length of codeword */
#define tt 2 /* number of errors that can be corrected (可以校正的错误数)*/ //lq
#define kk 11 /* kk = nn-2*tt */  //lq





void encode_rs(int recd[nn], int data[kk], int  bb[nn-kk]);
void decode_rs(int recd[nn], int dataout[kk]);
void encode_rs_8(unsigned char recd[nn], unsigned char data[kk], unsigned char  bb[nn-kk]);
void decode_rs_8(unsigned char recd[nn], unsigned char dataout[kk]);
const int pp[mm+1] = { 1, 1, 0, 0, 1} ; /* specify irreducible polynomial coeffts */
const int alpha_to[nn+1] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9, 0};
const int gg[nn-kk+1] = {10, 3, 6, 13, 0};
const int index_of[nn+1] = {-1, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12};

#endif