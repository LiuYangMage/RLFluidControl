/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 */

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: ftwi_9.c,v 1.1 2004/09/26 10:32:35 ssherw Exp $ */

/* This function contains 108 FP additions and 72 FP multiplications */

void fftwi_twiddle_9(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 8) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_2_0 = (tr * twr) + (ti * twi);
		    tim1_2_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    tre0_1_0 = tre2_0_0 + tre2_1_0;
		    tre0_2_0 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    tim0_1_0 = tim2_0_0 + tim2_1_0;
		    tim0_2_0 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_2_0 = (tr * twr) + (ti * twi);
		    tim1_2_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    tre0_1_1 = tre2_0_0 + tre2_1_0;
		    tre0_2_1 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    tim0_1_1 = tim2_0_0 + tim2_1_0;
		    tim0_2_1 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[8 * stride]);
		    ti = c_im(inout[8 * stride]);
		    twr = c_re(W[7]);
		    twi = c_im(W[7]);
		    tre1_2_0 = (tr * twr) + (ti * twi);
		    tim1_2_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    tre0_1_2 = tre2_0_0 + tre2_1_0;
		    tre0_2_2 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    tim0_1_2 = tim2_0_0 + tim2_1_0;
		    tim0_2_2 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_2 - tim0_0_1);
	       c_re(inout[3 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[6 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_1 - tre0_0_2);
	       c_im(inout[3 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[6 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tre0_1_1) - (((FFTW_REAL) FFTW_K642787609) * tim0_1_1);
	       tim1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tim0_1_1) + (((FFTW_REAL) FFTW_K642787609) * tre0_1_1);
	       tre1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_1_2) - (((FFTW_REAL) FFTW_K984807753) * tim0_1_2);
	       tim1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_1_2) + (((FFTW_REAL) FFTW_K984807753) * tre0_1_2);
	       c_re(inout[stride]) = tre0_1_0 + tre1_1_0 + tre1_2_0;
	       c_im(inout[stride]) = tim0_1_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_1_0;
		    c_re(inout[7 * stride]) = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_1_0;
		    c_im(inout[7 * stride]) = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_2_1) - (((FFTW_REAL) FFTW_K984807753) * tim0_2_1);
	       tim1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_2_1) + (((FFTW_REAL) FFTW_K984807753) * tre0_2_1);
	       tre1_2_0 = (((FFTW_REAL) FFTW_K939692620) * tre0_2_2) + (((FFTW_REAL) FFTW_K342020143) * tim0_2_2);
	       tim1_2_0 = (((FFTW_REAL) FFTW_K342020143) * tre0_2_2) - (((FFTW_REAL) FFTW_K939692620) * tim0_2_2);
	       c_re(inout[2 * stride]) = tre0_2_0 + tre1_1_0 - tre1_2_0;
	       c_im(inout[2 * stride]) = tim0_2_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K499999999) * (tre1_2_0 - tre1_1_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_1_0;
		    c_re(inout[8 * stride]) = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim0_2_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 + tre1_2_0);
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_1_0;
		    c_im(inout[8 * stride]) = tim2_0_0 - tim2_1_0;
	       }
	  }
     }
}
