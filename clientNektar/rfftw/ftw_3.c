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

/* Generated by $Id: ftw_3.c,v 1.1 2004/09/26 10:32:35 ssherw Exp $ */

/* This function contains 18 FP additions and 12 FP multiplications */

void fftw_twiddle_3(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 2) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) - (ti * twi);
	       tim0_1_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) - (ti * twi);
	       tim0_2_0 = (tr * twi) + (ti * twr);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_0 + tre0_2_0));
	       tre1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_0 - tim0_2_0);
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[2 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_0 + tim0_2_0));
	       tim1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_2_0 - tre0_1_0);
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[2 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
