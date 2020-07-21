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

/* Generated by $Id: fn_4.c,v 1.1 2004/09/26 10:32:35 ssherw Exp $ */

/* This function contains 16 FP additions and 0 FP multiplications */

void fftw_no_twiddle_4(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[2 * istride]);
	  tim1_1_0 = c_im(in[2 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1;
     c_im(out[0]) = tim0_0_0 + tim0_0_1;
     c_re(out[2 * ostride]) = tre0_0_0 - tre0_0_1;
     c_im(out[2 * ostride]) = tim0_0_0 - tim0_0_1;
     c_re(out[ostride]) = tre0_1_0 + tim0_1_1;
     c_im(out[ostride]) = tim0_1_0 - tre0_1_1;
     c_re(out[3 * ostride]) = tre0_1_0 - tim0_1_1;
     c_im(out[3 * ostride]) = tim0_1_0 + tre0_1_1;
}
