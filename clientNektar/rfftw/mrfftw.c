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

#include "rfftw.h"

rfftw_plan rfftw_create_plan(int n, fftw_direction dir,
			     int flags, rfftw_type type)
{
    rfftw_plan p;

    if ((n & 1) || n <= 0)
	return 0;		/* n must be even for 
				   real-complex FFT's! 
				   (and n must be > 0 for fftw) */

    p = fftw_malloc(sizeof(rfftw_plan_struct));
    if (!p)
	return 0;

    p->which_case = 0;

    if (dir == -1)
	p->which_case += 1;

    if (flags & FFTW_IN_PLACE)
	p->which_case += 2;

    if (type == COMPLEX_TO_REAL)
	p->which_case += 4;

    p->tw = 0;
    p->plan = 0;    /* initialize to NULL in case an error occurs */

    /* p->tw contains an array of the roots of unity:
       p->tw->tw_array[k] = exp(-2*pi*i * k/n) */
    p->tw = fftw_create_twiddle(n,2,n);
    if (!p->tw) {
	rfftw_destroy_plan(p);
	return 0;
    }
    p->plan = fftw_create_plan(n / 2, dir, flags);
    if (!p->plan) {
	rfftw_destroy_plan(p);
	return 0;
    }
    return p;
}

void rfftw_destroy_plan(rfftw_plan plan)
{
    if (plan) {
	if (plan->tw)
	    fftw_destroy_twiddle(plan->tw);
	if (plan->plan)
	    fftw_destroy_plan(plan->plan);
	fftw_free(plan);
    }
}

void rfftw(rfftw_plan plan, int howmany,
	   FFTW_COMPLEX * in, int istride, int idist,
	   FFTW_COMPLEX * out, int ostride, int odist)
{
    int iter, k;
    int n;
    FFTW_COMPLEX *omega;
    FFTW_REAL omega_re, omega_im;
    FFTW_REAL even_re, even_im, odd_re, odd_im;

#ifdef SCALING
    FFTW_REAL scal = 0.5/plan->plan->n;
#endif
    /* cache some values in local variables: */

    n = plan->plan->n;		/* 1/2 the size of the real array */
    omega = plan->tw->twarray;	/* roots of unity */

    /* Note: we put the "howmany" loop here, rather than using the
       howmany param. of fftw, because we have to do computations
       with the array after transforming it and we want to maximize
       locality: */

    /* Note also that we have this big "switch" statement so that
       we can pull as many "if" statements and multiplications out
       of the loop as possible. */

    switch (plan->which_case) {
    case 0:			/* REAL_TO_COMPLEX, FFTW_OUT_OF_PLACE, FFTW_BACKWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    fftw(plan->plan, 1,
		 in + iter * idist, istride, 0,
		 out + iter * odist, ostride, 0);

	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(out[k * ostride + iter * odist]) +
		    c_re(out[(n - k) * ostride + iter * odist]);
		odd_re = c_re(out[k * ostride + iter * odist]) -
		    c_re(out[(n - k) * ostride + iter * odist]);
		even_im = c_im(out[k * ostride + iter * odist]) -
		    c_im(out[(n - k) * ostride + iter * odist]);
		odd_im = c_im(out[k * ostride + iter * odist]) +
		    c_im(out[(n - k) * ostride + iter * odist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(out[k * ostride + iter * odist]) =
		    0.5 * (even_re +
			   odd_im * omega_re
			   - odd_re * omega_im);
		c_re(out[(n - k) * ostride + iter * odist]) =
		    0.5 * (even_re -
			   odd_im * omega_re
			   + odd_re * omega_im);

		c_im(out[k * ostride + iter * odist]) =
		    0.5 * (even_im -
			   odd_im * omega_im
			   - odd_re * omega_re);
		c_im(out[(n - k) * ostride + iter * odist]) =
		    0.5 * (-odd_im * omega_im -
			   even_im
			   - odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover 
				   k=n/2 when n is even */
		c_im(out[(n >> 1) * ostride + iter * odist]) *=
		    -c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
	    even_re =
		c_re(out[iter * odist]) + c_im(out[iter * odist]);
#ifndef COMPACT
	    c_re(out[n * ostride + iter * odist]) =
		c_re(out[iter * odist]) - c_im(out[iter * odist]);
	    c_im(out[n * ostride + iter * odist]) = 0.0;
	    c_re(out[iter * odist]) = even_re;
	    c_im(out[iter * odist]) = 0.0;
#else
	    c_im(out[iter * odist]) =
	        c_re(out[iter * odist]) - c_im(out[iter * odist]);
	    c_re(out[iter * odist]) = even_re;
#endif
	}
	break;

    case 1:			/* REAL_TO_COMPLEX, FFTW_OUT_OF_PLACE, FFTW_FORWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    fftw(plan->plan, 1, in + iter * idist, istride, 0,
		 out + iter * odist, ostride, 0);

	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(out[k * ostride + iter * odist]) +
		    c_re(out[(n - k) * ostride + iter * odist]);
		odd_re = c_re(out[k * ostride + iter * odist]) -
		    c_re(out[(n - k) * ostride + iter * odist]);
		even_im = c_im(out[k * ostride + iter * odist]) -
		    c_im(out[(n - k) * ostride + iter * odist]);
		odd_im = c_im(out[k * ostride + iter * odist]) +
		    c_im(out[(n - k) * ostride + iter * odist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(out[k * ostride + iter * odist]) =
		    0.5 * (even_re +
			   odd_im * omega_re
			   + odd_re * omega_im);
		c_re(out[(n - k) * ostride + iter * odist]) =
		    0.5 * (even_re -
			   odd_im * omega_re
			   - odd_re * omega_im);

		c_im(out[k * ostride + iter * odist]) =
		    0.5 * (even_im +
			   odd_im * omega_im
			   - odd_re * omega_re);
		c_im(out[(n - k) * ostride + iter * odist]) =
		    0.5 * (odd_im * omega_im -
			   even_im
			   - odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover 
				   k=n/2 when n is even */
		c_im(out[(n >> 1) * ostride + iter * odist]) *=
		    c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
	    even_re =
		c_re(out[iter * odist]) + c_im(out[iter * odist]);
#ifndef COMPACT
	    c_re(out[n * ostride + iter * odist]) =
		c_re(out[iter * odist]) - c_im(out[iter * odist]);
	    c_im(out[n * ostride + iter * odist]) = 0.0;
	    c_re(out[iter * odist]) = even_re;
	    c_im(out[iter * odist]) = 0.0;
#else
	    c_im(out[iter * odist]) =
	        c_re(out[iter * odist]) - c_im(out[iter * odist]);
	    c_re(out[iter * odist]) = even_re;
#endif
	}
	break;

    case 2:			/* REAL_TO_COMPLEX, FFTW_IN_PLACE, FFTW_BACKWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    fftw(plan->plan, 1, in + iter * idist, istride, 0, out, 1, 0);

	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(in[k * istride + iter * idist]) +
		    c_re(in[(n - k) * istride + iter * idist]);
		odd_re = c_re(in[k * istride + iter * idist]) -
		    c_re(in[(n - k) * istride + iter * idist]);
		even_im = c_im(in[k * istride + iter * idist]) -
		    c_im(in[(n - k) * istride + iter * idist]);
		odd_im = c_im(in[k * istride + iter * idist]) +
		    c_im(in[(n - k) * istride + iter * idist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(in[k * istride + iter * idist]) =
		    0.5 * (even_re +
			   odd_im * omega_re
			   - odd_re * omega_im);
		c_re(in[(n - k) * istride + iter * idist]) =
		    0.5 * (even_re -
			   odd_im * omega_re
			   + odd_re * omega_im);

		c_im(in[k * istride + iter * idist]) =
		    0.5 * (even_im -
			   odd_im * omega_im
			   - odd_re * omega_re);
		c_im(in[(n - k) * istride + iter * idist]) =
		    0.5 * (-odd_im * omega_im -
			   even_im
			   - odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover k=n/2 
				   when n is even */
		c_im(in[(n >> 1) * istride + iter * idist]) *=
		    -c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
	    even_re =
		c_re(in[iter * idist]) + c_im(in[iter * idist]);
#ifndef COMPACT
	    c_re(in[n * istride + iter * idist]) =
		c_re(in[iter * idist]) - c_im(in[iter * idist]);
	    c_im(in[n * istride + iter * idist]) = 0.0;
	    c_re(in[iter * idist]) = even_re;
	    c_im(in[iter * idist]) = 0.0;
#else
	    c_im(in[iter * idist]) =
	        c_re(in[iter * idist]) - c_im(in[iter * idist]);
	    c_re(in[iter * idist]) = even_re;
#endif
	}
	break;

    case 3:			/* REAL_TO_COMPLEX, FFTW_IN_PLACE, FFTW_FORWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    fftw(plan->plan, 1, in + iter * idist, istride, 0, out, 1, 0);

	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(in[k * istride + iter * idist]) +
		    c_re(in[(n - k) * istride + iter * idist]);
		odd_re = c_re(in[k * istride + iter * idist]) -
		    c_re(in[(n - k) * istride + iter * idist]);
		even_im = c_im(in[k * istride + iter * idist]) -
		    c_im(in[(n - k) * istride + iter * idist]);
		odd_im = c_im(in[k * istride + iter * idist]) +
		    c_im(in[(n - k) * istride + iter * idist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(in[k * istride + iter * idist]) =
		    0.5 * (even_re +
			   odd_im * omega_re
			   + odd_re * omega_im);
		c_re(in[(n - k) * istride + iter * idist]) =
		    0.5 * (even_re -
			   odd_im * omega_re
			   - odd_re * omega_im);

		c_im(in[k * istride + iter * idist]) =
		    0.5 * (even_im +
			   odd_im * omega_im
			   - odd_re * omega_re);
		c_im(in[(n - k) * istride + iter * idist]) =
		    0.5 * (odd_im * omega_im -
			   even_im
			   - odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover k=n/2 
				   when n is even */
		c_im(in[(n >> 1) * istride + iter * idist]) *=
		    c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
	    even_re =
		c_re(in[iter * idist]) + c_im(in[iter * idist]);
#ifndef COMPACT
	    c_re(in[n * istride + iter * idist]) =
		c_re(in[iter * idist]) - c_im(in[iter * idist]);
	    c_im(in[n * istride + iter * idist]) = 0.0;
	    c_re(in[iter * idist]) = even_re;
	    c_im(in[iter * idist]) = 0.0;
#else
	    c_im(in[iter * idist]) =
	        c_re(in[iter * idist]) - c_im(in[iter * idist]);
	    c_re(in[iter * idist]) = even_re;
#endif
#ifdef SCALING
	    for (k = 0; k < n; k++) {
	      c_re(in[k * istride + iter * idist]) *= scal;
	      c_im(in[k * istride + iter * idist]) *= scal;
	    }
#endif
	}
	break;

/**************** Complex -> Real Transforms ******************/
	/* Here, we are inverting the real->complex transforms above.
	   The process is algebraically very similar.  However, some
	   signs are different and everything is multiplied by a factor
	   of two!  (The factor of two is to fix the normalization so
	   that rfft followed by its inverse will yield N times the
	   original array, where N is the number of real data points.
	   This way, the normalization is the same as if we had used
	   the complex-complex FFT.) (Actually, multiplying by two
	   saves us some multiplications rather than costing anything!) */

    case 4:			/* COMPLEX_TO_REAL, FFTW_OUT_OF_PLACE, FFTW_BACKWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(in[k * istride + iter * idist]) +
		    c_re(in[(n - k) * istride + iter * idist]);
		odd_re = c_re(in[k * istride + iter * idist]) -
		    c_re(in[(n - k) * istride + iter * idist]);
		even_im = c_im(in[k * istride + iter * idist]) -
		    c_im(in[(n - k) * istride + iter * idist]);
		odd_im = c_im(in[k * istride + iter * idist]) +
		    c_im(in[(n - k) * istride + iter * idist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(in[k * istride + iter * idist]) =
		    (even_re - odd_im * omega_re
		     + odd_re * omega_im);
		c_re(in[(n - k) * istride + iter * idist]) =
		    (even_re + odd_im * omega_re
		     - odd_re * omega_im);

		c_im(in[k * istride + iter * idist]) =
		    (even_im + odd_im * omega_im
		     + odd_re * omega_re);
		c_im(in[(n - k) * istride + iter * idist]) =
		    (odd_im * omega_im - even_im
		     + odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover k=n/2
				   when n is even */
		c_re(in[(n >> 1) * istride +
			iter * idist]) *= 2.0;
		c_im(in[(n >> 1) * istride + iter * idist]) *=
		    2.0 * c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
#ifndef COMPACT
	    even_re =
		(c_re(in[iter * idist]) +
		 c_re(in[n * istride + iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_re(in[n * istride + iter * idist]));
#else
	    even_re =
		(c_re(in[iter * idist]) +
		 c_im(in[iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_im(in[iter * idist]));
#endif
	    c_re(in[iter * idist]) = even_re;

	    fftw(plan->plan, 1, in + iter * idist, istride, 0,
		 out + iter * odist, ostride, 0);
	}
	break;

    case 5:			/* COMPLEX_TO_REAL, FFTW_OUT_OF_PLACE, FFTW_FORWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(in[k * istride + iter * idist]) +
		    c_re(in[(n - k) * istride + iter * idist]);
		odd_re = c_re(in[k * istride + iter * idist]) -
		    c_re(in[(n - k) * istride + iter * idist]);
		even_im = c_im(in[k * istride + iter * idist]) -
		    c_im(in[(n - k) * istride + iter * idist]);
		odd_im = c_im(in[k * istride + iter * idist]) +
		    c_im(in[(n - k) * istride + iter * idist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(in[k * istride + iter * idist]) =
		    (even_re - odd_im * omega_re
		     - odd_re * omega_im);
		c_re(in[(n - k) * istride + iter * idist]) =
		    (even_re + odd_im * omega_re
		     + odd_re * omega_im);

		c_im(in[k * istride + iter * idist]) =
		    (even_im - odd_im * omega_im
		     + odd_re * omega_re);
		c_im(in[(n - k) * istride + iter * idist]) =
		    (-odd_im * omega_im - even_im
		     + odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover k=n/2
				   when n is even */
		c_re(in[(n >> 1) * istride +
			iter * idist]) *= 2.0;
		c_im(in[(n >> 1) * istride + iter * idist]) *=
		    -2.0 * c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
#ifndef COMPACT
	    even_re =
		(c_re(in[iter * idist]) +
		 c_re(in[n * istride + iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_re(in[n * istride + iter * idist]));
#else
	    even_re =
		(c_re(in[iter * idist]) +
		 c_im(in[iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_im(in[iter * idist]));
#endif
	    c_re(in[iter * idist]) = even_re;

	    fftw(plan->plan, 1, in + iter * idist, istride, 0,
		 out + iter * odist, ostride, 0);
	}
	break;

    case 6:			/* COMPLEX_TO_REAL, FFTW_IN_PLACE, FFTW_BACKWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(in[k * istride + iter * idist]) +
		    c_re(in[(n - k) * istride + iter * idist]);
		odd_re = c_re(in[k * istride + iter * idist]) -
		    c_re(in[(n - k) * istride + iter * idist]);
		even_im = c_im(in[k * istride + iter * idist]) -
		    c_im(in[(n - k) * istride + iter * idist]);
		odd_im = c_im(in[k * istride + iter * idist]) +
		    c_im(in[(n - k) * istride + iter * idist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(in[k * istride + iter * idist]) =
		    (even_re - odd_im * omega_re
		     + odd_re * omega_im);
		c_re(in[(n - k) * istride + iter * idist]) =
		    (even_re + odd_im * omega_re
		     - odd_re * omega_im);

		c_im(in[k * istride + iter * idist]) =
		    (even_im + odd_im * omega_im
		     + odd_re * omega_re);
		c_im(in[(n - k) * istride + iter * idist]) =
		    (odd_im * omega_im - even_im
		     + odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover k=n/2 
				   when n is even */
		c_re(in[(n >> 1) * istride +
			iter * idist]) *= 2.0;
		c_im(in[(n >> 1) * istride + iter * idist]) *=
		    2.0 * c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
#ifndef COMPACT
	    even_re =
		(c_re(in[iter * idist]) +
		 c_re(in[n * istride + iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_re(in[n * istride + iter * idist]));
#else
	    even_re =
		(c_re(in[iter * idist]) +
		 c_im(in[iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_im(in[iter * idist]));
#endif
	    c_re(in[iter * idist]) = even_re;

	    fftw(plan->plan, 1, in + iter * idist, istride, 0, out, 1, 0);
	}
	break;

    case 7:			/* COMPLEX_TO_REAL, FFTW_IN_PLACE, FFTW_FORWARD */
	for (iter = 0; iter < howmany; ++iter) {
	    for (k = 1; k < ((n + 1) >> 1); ++k) {
		even_re = c_re(in[k * istride + iter * idist]) +
		    c_re(in[(n - k) * istride + iter * idist]);
		odd_re = c_re(in[k * istride + iter * idist]) -
		    c_re(in[(n - k) * istride + iter * idist]);
		even_im = c_im(in[k * istride + iter * idist]) -
		    c_im(in[(n - k) * istride + iter * idist]);
		odd_im = c_im(in[k * istride + iter * idist]) +
		    c_im(in[(n - k) * istride + iter * idist]);

		omega_re = c_re(omega[k]);
		omega_im = c_im(omega[k]);

		c_re(in[k * istride + iter * idist]) =
		    (even_re - odd_im * omega_re
		     - odd_re * omega_im);
		c_re(in[(n - k) * istride + iter * idist]) =
		    (even_re + odd_im * omega_re
		     + odd_re * omega_im);

		c_im(in[k * istride + iter * idist]) =
		    (even_im - odd_im * omega_im
		     + odd_re * omega_re);
		c_im(in[(n - k) * istride + iter * idist]) =
		    (-odd_im * omega_im - even_im
		     + odd_re * omega_re);
	    }

	    if (!(n & 1)) {	/* handle leftover k=n/2 
				   when n is even */
		c_re(in[(n >> 1) * istride
			+ iter * idist]) *= 2.0;
		c_im(in[(n >> 1) * istride + iter * idist]) *=
		    -2.0 * c_im(omega[n >> 1]);
	    }
	    /* Do k = 0 and k = n/2 cases: */
#ifndef COMPACT
	    even_re =
		(c_re(in[iter * idist]) +
		 c_re(in[n * istride + iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_re(in[n * istride + iter * idist]));
#else
	    even_re =
		(c_re(in[iter * idist]) +
		 c_im(in[iter * idist]));
	    c_im(in[iter * idist]) =
		(c_re(in[iter * idist]) -
		 c_im(in[iter * idist]));
#endif
	    c_re(in[iter * idist]) = even_re;

	    fftw(plan->plan, 1, in + iter * idist, istride, 0, out, 1, 0);
	}
	break;

    default:
	fftw_die("Impossible plan in rfftw!\n");
    }
}
