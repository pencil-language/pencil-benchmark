/*
 * Copyright (c) 2013-2014, ARM Limited
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* Some macros to measure wall clock time.
 *
 * Usage:
 *
 *   START_MEASURE_TIME
 *   ... do the work ...
 *   STOP_MEASURE_TIME
 *
 *   PRINT_TIME
 */

#ifndef MEASURE_TIME_H
#define MEASURE_TIME_H

#define START_MEASURE_TIME prl_prof_reset(); prl_prof_start();
#define STOP_MEASURE_TIME prl_prof_stop();
#define PRINT_TIME prl_prof_dump();

#include <sys/time.h>

#define USEC_PER_SEC 1e6

#define START_MEASURE_TIME_HOST struct timeval time_start_; \
	                           gettimeofday(&time_start_, NULL);

#define STOP_MEASURE_TIME_HOST  struct timeval time_stop_; \
	                           gettimeofday(&time_stop_, NULL);

#define PRINT_TIME_HOST printf("Duration of reference: %lf\n", (double)(time_stop_.tv_sec - time_start_.tv_sec) + ((double)(time_stop_.tv_usec - time_start_.tv_usec))/USEC_PER_SEC);

#endif
