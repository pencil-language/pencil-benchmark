/*
 * File:   math.h
 * Author: Robert David
 *
 * Experimental PENCIL standard library header.
 * Created on June 2, 2014, 10:21 AM
 */

#ifndef __PENCIL_MATH_H__
#define __PENCIL_MATH_H__

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef __PENCIL__
    /* Almost all of the functions are available as OpenCL builtins */

#define __PENCIL_IMPL_DEFINE_OVERLOADS_1_PARAM(functionname) \
    float  __attribute((overloadable)) functionname(float );   \
    double __attribute((overloadable)) functionname(double);   \

#define __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS(functionname) \
    float  __attribute((overloadable)) functionname(float , float );   \
    double __attribute((overloadable)) functionname(double, double);   \

#define __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS_2ND_POINTER(functionname) \
    float  __attribute((overloadable)) functionname(float , float *);   \
    double __attribute((overloadable)) functionname(double, double*);   \

#define __PENCIL_IMPL_DEFINE_OVERLOADS_3_PARAMS(functionname) \
    float  __attribute((overloadable)) functionname(float , float , float );   \
    double __attribute((overloadable)) functionname(double, double, double);   \

#define __PENCIL_IMPL_DEFINE_INTEGER_OVERLOADS_2_PARAMS(functionname) \
      signed     short int __attribute((overloadable)) functionname(  signed short int,   signed short int); \
      signed           int __attribute((overloadable)) functionname(  signed       int,   signed       int); \
      signed      long int __attribute((overloadable)) functionname(  signed  long int,   signed  long int); \
    unsigned     short int __attribute((overloadable)) functionname(unsigned short int, unsigned short int); \
    unsigned           int __attribute((overloadable)) functionname(unsigned       int, unsigned       int); \
    unsigned      long int __attribute((overloadable)) functionname(unsigned  long int, unsigned  long int); \
             char __attribute((overloadable)) functionname(         char,          char); \
      signed char __attribute((overloadable)) functionname(  signed char,   signed char); \
    unsigned char __attribute((overloadable)) functionname(unsigned char, unsigned char); \

#define __PENCIL_IMPL_DEFINE_INTEGER_OVERLOADS_3_PARAMS(functionname) \
    signed     short int __attribute((overloadable)) functionname(  signed short int,   signed short int,   signed short int); \
    signed           int __attribute((overloadable)) functionname(  signed       int,   signed       int,   signed       int); \
    signed      long int __attribute((overloadable)) functionname(  signed  long int,   signed  long int,   signed  long int); \
  unsigned     short int __attribute((overloadable)) functionname(unsigned short int, unsigned short int, unsigned short int); \
  unsigned           int __attribute((overloadable)) functionname(unsigned       int, unsigned       int, unsigned       int); \
  unsigned      long int __attribute((overloadable)) functionname(unsigned  long int, unsigned  long int, unsigned  long int); \
           char __attribute((overloadable)) functionname(         char,          char,          char); \
    signed char __attribute((overloadable)) functionname(  signed char,   signed char,   signed char); \
  unsigned char __attribute((overloadable)) functionname(unsigned char, unsigned char, unsigned char); \

    __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS(min)
    __PENCIL_IMPL_DEFINE_INTEGER_OVERLOADS_2_PARAMS(min)
    __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS(max)
    __PENCIL_IMPL_DEFINE_INTEGER_OVERLOADS_2_PARAMS(max)
    __PENCIL_IMPL_DEFINE_OVERLOADS_3_PARAMS(clamp)
    __PENCIL_IMPL_DEFINE_INTEGER_OVERLOADS_3_PARAMS(clamp)

    __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS(atan)
    __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS(atan2)
    __PENCIL_IMPL_DEFINE_OVERLOADS_1_PARAM(tan)
    __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS(hypot)

    __PENCIL_IMPL_DEFINE_OVERLOADS_1_PARAM(exp)

    __PENCIL_IMPL_DEFINE_OVERLOADS_1_PARAM(ceil)
    __PENCIL_IMPL_DEFINE_OVERLOADS_1_PARAM(floor)
    __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS_2ND_POINTER(fract)

    /* TODO: add definitions for more OpenCL builtins */

#undef __PENCIL_IMPL_DEFINE_OVERLOADS_1_PARAM
#undef __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS
#undef __PENCIL_IMPL_DEFINE_OVERLOADS_3_PARAMS
#undef __PENCIL_IMPL_DEFINE_OVERLOADS_2_PARAMS_2ND_POINTER
#undef __PENCIL_IMPL_DEFINE_INTEGER_OVERLOADS_2_PARAMS

    /* A few C standard library math functions are not available as OpenCL builtin */
    /* TODO:
     *  div/ldiv/lldiv,
     *  fpclassify,
     *  lrint, llrint, lround, llround,
     *  nearbyint, nexttoward,
     *  scalbn, scalbln
     */

#else
    /* Definitions for functions to be used when compiled as C code */

    /* Lot of stuff defined as type-generic macros */
#include <tgmath.h>

    /* Some functions can be defined with a simple equation */
            
    /* This solution uses a GCC extension to avoid multiple evaluation of macro parameters.
     * It wouls be better to use C11 _Generic macros, but that's only supported from GCC 4.9 and nVidia OpenCL only works with GCC 4.6 */
#define min(a,b) ({ \
    __typeof__(a) _a_temp_ = (a); \
    __typeof__(b) _b_temp_ = (b); \
    _a_temp_ <= _b_temp_ ? _a_temp_ : _b_temp_; \
    })
#define max(a,b) ({ \
    __typeof__(a) _a_temp_ = (a); \
    __typeof__(b) _b_temp_ = (b); \
    _a_temp_ >= _b_temp_ ? _a_temp_ : _b_temp_; \
    })
#define clamp(val, min, max) ({ \
    __typeof__(val) _val_temp_ = (val); \
    __typeof__(min) _min_temp_ = (min); \
    __typeof__(max) _max_temp_ = (max); \
    (_val_temp_ < _min_temp_) ? _min_temp_ : (_val_temp_ > _max_temp_) ? _max_temp_ : _val_temp_; \
    })
#define fract(a, b) ({ \
    __typeof__(a) _a_temp_ = (a); \
    __typeof__(b) _b_temp_ = (b); \
    *_b_temp_ = floor(_a_temp_); \
    fmin( _a_temp_ - *_b_temp_, nexttoward((__typeof__(a))1.0, (__typeof__(a))0.0 )); \
    })

    /* TODO:
     *  acospi, asinpi, atan2pi, atanpi, cospi, sinpi, tanpi,
     *  exp10, pown, powr, rootn,
     *  rsqrt,
     *  degrees, radians,
     *  isequal, isnotequal, isordered,
     *  lgamma_r,
     *  select, bitselect,
     *  smoothstep, step,
     *
     * ... (rest don't have clear formulas, figure out what to do)
     *  */
#endif

    /* Some added functions */

#ifdef	__cplusplus
}
#endif

#endif	/* __PENCIL_MATH_H__
*/

