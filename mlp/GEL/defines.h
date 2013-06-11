/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

/*! \file defines.h
    \brief Provides common typedefs, and macros.
    \see GEL_EXPORT
*/

#ifndef CORE_DEFINES_H
#define CORE_DEFINES_H

#include <iostream>

/*!
 \def M_PI
 \brief Define Pi if needed.
*/
#ifndef M_PI
#  define M_PI 3.141592653589793238462643
#endif

/*!
 \def M_LN2
 \brief Define ln(2) if needed.
*/
#ifndef M_LN2
#  define M_LN2 0.69314718055994530942
#endif
/*!
 \def GEL_EXPORT
 \brief GEL_EXPORT macro for exporting classes and functions.

 On Windows systems you need to decorate your functions or classes with __declspec(dllexport) if you wanted to use outside of the dll,
 and you need to use __declspec(dllimport) when you include the header file from other projects.

 This macro helps you to do this.
 \code
 // Exported function:
 GEL_EXPORT void Function1();

 // Exported class:
 class GEL_EXPORT Class1
 {
 ...
 };
 \endcode
*/
#ifdef _WIN32
#  ifdef GEL_EXPORTS
    /* building the lib */
#    define GEL_EXPORT __declspec(dllexport)
#  else
    /* building user code */
//#    define GEL_EXPORT __declspec(dllimport)
#    define GEL_EXPORT
#  endif
#else
#  define GEL_EXPORT
#endif

#ifdef _MSC_VER
#  pragma warning(disable: 4251)
#endif

/*! \namespace gel
 \brief Namespace for Gaze and Emotions Library
*/

typedef double data_t;

// this is a macro which prints a variable name together with its value
// it is useful for debugging
//#define PRINT(x) std::cout << "_debug: " << #x << " = " << x << "\n"

/*!
 \def DEFINE_ENUM_VALUES
 \brief Helper function to create enum values from list of EDEF's
 \code
 // Define an EDEF list:
 #define FOO_ENUM(EDEF)\\
    EDEF(ELEMENT1,1)\\
    EDEF(ELEMENT2,2)

 // Define the enum:
 enum FooEnum { FOO_ENUM(DEFINE_ENUM_VALUES) };
 \endcode
*/
#define DEFINE_ENUM_VALUES(x,y) x = y,

#endif /* CORE_DEFINES_H */
