/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012-2014
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

#ifndef NEL_FAST_APPROXIMATE_MATH_H
#define NEL_FAST_APPROXIMATE_MATH_H

#include <cmath>
#include <type_traits>

namespace nel {
    // Approximates atan2(y, x) normalized to the [0,4) range
    // with a maximum error of 0.1620 degrees
    inline float normalized_atan2(float y, float x) {
        static const float b = 0.596227f;

        if (x == 0.0f && y == 0.0f)
            return 0.0f;

        // Extract the sign bits
        bool ux_s = x < 0.0f;
        bool uy_s = y < 0.0f;

        // Determine the quadrant offset
        float q = ux_s ? 2.0f : (uy_s ? 4.0f : 0.0f);

        // Calculate the arctangent in the first quadrant
        float bxy_a = std::abs(b * x * y);
        float num = bxy_a + y * y;
        float atan_1q = num / (x * x + bxy_a + num);

        // Translate it to the proper quadrant
        float uatan_2q = (ux_s != uy_s) ? -atan_1q : atan_1q;
        return q + uatan_2q;
    }

    // Approximates atan(x) normalized to the [-1,1] range
    // with a maximum error of 0.1620 degrees.
    inline float normalized_atan( float x )
    {
        static const float b = 0.596227f;

        // Extract the sign bit
        bool ux_s = x < 0.0f;

        // Calculate the arctangent in the first quadrant
        float bx_a = std::abs(b * x);
        float num = bx_a + x * x;
        float atan_1q = num / (1.0f + bx_a + num);

        // Restore the sign bit
        float atan_2q = ux_s ? -atan_1q : atan_1q;
        return atan_2q;
    }

    template<typename T>
    inline int fast_floor(T f) {
        static_assert(std::is_floating_point<T>::value, "fast_floor: Parameter must be floating point.");
        return static_cast<int>(f)-(f < T(0));
    }

    template<typename T>
    inline int fast_ceil(T f) {
        static_assert(std::is_floating_point<T>::value, "fast_floor: Parameter must be floating point.");
        return static_cast<int>(f)+(f > T(0));
    }
}

#endif /* NEL_FAST_APPROXIMATE_MATH_H */
