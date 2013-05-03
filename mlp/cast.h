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

// UjoImro, 2012 

// this file performs semantic casting

#ifndef CAST_H_
#define CAST_H_

#include <sstream>

namespace gel
{

    template <class RT, class T0>
    RT cast( const T0 & from )
    {
        std::stringstream stringstream;
        RT to;    
        stringstream << from;
        stringstream >> to;
        return to;    
    } // gel_cast
     
        
} // namespace gel


#endif /* CAST_H_ */

// LuM end of file
