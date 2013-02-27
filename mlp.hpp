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

//  UjoImro, 2012

#ifndef CLM_MLP_H
#define CLM_MLP_H

#include <vector>
#include <fstream>
#include <iostream>
#include <stdint.h>

#include <opencv2/core/core.hpp>

#include "Core/defines.h"

namespace gel
{

/*!
\brief Multi-Layer Perceptron based classifier class

*/
template <typename T0>
class MLP
{
public:
    typedef T0 pixel_type;

    template <class MT0>
    void
    serialize( MT0 & archiver, unsigned int ) {
        GEL_EXPORT_VAR( archiver, m_patchSize );
        GEL_EXPORT_VAR( archiver, m_wIn );
        GEL_EXPORT_VAR( archiver, m_wOut );
        GEL_EXPORT_VAR( archiver, m_U );
        GEL_EXPORT_VAR( archiver, hidden_num );
        GEL_EXPORT_VAR( archiver, rho2 );
    } // serialize

    
private:

    int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
    cv::Mat_<pixel_type> m_wIn; /*!< \brief */
    cv::Mat_<pixel_type> m_wOut; /*!< \brief  */
    cv::Mat_<pixel_type> m_U; /*!< \brief */
    int hidden_num;
    double rho2;

    enum NormalizationMethod { none, maxAbs, meanStd };

private:
    NormalizationMethod preSVDNormalizationMethod;
    NormalizationMethod postSVDNormalizationMethod;
    
    // we store the data for the response map generation if the
    // mapSize parameter does not change (and there is no reason
    // why it should in the same execution), than we can reuse the
    // matrices for the response map generation
    mutable int m_currentMapSize;
    mutable int m_number_of_patches_per_line;
    mutable int m_patch_line_size;
    mutable cv::Mat_<pixel_type> m_all_patches; /*!< \brief This matrix either holds the normalized patches, or the training samples (aka m_U.t()*patches - if post method = none)*/
    mutable cv::Mat_<pixel_type> m_all_bIns;
    mutable cv::Mat_<pixel_type> m_all_bOuts;
    mutable cv::Mat_<pixel_type> m_wIn_gemm; /*!< \brief This matrix either holds m_wIn * m_U.t() (if post method = none) or m_wIn*/
    mutable cv::Mat_<pixel_type> m_wOut_gemm;
};

}  

#endif /* CLM_MLP_H */
