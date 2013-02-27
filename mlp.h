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

    /*!
      \brief Constructor

    */
    MLP();

    /*!
      \brief Copy Constructor

      Creates a deep copy of the object.

      \param other
    */
    MLP(const MLP & other);

    /*!  \brief Gets radius like patch size, the true size of the
      patch is [(2*patchSize+1) x (2*patchSize+1)]
    */
    int patchSize() const;

    /*!
      \brief Generates response map around a given coordinate

      \param image Image, [n x m], uint8
      \param center Center coordinates of the response map, [2 x 1], integer
      \param mapSize Radius like size of the response map, integer
      \return Response of the classifier, [(2*mapSize+1) x (2*mapSize+1)] double
    */
    cv::Mat_<pixel_type> generateResponseMap(const cv::Mat_<uint8_t>& image, const cv::Point2i& center, int mapSize ) const;

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

    /*!
       This function updates the mutable variables which are used
       for the response map generation. This function has to be
       called when a new potential mapSize comes in the function,
       but it recalculates the matrices only if the mapSize of the
       requested response map changes.

       IMPORTANT: This function is NOT thread safe. The
       responseMaps are generated separately for each control
       point. This function however can be made thread safe for
       the cost of the synchronization overhead.

       \param new_mapSize
     */
    void update(int newMapSize ) const;

    /**
     * Generates a matrix with normalized patches.
     *
     * @param sample
     * @param patch_size
     *
     * @return
     */
    cv::Mat_<pixel_type> generatePatch(const cv::Mat_<uint8_t>& sample,int patch_size) const;

    /*!
      \brief Evaluate normalized samples
      
      \return The values
    */
    cv::Mat_<pixel_type> evaluateSamples() const;

    int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
    cv::Mat_<pixel_type> m_wIn; /*!< \brief */
    cv::Mat_<pixel_type> m_wOut; /*!< \brief  */
    cv::Mat_<pixel_type> m_U; /*!< \brief */
    int hidden_num;
    double rho2;

public:    
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

#include "mlp_impl.h"

#endif /* CLM_MLP_H */
