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
#include <boost/preprocessor.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "GEL/histo.h"

const int gel_xml_export_version = 1;

namespace boost
{
    namespace serialization
    {

        /**
        \brief Instantiated explicit export of the cv::Mat_ class.
        This function handles the export of the cv::Mat_ class into the XML library.

        \param archiver The serializator object (created by boost)
        \param matrix The cv::Mat_ object to serialize
        \param int unused
        */
        template<class T0, class T1>
        void save( T0 & archiver, const cv::Mat_<T1> & matrix, unsigned int );

        /**
           \brief Instantiated explicit export of the cv::Mat_ class.
           This function handles the export of the cv::Mat_ class into the XML library.

           \param archiver The serializator object (created by boost)
           \param matrix The cv::Mat_ object to serialize
           \param int unused
        */
        template<class T0, class T1>
        void load( T0 & archiver, cv::Mat_<T1> & matrix, unsigned int );

        // splitting the serialization into two functions
        template<class T0, class T1> inline
        void serialize( T0 & archiver, cv::Mat_<T1> & matrix, const unsigned int version );

        template<class T0>
        void
        serialize( T0 & archiver, struct tm & tm, const unsigned int );

    } /* namespace serialization */
} /* namespace boost */

template<class T0, class T1>
void boost::serialization::save( T0 & archiver, const cv::Mat_<T1> & matrix, unsigned int )
{
    // this file has been prepared for version 1
    archiver << BOOST_SERIALIZATION_NVP(gel_xml_export_version);

    int rows = matrix.rows;
    int cols = matrix.cols;
    int type = matrix.type();

    archiver << BOOST_SERIALIZATION_NVP(rows);
    archiver << BOOST_SERIALIZATION_NVP(cols);
    archiver << BOOST_SERIALIZATION_NVP(type);

    const T1 * tmp;
    for ( int row = 0; row < matrix.rows; row++ )
    {
        archiver << BOOST_SERIALIZATION_NVP(row);

        tmp = reinterpret_cast<const T1*>(matrix.ptr(row));
        archiver << boost::serialization::make_array( const_cast<T1*>(tmp), matrix.cols );
    }
}

template<class T0, class T1>
void boost::serialization::load( T0 & archiver, cv::Mat_<T1> & matrix, unsigned int )
{
    // this file has been prepared for version 1

    int gel_xml_export_version;
    archiver >> BOOST_SERIALIZATION_NVP(gel_xml_export_version);

    // here we list all the known versions
    switch (gel_xml_export_version)
    {
    case 1:
        // this is the version 1 file reading

        int rows;
        int cols;
        int type;

        archiver >> BOOST_SERIALIZATION_NVP(rows);
        archiver >> BOOST_SERIALIZATION_NVP(cols);
        archiver >> BOOST_SERIALIZATION_NVP(type);

        if (matrix.type()!=type)
            throw std::runtime_error( std::string("serialization error: the matrix data type (") + std::to_string(matrix.type()) + ") in the file is different from the class type (" + std::to_string(type) + ")" );

        matrix.create(rows, cols);

        int row;
        for ( int q = 0; q < matrix.rows; q++ )
        {
            archiver >> BOOST_SERIALIZATION_NVP(row);

            archiver >> boost::serialization::make_array( reinterpret_cast<T1*>(matrix.ptr(q)), matrix.cols );
        }
        break;

    default:
        throw std::runtime_error( std::string("serialization error: file version (") + std::to_string(gel_xml_export_version) + ") not recognized " );
    }
}

// splitting the serialization into two functions
template<class T0, class T1> inline void boost::serialization::serialize( T0 & archiver, cv::Mat_<T1> & matrix, const unsigned int version )
{
    boost::serialization::split_free( archiver, matrix, version );
}

#define GEL_EXPORT_VAR( archiver, var ) archiver & boost::serialization::make_nvp( BOOST_PP_STRINGIZE(var), var);

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

        /*!
          \brief Generates response map around a given coordinate

          \param image Image, [n x m], uint8
          \param center Center coordinates of the response map, [2 x 1], integer
          \param mapSize Radius like size of the response map, integer
          \return Response of the classifier, [(2*mapSize+1) x (2*mapSize+1)] double
        */
        cv::Mat_<pixel_type> generateResponseMap(const cv::Mat_<uint8_t>& image, const cv::Point2i& center, int mapSize ) const;

    public:

        int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        cv::Mat_<pixel_type> m_wIn; /*!< \brief */
        cv::Mat_<pixel_type> m_wOut; /*!< \brief  */
        cv::Mat_<pixel_type> m_U; /*!< \brief */
        int hidden_num;
        double rho2;

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

} // namespace gel

/************************************************ implementation of GEL functions ************************************************/

template <typename T0>
cv::Mat_<typename gel::MLP<T0>::pixel_type> gel::MLP<T0>::generateResponseMap( const cv::Mat_<uint8_t>& image,
                             const cv::Point2i& center,
                             int mapSize ) const
{
    // make sure that we have the necessary matrices
    // calculated (this is a method; it is NOT thread safe!), but it is called for each classifier once
    update(mapSize);

    for ( int ncy = 0, cy = center.y - mapSize; cy < center.y + mapSize + 1; ++ncy, ++cy ) {
        cv::Mat_<uint8_t> sample = image( cv::Range( cy - m_patchSize, cy + m_patchSize + 1 ), cv::Range( center.x - mapSize - m_patchSize, center.x + mapSize + m_patchSize + 2 ) );
        cv::Mat_<pixel_type> pack_patches = generatePatch(sample, m_patch_line_size);

        cv::Mat_<pixel_type> target = m_all_patches( cv::Range::all(), cv::Range( ncy * m_number_of_patches_per_line, (ncy + 1) * m_number_of_patches_per_line ) );
        pack_patches.copyTo(target);
    }

    cv::Mat_<pixel_type> result = evaluateSamples();

    return result.cv::Mat::reshape(0,m_number_of_patches_per_line);
}

template <typename T0>
void gel::MLP<T0>::update(int newMapSize) const
{
    // this function is interesting only if the mapSize has
    // changed
    if ( newMapSize != m_currentMapSize )
    {
        cv::Mat_<pixel_type> target;

        m_currentMapSize = newMapSize;

        // extract to the creation
        m_wIn_gemm = m_wIn.colRange(0,m_wIn.cols-1)*m_U.t();
        cv::Mat_<pixel_type> bIn = m_wIn.col(m_wIn.cols-1);
        m_wOut_gemm = m_wOut.colRange(0,m_wOut.cols-1);

        pixel_type bOut = m_wOut(0,m_wOut.cols-1);

        m_number_of_patches_per_line = 2 * m_currentMapSize + 1;
        // here we assume that the patches are always square-like
        int number_of_all_patches  = m_number_of_patches_per_line * m_number_of_patches_per_line;
        m_patch_line_size = 2 * m_patchSize + 1;

        m_all_bIns.create ( bIn.rows, number_of_all_patches );
        m_all_bOuts.create( 1, number_of_all_patches );
        m_all_patches.create( m_wIn_gemm.cols, number_of_all_patches );

        for (int q=0; q < number_of_all_patches; q++)
        {
            target = m_all_bIns( cv::Range::all(), cv::Range(q, q+1) );
            bIn.copyTo(target);
            m_all_bOuts( 0, q ) = bOut;
        }
    }
}

template <typename T0>
cv::Mat_<typename gel::MLP<T0>::pixel_type> gel::MLP<T0>::generatePatch(const cv::Mat_<uint8_t>& sample,int patch_size) const
{
    cv::Mat_<pixel_type> result( sample.rows * patch_size, sample.cols - patch_size );

    // the first patch is treated specially
    cv::Mat_<uint8_t> patch_image = sample( cv::Range::all(), cv::Range(0, patch_size) );

    histo<pixel_type> stat;
    stat.insert( patch_image );

    // we are now processing line by line
    for ( int q = patch_size, col = 0;
        q < sample.cols;
        q++, col++ )
    {
        pixel_type sample_mean = stat.mean();
        pixel_type sample_max = stat.max() - sample_mean;
        pixel_type sample_min = stat.min() - sample_mean;

        sample_max = std::max( std::abs( sample_max ), std::abs( sample_min ) );
        if (sample_max==0) sample_max = 1;

        patch_image = sample( cv::Range::all(), cv::Range( col, col + patch_size ) );

        cv::Mat_<pixel_type> normalized_patch;
        patch_image.convertTo (
            normalized_patch, normalized_patch.type(),
            1.0/sample_max, -(1.0/sample_max)*sample_mean );

        cv::Mat_<pixel_type> patch_line = normalized_patch.cv::Mat::reshape(0,normalized_patch.rows * normalized_patch.cols);
        cv::Mat_<pixel_type> target = result( cv::Range::all(), cv::Range( col, col + 1) );

        patch_line.copyTo( target );

        cv::Mat_<uint8_t> patch_line_in  = sample( cv::Range::all(), cv::Range( q, q + 1 ) );
        cv::Mat_<uint8_t> patch_line_out = sample( cv::Range::all(), cv::Range( col, col + 1 ) );
        stat.insert( patch_line_in );
        stat.remove( patch_line_out );
    }

    return result;
}

template <typename T0>
cv::Mat_<typename gel::MLP<T0>::pixel_type> gel::MLP<T0>::evaluateSamples() const
{
    cv::Mat_<pixel_type> xOut;
    cv::gemm( m_wIn_gemm, m_all_patches, -1., m_all_bIns, -1., xOut );
    cv::Mat_<pixel_type> e;
    cv::exp( xOut, e );
    cv::add( e, 1.0, xOut );
    cv::divide( 2.0, xOut, e );
    cv::subtract( e, 1.0, xOut );

    cv::Mat_<pixel_type> dot;
    cv::gemm( m_wOut_gemm, xOut, -1., m_all_bOuts, -1., dot );

    cv::Mat_<pixel_type> result;
    cv::exp( dot, result );
    cv::add( result, 1.0, dot );
    cv::divide( 1.0, dot, result );

    return result;
}



#endif /* CLM_MLP_H */
