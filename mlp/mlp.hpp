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

#include "cast.h"

const int gel_xml_export_version = 1;

class exception : public std::exception {
private:
    std::string m_message;
    
public:
    
    exception( const std::string & message ) throw() : m_message(message) { }

    const char* what() const throw() {return m_message.c_str(); }

    ~exception() throw() {}; 
}; // exception

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
    return;            

} /* save */

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
            throw ::exception( std::string("serialization error: the matrix data type (") + gel::cast<std::string>(matrix.type()) + ") in the file is different from the class type (" + gel::cast<std::string>(type) + ")" );

        matrix.create(rows, cols);
                
        int row;
        for ( int q = 0; q < matrix.rows; q++ )
        {
            archiver >> BOOST_SERIALIZATION_NVP(row);                
                    
            archiver >> boost::serialization::make_array( reinterpret_cast<T1*>(matrix.ptr(q)), matrix.cols );
        }
        break;
                
    default:
        throw ::exception ( std::string("serialization error: file version (") + gel::cast<std::string>(gel_xml_export_version) + ") not recognized " );
    } /* switch */

    return;            
} /* load */

// splitting the serialization into two functions
template<class T0, class T1> inline void boost::serialization::serialize( T0 & archiver, cv::Mat_<T1> & matrix, const unsigned int version )
{
    boost::serialization::split_free( archiver, matrix, version );
}

#define GEL_EXPORT_VAR( archiver, var )                                 \
    try                                                                 \
    {                                                                   \
        archiver & boost::serialization::make_nvp( BOOST_PP_STRINGIZE(var), var); \
    }                                                                   \
    catch (boost::archive::archive_exception exception)                 \
    {                                                                   \
        throw ::exception( std::string("error: importing serialization; at member variable '") + BOOST_PP_STRINGIZE(var) + "' (message: " + exception.what() + ")" ); \
    }


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

    
    public:

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

} // namespace gel  

#endif /* CLM_MLP_H */
