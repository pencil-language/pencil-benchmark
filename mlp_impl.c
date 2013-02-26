// UjoImro, 2013
// Experimental code of the CARP Project
// Copyright (c) RealEyes, 2013

#include <stdint.h>

const int MAX_INT = 1 << sizeof(int) - 1;

typedef struct {
    int rows;
    int cols;
    int step;
    float * data;
} MatFloat; // struct MatFloat

typedef struct {
    int rows;
    int cols;
    int step;
    uint8_t * data;    
} MatChar; // struct MatChar

typedef struct {
    int x;
    int y;
} Point2i;

MatChar GetBlockChar( MatChar self, int row_from, int row_to, int col_from, int col_to )
{
    
}

MatFloat GetBlockFloat( MatFloat self, int row_from, int row_to, int col_from, int col_to )
{
    
}

MatFloat
generatePatch( MatChar sample, int m_patch_line_size )
{

}

MatFloat
reshape( MatFloat self, int ignored, int new_col_size )
{

}

MatFloat evaluateSamples()
{

}

void update( int mapSize )
{

}

void copyToFloat( MatFloat src, MatFloat dst )
{

}

MatFloat
generateResponseMap( 
    const MatChar image,
    const Point2i center,
    int mapSize, 
    int m_patchSize,
    int m_patch_line_size,
    MatFloat m_all_patches,
    int m_number_of_patches_per_line
    )
{
    // make sure that we have the necessary matrices
    // calculated (this is a method; it is NOT thread safe!), but it is called for each classifier once
    update(mapSize);

    int ncy, cy;
    for ( ncy = 0, cy = center.y - mapSize; cy < center.y + mapSize + 1; ++ncy, ++cy ) {
	MatChar  sample = GetBlockChar( image, cy - m_patchSize, cy + m_patchSize + 1 , center.x - mapSize - m_patchSize, center.x + mapSize + m_patchSize + 2 );
        MatFloat pack_patches = generatePatch( sample, m_patch_line_size );

        MatFloat target = GetBlockFloat( m_all_patches, -MAX_INT, MAX_INT, ncy * m_number_of_patches_per_line, (ncy + 1) * m_number_of_patches_per_line );
        copyToFloat( pack_patches, target );
    }

    MatFloat result = evaluateSamples();

    return reshape( result, 0, m_number_of_patches_per_line );

} // generateResponseMap


// LuM end of file
