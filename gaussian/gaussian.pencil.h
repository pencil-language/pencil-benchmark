#ifdef __cplusplus
extern "C" {
#endif
void pencil_gaussian( const int rows
                    , const int cols
                    , const int step
                    , const float src[]
                    , const int kernelX_length
                    , const float kernelX[]
                    , const int kernelY_length
                    , const float kernelY[]
                    , float conv[]
                    );
#ifdef __cplusplus
} // extern "C"
#endif
