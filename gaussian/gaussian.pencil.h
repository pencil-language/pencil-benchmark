#ifdef __cplusplus
extern "C" {
#endif
void pencil_gaussian( const int rows
                    , const int cols
                    , const int step
                    , const float src[]
                    , const int kernelX_rows
                    , const int kernelX_cols
                    , const int kernelX_step
                    , const float kernelX[]
                    , const int kernelY_rows
                    , const int kernelY_cols
                    , const int kernelY_step
                    , const float kernelY[]
                    , float conv[]
                    );
#ifdef __cplusplus
} // extern "C"
#endif
