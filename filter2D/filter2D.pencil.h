#ifdef __cplusplus
extern "C" {
#endif
    void pencil_filter2D( int rows
                        , int cols
                        , int step
                        , float src[]
                        , int kernel_rows
                        , int kernel_cols
                        , int kernel_step
                        , float kernel[]
                        , float conv[]
                        );
#ifdef __cplusplus
} // extern "C"
#endif
