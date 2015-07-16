#ifdef __cplusplus
extern "C" {
#endif
    void pencil_filter2D( const int        rows, const int        cols, const int        step, const float src[]
                        , const int kernel_rows, const int kernel_cols, const int kernel_step, const float kernel[]
                        , float conv[]
                        );
#ifdef __cplusplus
} // extern "C"
#endif
