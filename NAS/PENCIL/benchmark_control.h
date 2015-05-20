#ifndef __BENCHMARK_CONTROL_H
#define __BENCHMARK_CONTROL_H

/* Used to control whether the original OpenMP directives should be kept.  */
#define USE_ORIGINAL_OMP_DIRECTIVES 	1

/* Without inlining, one can only optimize innermost SCOPs.  After inlining one
 * can use outermost SCOPs.  This variable can be used to control which SCOPs
 * should be considered.  */
#define USE_INNERMOST_SCOPS		1

#endif /* __BENCHMARK_CONTROL_H */
