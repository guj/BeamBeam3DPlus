/*
 *      Author: Junmin Gu
 */

#ifndef SENSEI_BEAM_BINDINGS_C_TYPE_H_
#define SENSEI_BEAM_BINDINGS_C_TYPE_H_


#include <stdint.h>
#include <mpi.h>

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" {
#endif

// From the C side, we use an opaque pointer.
  typedef struct SenseiHandleWriter SenseiHandleWrier;
  typedef struct SenseiHandleReader SenseiHandleReader;

#ifdef __cplusplus
} // end extern C
#endif


#endif /* SENSEI_BEAM_BINDINGS_C_TYPE_H_ */
