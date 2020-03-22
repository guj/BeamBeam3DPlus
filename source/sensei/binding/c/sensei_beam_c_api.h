/*
 *      Author: Junmin Gu
 */

#ifndef SENSEI_BEAM_BINDINGS_C_H_
#define SENSEI_BEAM_BINDINGS_C_H_

#include "sensei_beam_c_type.h"

#include <stdint.h>
#include <mpi.h>

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" {
#endif


// Constructor
//VTKSenseiHandle* SenseiHandle_create(int a, int b);
SenseiHandleWriter* SenseiHandleWriter_createFile(MPI_Comm comm, const char* name);

// Destructor
void SenseiHandleWriter_remove(SenseiHandleWriter* foo);

// The const qualificators maps from the member function to pointers to the
// class instances.
int  SenseiHandleWriter_initXML(SenseiHandleWriter* foo, const char* name); 				 

int  SenseiHandleWriter_turnStart(SenseiHandleWriter* foo, int turnID);
int  SenseiHandleWriter_setupPtl(SenseiHandleWriter* foo,  int nBunchs, uint64_t start, uint64_t count, int outRate); 
int  SenseiHandleWriter_writeAttr(SenseiHandleWriter* foo, int currBunch, const char* name, double* attrData);

int  SenseiHandleWriter_writePtl(SenseiHandleWriter* foo, int currBunch, double* attrData);

int  SenseiHandleWriter_turnEnd(SenseiHandleWriter* foo);


  //
  // Constructor reader
  //
SenseiHandleReader* SenseiHandleReader_createFile(MPI_Comm comm);
  //
  // Destructor reader
void SenseiHandleReader_remove(SenseiHandleReader* foo);

int  SenseiHandleReader_initXML(SenseiHandleReader* foo, const char* name);
int  SenseiHandleReader_getInput(SenseiHandleReader* foo, uint64_t* nturn, uint64_t* nptls, uint32_t* bunchID);
int  SenseiHandleReader_getData(SenseiHandleReader* foo, double* data1, double* data2); 
uint64_t  SenseiHandleReader_selSize(SenseiHandleReader* foo); 
uint64_t  SenseiHandleReader_nPtls(SenseiHandleReader* foo); 
uint64_t  SenseiHandleReader_nTurns(SenseiHandleReader* foo); 
//void foo_speaker(const char* s);

#ifdef __cplusplus
}
#endif

#endif /* SENSEI_BEAM_BINDINGS_C_H_ */
