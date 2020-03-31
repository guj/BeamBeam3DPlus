/*
 */

#include "sensei_f2c_beam.h"

// cori needs  it
#include <cstring>

#ifdef __cplusplus
extern "C" {
#endif



// this function is not exposed in the public APIs
void FC_GLOBAL(senseihandlewriter_createfile_f2c,
               SENSEIHANDLEWRITER_CREATEFILE_F2C)(SenseiHandleWriter ** sensei, 
						  MPI_Fint *comm,
						  const char *name,
						  int* ierr)
{
  *sensei = SenseiHandleWriter_createFile(MPI_Comm_f2c(*comm), name);

  *ierr =  (*sensei == NULL)? -1 : 0;
    
}


void FC_GLOBAL(senseihandlewriter_remove_f2c,
	       SENSEIHANDLEWRITER_REMOVE_F2C)(SenseiHandleWriter ** sensei,
					      int* ierr)
{
  if (*sensei != NULL) 
    SenseiHandleWriter_remove(*sensei);

  *ierr = 0;
}

void FC_GLOBAL(senseihandlewriter_setupptl_f2c,
	       SENSEIHANDLEWRITER_SETUPPTL_F2C)(SenseiHandleWriter ** sensei,
						int* nBunchs, uint64_t* start, uint64_t* count, int* outRate,
						int* ierr)
{
  if (*sensei == NULL) {
    *ierr = -1;
    return;
  }

  *ierr = SenseiHandleWriter_setupPtl(*sensei, *nBunchs, *start, *count, *outRate); 
}

void FC_GLOBAL(senseihandlewriter_turnstart_f2c,
	       SENSEIHANDLEWRITER_TURNSTART_F2C)(SenseiHandleWriter ** sensei,
						 int* turnID,
						 int* ierr)
{
  if (*sensei == NULL) {
    *ierr = -1;
    return;
  }
  *ierr = SenseiHandleWriter_turnStart(*sensei, *turnID);
}


void FC_GLOBAL(senseihandlewriter_turnend_f2c,
	       SENSEIHANDLEWRITER_TURNEND_F2C)(SenseiHandleWriter ** sensei, 
					       int* ierr)
{
  if (*sensei == NULL) {
    *ierr = -1;
    return;
  }
  *ierr = SenseiHandleWriter_turnEnd(*sensei);
}


void FC_GLOBAL(senseihandlewriter_attr_f2c,
	       SENSEIHANDLEWRITER_ATTR_F2C)(SenseiHandleWriter ** sensei, 
					    int* bunchID,
					    const char *attrName,
					    double* data,
					    int* ierr)
{
  if (*sensei == NULL) {
    *ierr = -1;
    return;
  }
  *ierr = SenseiHandleWriter_writeAttr(*sensei, *bunchID, attrName, data);
}


  //
  //
  //

void FC_GLOBAL(senseihandlereader_createfile_f2c,
               SENSEIHANDLEREADER_CREATEFILE_F2C)(SenseiHandleReader ** sensei, 
						  MPI_Fint *comm,  
						  const char *name,
						  int* ierr)
{
  *sensei = SenseiHandleReader_createFile(MPI_Comm_f2c(*comm));
  *ierr =  (*sensei == NULL)? -1 : 0;
  if (*sensei != NULL) {
    if ((name != NULL)  &&  (strlen(name) > 0))
      SenseiHandleReader_initXML(*sensei, name);  
    else
      SenseiHandleReader_initXML(*sensei, "sensei.xml");
  }
}

  /*
void FC_GLOBAL(senseihandlewriter_initxml_f2c,
               SENSEIHANDLEWRITER_INITXML_F2C)(SenseiHandleReader ** sensei, 
					       const char *name,
					       int* ierr)
{
  *ierr = SenseiHandleReader_initXML(*sensei, name);  
}
  */

void FC_GLOBAL(senseihandlereader_get_input_f2c,
               SENSEIHANDLEREADER_GET_INPUT_F2C)(SenseiHandleReader ** sensei, 
						 uint64_t* nturns,
						 uint64_t* nptls,
						 uint32_t* bunchID,
						 int* ierr)
{
  SenseiHandleReader_getInput(*sensei, nturns, nptls, bunchID);  
}


void FC_GLOBAL(senseihandlereader_remove_f2c,
	       SENSEIHANDLEREADER_REMOVE_F2C)(SenseiHandleReader ** sensei,
					      int* ierr)
{
  if (*sensei != NULL) 
    SenseiHandleReader_remove(*sensei);

  *ierr = 0;
}

void FC_GLOBAL(senseihandlereader_get_num_turns_f2c,
	       SENSEIHANDLEREADER_GET_NUM_TURNS_F2C)(SenseiHandleReader ** sensei,
						     uint64_t* result,
						     int* ierr)
{
  if (*sensei == NULL) {
    printf("FC Error: NULL SENSEI PTR");
    *ierr = -1;
    return;
  }
  *result = SenseiHandleReader_nTurns(*sensei);
  *ierr = 0;
}


void FC_GLOBAL(senseihandlereader_get_num_particles_f2c,
	       SENSEIHANDLEREADER_GET_NUM_PARTICLES_F2C)(SenseiHandleReader ** sensei,
							 uint64_t* result,
							 int* ierr)
{
  if (*sensei == NULL) {
    *ierr = -1;
    return;
  }
  *result = SenseiHandleReader_nPtls(*sensei);
  *ierr = 0;
}


void FC_GLOBAL(senseihandlereader_get_data_f2c,
	       SENSEIHANDLEREADER_GET_data_F2C)(SenseiHandleReader ** sensei,
						double* data1,
						double* data2,
						int* ierr)
{
  if (*sensei == NULL) {
    *ierr = -1;
    return;
  }
  
  *ierr = SenseiHandleReader_getData(*sensei, data1, data2);
}


#ifdef __cplusplus
}
#endif
