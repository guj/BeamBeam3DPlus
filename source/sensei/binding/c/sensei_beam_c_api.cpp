#include "sensei_beam_c_api.h"
#include "VTKBeamBeam3D.h"

#include <iostream>
using namespace std;

//
// READER
//
SenseiHandleReader* SenseiHandleReader_new(MPI_Comm comm)
{
  SenseiHandleReader*  reader = nullptr;
    try
      {
        reader = reinterpret_cast<SenseiHandleReader*>(new VTKBeamBeam3DReader(comm));
      }
    catch (...)
      {
	std::cerr<<"Unable to create new reader"<<std::endl;
      }
    return reader;
}

SenseiHandleReader* SenseiHandleReader_createFile(MPI_Comm comm)
{
  return SenseiHandleReader_new(comm);
}


int SenseiHandleReader_initXML(SenseiHandleReader* reader, const char* name)
{
  return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->InitializeXML(name);
  //return reader->InitializeXML(name);
}

int SenseiHandleReader_getInput(SenseiHandleReader* reader,  uint64_t* nturns, uint64_t* nptls, uint32_t* bunchID)
{
  if (reader == NULL)
    return -1;

  if (nptls != NULL)
    *nptls =   reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetNumParticles();

  if (nturns != NULL)
    *nturns =  reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetNumTurns();

  if (bunchID != NULL)
    *bunchID = reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetBunchID();

  std::cout<<*nptls<<"  "<<*nturns<<"  "<<*bunchID<<std::endl;
  return 0;
  //return reader->InitializeXML(name);
}


uint64_t SenseiHandleReader_selSize(SenseiHandleReader* reader)
{
  return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetSelectionSize();
  //return reader->GetSelectionSize();
}

//
// get how many particles read per turn
//
uint64_t SenseiHandleReader_nPtls(SenseiHandleReader* reader)
{
  //return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetPtlsReadPerTurn();
  return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetNumParticles();
  //return reader->GetPtlsReadPerTurn();
}

//
// get how many valid turns read in
//
uint64_t SenseiHandleReader_nTurns(SenseiHandleReader* reader)
{
  if (reader == NULL) {
    printf("ERROR, NULL READER\n");
    return 0;
  }

  printf("##DEBUG: Num PARTICLES: %llu\n", reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetNumParticles());
  printf("##DEBUG: Num TURNS: %llu\n", reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetNumTurns());
  //return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetTurnsRead();
  return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->GetNumTurns();
  //return reader->GetTurnsRead();
}

int  SenseiHandleReader_getData(SenseiHandleReader* reader, double* data1, double* data2)
{
  return reinterpret_cast<VTKBeamBeam3DReader*>(reader)->readPtl(data1, data2);;
  //return reader->readPtl(data1, data2);
}


void SenseiHandleReader_remove(SenseiHandleReader* Sensei)
{
  if (Sensei != NULL) 
    delete reinterpret_cast<VTKBeamBeam3DReader *>(Sensei);

}
//
// WRITER
//
SenseiHandleWriter* SenseiHandleWriter_createFile(MPI_Comm comm, const char* name)
{
  SenseiHandleWriter*  writer = nullptr;
  try
    {
      if ((name == NULL) || (strlen(name) == 0)) {
	VTKBeamBeam3DWriter* r = new VTKBeamBeam3DWriter(comm);
	r->InitializeXML("sensei.xml");
	writer = reinterpret_cast<SenseiHandleWriter*>(r);
      } else {
	writer = reinterpret_cast<SenseiHandleWriter*>(new VTKBeamBeam3DWriter(comm, name));
      }
    }
  catch (...)
    {
      std::cerr<<"Unable to create new writer"<<std::endl;
    }
  
  return writer;
}


void SenseiHandleWriter_remove(SenseiHandleWriter* Sensei)
{
  cout << "C API, delete_Sensei " << Sensei<<endl;
  if (Sensei != NULL) 
    delete reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei);
}

int  SenseiHandleWriter_initXML(SenseiHandleWriter* Sensei, const char* name)
{
  return reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei)->InitializeXML(name);
  //return Sensei->InitializeXML(name);
}

int  SenseiHandleWriter_turnStart(SenseiHandleWriter* Sensei, int turnID)
{
  return reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei)->turnStart(turnID);
  //return Sensei->turnStart(turnID);
}

int  SenseiHandleWriter_setupPtl(SenseiHandleWriter* Sensei,  int nBunchs, uint64_t start, uint64_t count, int outRate)
{
  return reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei)->setupPtl(nBunchs, start, count, outRate);
  //return Sensei->setupPtl(nBunchs, start, count, outRate);
}

int  SenseiHandleWriter_writeAttr(SenseiHandleWriter* Sensei, int currBunch, const char* name, double* attrData)
{
  return reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei)->writeAttr(currBunch, name, attrData);
  //return Sensei->writeAttr(currBunch, name, attrData);
}

int  SenseiHandleWriter_writePtl(SenseiHandleWriter* Sensei, int currBunch, double* attrData)
{
  return reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei)->writePtl(currBunch, attrData);
  //return Sensei->writePtl(currBunch,  attrData);
}


int  SenseiHandleWriter_turnEnd(SenseiHandleWriter* Sensei)
{
  return reinterpret_cast<VTKBeamBeam3DWriter*>(Sensei)->turnEnd();
  //return Sensei->turnEnd();
}




