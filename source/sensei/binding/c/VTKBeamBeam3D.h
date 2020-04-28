#include <string>
#include <mpi.h>

#include "ConfigurableInTransitDataAdaptor.h"
#include "ConfigurableAnalysis.h"

#include "VTKDataAdaptor.h"
//#include "HDF5DataAdaptor.h"
//#include "HDF5AnalysisAdaptor.h"
#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkCompositeDataIterator.h>
#include "XMLUtils.h"

//using H5AnalysisAdaptorPtr = vtkSmartPointer<sensei::HDF5AnalysisAdaptor>;
//using H5DataAdaptorPtr = vtkSmartPointer<sensei::HDF5DataAdaptor>;
using SenseiAnalysisAdaptorPtr = vtkSmartPointer<sensei::ConfigurableAnalysis>;
using SenseiDataAdaptorPtr = vtkSmartPointer<sensei::ConfigurableInTransitDataAdaptor>;

class VTKBeamBeam3DReader 
{
public:
  VTKBeamBeam3DReader(MPI_Comm comm);
  VTKBeamBeam3DReader(MPI_Comm comm, const char* name);
  ~VTKBeamBeam3DReader();

  bool InitializeXML(const char* name);
  //bool setAttr(uint32_t pos1, uint32_t pos2);
  bool setAttr(const char* n1, const char* n2);
  bool setBunchID(uint32_t bid);
  bool setTurnRange(uint32_t tStart, uint32_t tSize);
  bool setParticleRange(uint64_t pStart, uint64_t pSize);
  bool readPtl(double* data1, double* data2);
  bool readPtlSerial(); 

  void showError(const std::string& msg);

  uint32_t GetBunchID() {return m_SelectionBunchID;}
  uint64_t GetSelectionSize() {return m_SelectionParticleCount * m_SelectionTurnCount;}

  uint64_t GetNumTurns() {return m_SelectionTurnCount;}
  uint64_t GetNumParticles() {return m_SelectionParticleCount;}
  //uint64_t GetNumTurns() {return m_TurnsRead;}
  //uint64_t GetNumParticles() {return m_NumPtlPerTurn;}

  uint64_t GetTurnsRead() {return m_TurnsRead;}
  uint64_t GetPtlsReadPerTurn() {return m_NumPtlPerTurn;}

  bool m_Debug = false;
  bool m_EndOfStream = false;
  uint64_t m_TurnsRead = 0;
  uint64_t m_NumPtlPerTurn = 0;
 private:

  //H5DataAdaptorPtr m_Reader;
  SenseiDataAdaptorPtr m_Reader;
  uint64_t readMultiBlock(vtkMultiBlockDataSet* ds, sensei::MeshMetadataPtr mmd);
  bool isValidTurn(uint32_t t, uint32_t base);
  MPI_Comm m_Comm;

  uint32_t m_SelectionBunchID=0; 

  uint64_t m_SelectionParticleStart=0;
  uint64_t m_SelectionParticleCount=0;

  uint32_t m_SelectionTurnStart=0;
  uint32_t m_SelectionTurnCount=0;
  uint32_t m_TurnRepeat = 0;
  uint32_t m_TurnRepeatCounter = 0; 

  //uint32_t m_SelectionPos1=0;
  //uint32_t m_SelectionPos2=0;

  std::string m_SelectionColName1; 
  std::string m_SelectionColName2;

  double* m_data1 = NULL;
  double* m_data2 = NULL;  

  int       m_Rank;
  int       m_Size;
  
};



class VTKBeamBeam3DWriter 
{
public:
  VTKBeamBeam3DWriter(MPI_Comm comm);
  VTKBeamBeam3DWriter(MPI_Comm comm, const char* name);
  ~VTKBeamBeam3DWriter();

  bool InitializeXML(const char* name);
  int setupPtl(int nBunchs, uint64_t start, uint64_t count, int outRate); 
  /*
  int setupPhaseout(int nBunchs, uint64_t total, uint64_t start, uint64_t count);
  int phaseout(double* data, int currbunch, int itic, int64_t iseedinit, double* close2g);
  */
  int turnStart(int turnID);
  int turnEnd();
  
  int writeAttr(int currBunch, const std::string& attrName, double* attrData);
  void printMe();

  int writePtl(int currBunch, double* data); // first attemp. not efficient for attr based access pattern
private:
  MPI_Comm m_Comm;
  
  //H5AnalysisAdaptorPtr m_Writer;
  SenseiAnalysisAdaptorPtr m_Writer;

  std::string m_MeshName;
  vtkSmartPointer<vtkMultiBlockDataSet> m_Mesh;
  vtkSmartPointer<sensei::VTKDataAdaptor> m_DataAtTurn;

  uint64_t  m_Start;
  uint64_t  m_Count;
  int       m_BunchTotal;

  int       m_Rank;
  int       m_Size;

  //int writePtl(int currBunch, double* data); // first attemp. not efficient for attr based access pattern
};


//void foo_speaker(std::string s);

//integer*8, dimension(3) :: vTotal, vStart, vCount    ! for var_handle
//integer*8, dimension(3) :: phTotal, phStart, phCount ! for phvar_handle
