#include "VTKBeamBeam3D.h"
#include <iostream>
#include <sstream> 
#include <vtkIndent.h>


using namespace std;

/*

class ParticlesStrucured
{
  vtkImageData*  SaveBlock(unsigned long bunchID, unsigned long start, unsigned long count)
  {
    vtkImageData *im = vtkImageData::New();
    im->SetExtent(bunchID, bunchID, ptlStart, ptlStart + ptlCount-1, 0, 5);
    unsigned long nx = 1;
    unsigned long ny = ptlCount1;
    unsigned long nz = 6;
    //unsigned long npts = (nx + 1) * (ny + 1) * (nz + 1);
    //get_data_arrays(npts, im->GetPointData());
    
    unsigned long ncells = nx * ny * nz;
    
    get_data_arrays(ncells, im->GetCellData());
    
    return im;
  }
  
};
*/

VTKBeamBeam3DWriter::VTKBeamBeam3DWriter(MPI_Comm comm)
  :m_Comm(comm),
   m_MeshName("particles")
{
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &m_Rank);
  MPI_Comm_size(comm, &m_Size);

  m_Writer = sensei::HDF5AnalysisAdaptor::New();
}


bool VTKBeamBeam3DWriter::InitializeXML(const char* xmlname)
{
  pugi::xml_document doc;
  if (sensei::XMLUtils::Parse(m_Comm, xmlname, doc))
    {
      SENSEI_ERROR("Failed to load, parse, and share XML configuration")
	MPI_Abort(m_Comm, -1);
      return false;
    }

  pugi::xml_node node = doc.child("sensei").child("writer");

  pugi::xml_attribute filename = node.attribute("filename");
  pugi::xml_attribute methodAttr = node.attribute("method");

  if (filename)
    m_Writer->SetStreamName(filename.value());
  else {
    SENSEI_ERROR("Failed to proceed without attribute: filename")
    return false;
  }

  if (methodAttr)
    {
      std::string method = methodAttr.value();

      if (method.size() > 0)
        {
          bool doStreaming = ('s' == method[0]);
          bool doCollectiveTxf = ((method.size() > 1) && ('c' == method[1]));

          m_Writer->SetStreaming(doStreaming);
          m_Writer->SetCollective(doCollectiveTxf);
        }
    }

  return true;
}


VTKBeamBeam3DWriter::VTKBeamBeam3DWriter(MPI_Comm comm, const char* name)
  :m_Comm(comm), 
   m_MeshName("particles")
{  
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  //MPI_Init(NULL, NULL);
  MPI_Comm_rank(comm, &m_Rank);
  MPI_Comm_size(comm, &m_Size);

  if (m_Rank == 0) 
    std::cout<<"name = "<<name<<std::endl;


  bool doStreaming = false;
  bool doCollective = true;
  //m_FID = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT,  H5P_DEFAULT); 
  m_Writer = sensei::HDF5AnalysisAdaptor::New();
  m_Writer->SetStreamName(name);
  m_Writer->SetStreaming(doStreaming);
  m_Writer->SetCollective(doCollective);

  m_Writer->SetCommunicator(comm);
}

//VTKBeamBeam3DWriter::VTKBeamBeam3DWriter(int _a, int _b): a(_a), b(_b){
//cout << "C++ side, constructor" << endl;
//}

VTKBeamBeam3DWriter::~VTKBeamBeam3DWriter()
{
  m_Writer->Finalize();
  m_Writer->Delete();
}

/*
*/

int VTKBeamBeam3DWriter::setupPtl(int nBunchs, uint64_t start, uint64_t count, int outRate)
{
  /*
  vStart(1) = 0; vStart(2) = start;  vStart(3) = 0;
  vCount(1) = 6; vCount(2) = count;  vCount(3) = 1;
  vTotal(1) = 6; vTotal(2) = global; vTotal(3) = nbunchs;
  pNbunch = nbunchs

    ! not constant dim. (fixed shape/start/count)                                                                                                                                                     
    call adios2_define_variable(var_handle, a2_io, "particles", adios2_type_dp, 3, &
				vTotal, vStart, vCount, .false., ierr)
  */
  m_Start = start;
  m_Count = count;
  m_BunchTotal = nBunchs;

  m_Mesh = vtkSmartPointer<vtkMultiBlockDataSet>::New(); 
  m_Mesh->SetNumberOfBlocks(m_Size * m_BunchTotal);

  //if (m_Rank == 0) std::cout<<" num blocks = " <<m_Size * m_BunchTotal<<std::endl;
  return true;
}

int VTKBeamBeam3DWriter::writeAttr(int bunchID, const std::string& attrName, double* attrData)
{
  int bid = m_Rank * m_BunchTotal + bunchID;
  
  //if (m_Rank == 1) std::cout<<"bid="<<bid<<" bunch_id="<<bunchID<<"   "<<attrName<<"\t "<<attrData<<" \t"<<attrData[0]<<std::endl;

  unsigned long nx = 1;
  unsigned long ny = m_Count;
  unsigned long nz = 1;


  unsigned long nCells = nx * ny * nz;   
  vtkDoubleArray* temp = vtkDoubleArray::New();
  temp->SetName(attrName.c_str());
#ifdef NEVER // zero copy
  temp->SetArray(attrData, nCells, 1);
#else
  temp->SetNumberOfComponents(1);
  temp->SetNumberOfTuples(nCells);
  
  for (int i = 0; i < nCells; i++)
    {
      temp->SetTuple1(i, attrData[i]);
    }  
#endif
  vtkImageData *im = vtkImageData::SafeDownCast(m_Mesh->GetBlock(bid));
  bool blockexists = (im != nullptr);

  //if (m_Mesh->GetBlock(bid) != NULL) {
  if (im == nullptr) {        
    im = vtkImageData::New();
    im->SetExtent(bunchID, bunchID, m_Start, m_Start + m_Count, 0, 1);
  }
  
  im->GetCellData()->AddArray(temp);
  m_Mesh->SetBlock(bid, im);
    
  if (!blockexists) im->Delete();

  temp->Delete();
  //std::cout<<attrName<<"   block id="<<bid<<" data[0] = "<<attrData[0]<<std::endl;
  return true;
}


int VTKBeamBeam3DWriter::writePtl(int bunchID, double* data)
{
  int bid = m_Rank * m_BunchTotal + bunchID;
  
  vtkImageData *im = vtkImageData::New();
  im->SetExtent(bunchID, bunchID, m_Start, m_Start + m_Count, 0, 6);
  //std::cout<<m_Rank<<" bid="<<bid<<" bunchid="<<bunchID<<"   "<<m_Start<<"  "<<m_Start+m_Count<<std::endl;
  unsigned long nx = 1;
  unsigned long ny = m_Count;
  unsigned long nz = 6;
  //unsigned long npts = (nx + 1) * (ny + 1) * (nz + 1);
  //get_data_arrays(npts, im->GetPointData());
  
  unsigned long nCells = nx * ny * nz;
  
  vtkDoubleArray* temp = vtkDoubleArray::New();
  temp->SetName("ptldata");
#ifdef NEVER
  temp->SetArray(data, nCells, 1); // warning: zero copy...
#else
  temp->SetNumberOfComponents(1);
  temp->SetNumberOfTuples(nCells);
  
  for (int i = 0; i < nCells; i++)
    {
      temp->SetTuple1(i, data[i]);
    }  
#endif
  im->GetCellData()->AddArray(temp);
  m_Mesh->SetBlock(bid, im);

  // std::cout<<"   block id="<<bid<<" data[0] = "<<data[0]<<std::endl;
  temp->Delete();
  im->Delete();

  return true;
}



int VTKBeamBeam3DWriter::turnStart(int turnID)
{
    m_DataAtTurn = sensei::VTKDataAdaptor::New();
    m_DataAtTurn->SetDataTime(turnID);
    m_DataAtTurn->SetDataTimeStep(turnID);    

    return true;
}

void VTKBeamBeam3DWriter::printMe() {
  /*
  vtkCompositeDataIterator* it = m_Mesh->NewIterator();
  it->SetSkipEmptyNodes(0);
  it->InitTraversal();

  int num_blocks = 2; // only look at first 2 blocks
  for (unsigned int j = 0; j < num_blocks; ++j)
    {	
      std::cout<<".... block: "<<j <<std::endl;
      // pass to vtk
      vtkDataSet* ds = dynamic_cast<vtkDataSet*>(it->GetCurrentDataObject());
      if (!ds)
	{
	  SENSEI_ERROR("... printMe() Failed to get block " << j);
	  return;
	}

      vtkDataSetAttributes* dsa = dynamic_cast<vtkDataSetAttributes*>(ds->GetCellData());
	
      vtkDataArray* array = dsa->GetArray(0);	
    
      vtkDoubleArray* doubleA = vtkDoubleArray::SafeDownCast(dsa->GetArray(0));
      std::cout<<"   get value? "<<doubleA->GetValue(0)<<std::endl;
      //array->PrintSelf(std::cout, vtkIndent(10));

      void* data = array->GetVoidPointer(0);
      std::cout<<"   peek: "<<((double*)data)[0]<<", "<<((double*)data)[1]<<std::endl;
      // next block                                                                                                                                                                                     
      it->GoToNextItem();
    }
  
  it->Delete();
  */
  int numBlocks = m_Mesh->GetNumberOfBlocks(); 
  for (int i=0; i<numBlocks; i++) {
    //std::cout<<" block : "<<i<<std::endl;
    vtkDataObject* curr = m_Mesh->GetBlock(i);
    if (!curr)
      continue;

    vtkDataSet* ds = dynamic_cast<vtkDataSet*>(curr);
    //curr->PrintSelf(std::cout, vtkIndent(10));
    //ds->PrintSelf(std::cout, vtkIndent(10));
    vtkDataSetAttributes* dsa = dynamic_cast<vtkDataSetAttributes*>(ds->GetCellData()); 
    //dsa->PrintSelf(std::cout, vtkIndent(5)); 

    for (int s=0; s<6; s++) {
      vtkDoubleArray* doubleA = vtkDoubleArray::SafeDownCast(dsa->GetArray(s));
      if (doubleA == 0) 
	break;
      /*
      std::cout<<" ARRAY: "<<s<<"\t";
      for (int k=0; k<6; k++)
	std::cout<<" "<<doubleA->GetValue(k);            
      std::cout<<std::endl;
      */
    }
  }
}


int VTKBeamBeam3DWriter::turnEnd()
{
  m_DataAtTurn->SetDataObject(m_MeshName, m_Mesh);

  m_Writer->Execute(m_DataAtTurn);

  m_DataAtTurn->ReleaseData();
  m_DataAtTurn->Delete();

  return 1;
}


/*   
int VTKBeamBeam3DWriter::setupPhaseout(int nBunchs, uint64_t total, uint64_t start, uint64_t count)
{}

int VTKBeamBeam3DWriter::phaseout(double* data, int currbunch, int itic, int64_t iseedinit, double* close2g)
{}
*/



//
//
//
VTKBeamBeam3DReader::VTKBeamBeam3DReader(MPI_Comm comm)
  :m_Comm(comm)
{
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &m_Rank);
  MPI_Comm_size(comm, &m_Size);

  m_Reader = sensei::HDF5DataAdaptor::New();
}


VTKBeamBeam3DReader::VTKBeamBeam3DReader(MPI_Comm comm, const char* name)
  :m_Comm(comm)
{  
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  //MPI_Init(NULL, NULL);
  MPI_Comm_rank(comm, &m_Rank);
  MPI_Comm_size(comm, &m_Size);

  if (m_Rank == 0) 
    std::cout<<"name = "<<name<<std::endl;


  bool doStreaming = false;
  bool doCollective = false;
  m_Reader  = H5DataAdaptorPtr::New();

  m_Reader->SetCommunicator(comm);
  m_Reader->SetStreaming(doStreaming);
  m_Reader->SetCollective(doCollective);
  m_Reader->SetStreamName(name);
}

VTKBeamBeam3DReader::~VTKBeamBeam3DReader()
{
  m_Reader->Delete();
}

bool VTKBeamBeam3DReader::isValidTurn(uint32_t t, uint32_t base) 
{
  //return ((t >= m_SelectionTurnStart) && (t-m_SelectionTurnStart < m_SelectionTurnCount));
  uint32_t adjustedStart = m_SelectionTurnStart + base;
  return ((t >= adjustedStart ) && (t - adjustedStart < m_SelectionTurnCount));
}


bool  VTKBeamBeam3DReader::InitializeXML(const char* filename)
{
  pugi::xml_document doc;
  if (sensei::XMLUtils::Parse(m_Comm, filename, doc))
    {
      SENSEI_ERROR("Failed to load, parse, and share XML configuration")
	MPI_Abort(m_Comm, -1);
      return false;
    }

  pugi::xml_node root = doc.child("sensei").child("reader");
  
  if (!root) {
    if (m_Rank == 0)
      std::cout<<" !ERROR: Expecting xml node <reader> under root <sensei>"<<std::endl;
    return false;
  }

  if (m_Reader->Initialize(root))
    return false;

  pugi::xml_node cols = root.child("cols");
  if (cols) {    
    pugi::xml_attribute v1 = cols.attribute("value1");
    pugi::xml_attribute v2 = cols.attribute("value2");

    if (m_Rank == 0) std::cout<<"Input columns:\t"<<v1.value()<<" "<<v2.value()<<std::endl;
    setAttr(v1.value(), v2.value());
  }

  pugi::xml_node bunch = root.child("bunch");
  if (bunch) {
    pugi::xml_attribute v = bunch.attribute("value");
    std::stringstream peek(v.value());
    uint32_t bunchID=0; peek >> bunchID;
    if (m_Rank == 0) std::cout<<"Input bunch:\t"<<bunchID<<std::endl;
    setBunchID(bunchID);
  }

  pugi::xml_node turn = root.child("turn");
  if (turn) {
    pugi::xml_attribute vs = turn.attribute("start");
    pugi::xml_attribute vc = turn.attribute("count");
    std::stringstream p1(vs.value());
    std::stringstream p2(vc.value());
    uint32_t tStart=0; p1 >> tStart;
    uint32_t tCount=0; p2 >> tCount;

    if (m_Rank == 0) std::cout<<"Input turnID:\t"<<tStart<<" "<<tCount<<std::endl;
    setTurnRange(tStart, tCount);
  }

  pugi::xml_node ptl = root.child("ptl");
  if (ptl) {
    pugi::xml_attribute vs = ptl.attribute("start");
    pugi::xml_attribute vc = ptl.attribute("count");
    std::stringstream p1(vs.value());
    std::stringstream p2(vc.value());
    uint64_t pStart=0; p1 >> pStart;
    uint64_t pCount=0; p2 >> pCount;

    if (m_Rank == 0) std::cout<<"Input pID:\t"<<pStart<<" "<<pCount<<std::endl;
    setParticleRange(pStart, pCount);
  }

  pugi::xml_node debug = root.child("debug");
  if (debug)
    m_Debug = true;

  return true;
}


bool VTKBeamBeam3DReader::setParticleRange(uint64_t pStart, uint64_t pSize)
{
  m_SelectionParticleStart = pStart;
  m_SelectionParticleCount = pSize;

  return true;
}

bool VTKBeamBeam3DReader::setBunchID(uint32_t bunch)
{
  m_SelectionBunchID = bunch;

  return true;
}

bool VTKBeamBeam3DReader::setTurnRange(uint32_t tStart, uint32_t tSize)
{
  if (tSize < 1) 
    return false;

  m_SelectionTurnStart = tStart;
  m_SelectionTurnCount = tSize;

  return true;
}


bool VTKBeamBeam3DReader::setAttr(const char* n1, const char* n2)
{
  m_SelectionColName1 = n1;
  m_SelectionColName2 = n2;

  return true;
}


uint64_t VTKBeamBeam3DReader::readMultiBlock(vtkMultiBlockDataSet* ds, sensei::MeshMetadataPtr mmd)
{
  if (ds == NULL)
    return 0;

#ifdef MULTI_PROCESSOR
  // have to use iterator as numBlocks is does not reflect with correct behaviour
  // with multi-processor
  // num_blocks does not reflect what this rank has

  vtkCompositeDataIterator* it = ds->NewIterator();
  while (!it->IsDoneWithTraversal()) {
    vtkDataObject* curr = it->GetCurrentDataObject();

    vtkDataSet* bds = vtkDataSet::SafeDownCast(curr);
    if (bds == NULL) {
      std::cerr << " surprise! " << std::endl;
      break;
    }
    unsigned int n_arrays = mmd->NumArrays;
    for (unsigned int j = 0; j < n_arrays; j++) {
      vtkDataArray* array;
      if (mmd->ArrayCentering[j] == vtkDataObject::POINT) {
	array = bds->GetPointData()->GetArray(mmd->ArrayName[j].c_str());
      } else {
	array = bds->GetCellData()->GetArray(mmd->ArrayName[j].c_str());
      }

      if ((m_Rank == 0) && m_Debug)       
	std::cout<<"=> idx:"<<j<<" array name="<<mmd->ArrayName[j]<<" size="<<array->GetNumberOfTuples()<<std::endl;
	vtkDoubleArray* doubleA = vtkDoubleArray::SafeDownCast(array);
                                                                                                                             
	for (int k=0; k<6;k++)
	  std::cout<<" "<<doubleA->GetValue(k);            
	std::cout<<std::endl;	
      }
    }
    it->GoToNextItem();
  } // while it                                                                                                                                                                                       
  it->Delete();

#else

  int numBlocks = ds->GetNumberOfBlocks();
  //std::cout<<" NumBlocks = "<<numBlocks<<std::endl;

  int particleCounter = 0;

  bool firstHit = true;
  for (int i=0; i<numBlocks; i++) {
    //std::cout<<"=> block : "<<i<<" num cells:"<<mmd->BlockNumCells[i];
    auto ext = mmd->BlockExtents[i];
    //std::cout<<" bunch="<<ext[0]<<" particle range: "<<ext[2]<<" to "<<ext[3]<<std::endl;
    
    if (ext[0] != m_SelectionBunchID)
      continue;
    
    auto selEnd = m_SelectionParticleStart + m_SelectionParticleCount;
    if ((ext[3] <= m_SelectionParticleStart) || (ext[2] >= selEnd))
      continue;

    auto s = 0, t=0; // s=start index, t=end index (includsive)
    if (m_SelectionParticleStart >= ext[2]) 
      s = m_SelectionParticleStart - ext[2];
    else 
      s = 0;
    
    if (selEnd >= ext[3])
      t = mmd->BlockNumCells[i]-1;
    else 
      t = selEnd - ext[2] - 1;   
       
    //auto t = std::min(m_SelectionParticleCount, mmd->BlockNumCells[0]);

    vtkDataObject* curr = ds->GetBlock(i);
    if (!curr)
      return 0;
    
    vtkDataSet* ds = dynamic_cast<vtkDataSet*>(curr);
    vtkDataSetAttributes* dsa = dynamic_cast<vtkDataSetAttributes*>(ds->GetCellData()); 

    //for (int a=0; a<2; a++) { // b/c only 2 specified attrs are  read,
    vtkDoubleArray* doubleA = vtkDoubleArray::SafeDownCast(dsa->GetArray(0));
    vtkDoubleArray* doubleB = vtkDoubleArray::SafeDownCast(dsa->GetArray(1));

    if ((doubleA == NULL) || (doubleB == NULL)) {      
      std::cerr<<"ERROR:  Unable to get one of the arrays. "<<std::endl;
      return 0;
    }

    if (m_Debug) {
      for (auto k=0; k<mmd->BlockNumCells[i]; k++) {
	if (k < 5) {
	  std::cout<<k<<": "<<doubleA->GetValue(k);
	  std::cout<<" "<<doubleB->GetValue(k)<<std::endl;            
	} else {
	  std::cout<<"..."<<std::endl;
	  break;
	}
      }     
      std::cout<<std::endl;   
    }

    for (auto k=s; k<=t; k++) {
      *m_data1 = doubleA->GetValue(k); m_data1++;
      *m_data2 = doubleB->GetValue(k); m_data2++;
      /*
      if (!firstHit) {
	m_data1 ++; 
	m_data2 ++;
      } 
      *m_data1 = doubleA->GetValue(k);
      *m_data2 = doubleB->GetValue(k);
      firstHit = false;
      */
    }    
    particleCounter += (t-s+1);       
  }

  return particleCounter;

#endif

}

void VTKBeamBeam3DReader::showError(const std::string& msg)
{
  if (m_Rank == 0) 
    std::cout<<msg<<std::endl;  
}

bool VTKBeamBeam3DReader::readPtl(double* data1, double* data2)
{
  if ((data1 == NULL) || (data2 == NULL)) {
    showError("Please allocate data first.");
    return false;
  }

  if (m_SelectionColName1.compare(m_SelectionColName2) == 0) {
    showError("Please select which attrs to read.");
    return false;
  }
  if (0 == m_SelectionParticleCount) {
    showError("Please select number of particles to read.");
    return false;
  }
  if (0 == m_SelectionTurnCount) {
    showError("Please select number of turns to read.");
    return false;
  }
	
  if (1 == m_Size) {
    //return readPtlSerial(data1, data2);
    m_data1 = data1;
    m_data2 = data2;

    bool result = readPtlSerial();
    //std::cout<<"turns read: "<<m_TurnsRead<<std::endl;
    return result;
  } else {
    std::cout<<"\n\nWill need to assemble blocks from different ranks and map back to data"<<std::endl;
    return false;
  }

  return false;
}

bool VTKBeamBeam3DReader::readPtlSerial() //double* data1, double* data2)
{  
  m_Reader->OpenStream();
  uint32_t base = m_Reader->GetDataTimeStep();
  m_TurnsRead=0;
  while (true) {    
      double t = m_Reader->GetDataTime();
      int it = m_Reader->GetDataTimeStep();
      if ((m_Rank == 0) && m_Debug)
	std::cout << "\n===> Received step: " << it << " time: " << t<< ". Assigned turn range: "<<m_SelectionTurnStart<<", "<<m_SelectionTurnCount<<std::endl;      

      if (isValidTurn(it, base)) {      
	m_TurnsRead ++;
	unsigned int nMeshes;
	m_Reader->GetNumberOfMeshes(nMeshes);
	
	if (m_Debug) std::cout<<"\t has mesh: "<<nMeshes<<std::endl;
	unsigned int i = 0;
	while (i < nMeshes) { // uh, we only have one mesh with this app
	  sensei::MeshMetadataPtr mmd = sensei::MeshMetadata::New();
	  
	  if (m_Reader->GetMeshMetadata(i, mmd)) {
	    std::cerr << "Unable to get metadata" << std::endl;
	    break;
	  }
	  
	  std::string meshName = mmd->MeshName;
	  
	  vtkDataObject* mesh = nullptr;
	  m_Reader->GetMesh(meshName, false, mesh);
	  if (mesh == NULL) break;	  
	  
	  unsigned int n_arrays = mmd->NumArrays; 
	  
	  for (unsigned int j = 0; j < n_arrays; j++) {	  
	    std::string arrayName = mmd->ArrayName[j];
	    
	    if ((arrayName.compare(m_SelectionColName1) == 0) || (arrayName.compare(m_SelectionColName2) == 0))
	      {
		if (m_Debug) std::cout<<"reading in array: "<<arrayName<<std::endl;
		m_Reader->AddArray(mesh, meshName, mmd->ArrayCentering[j], arrayName);			       
	      }	
	  }
	
	  vtkMultiBlockDataSet* ds = vtkMultiBlockDataSet::SafeDownCast(mesh);	
	  
	  uint64_t numPtlRead = readMultiBlock(ds, mmd);
	  m_NumPtlPerTurn = numPtlRead;
	  
	  i += 1;
	  mesh->Delete();
	} // while mesh
	
	m_Reader->ReleaseData();	
      } // validTurn
      else {
	if (m_Debug) std::cout<<" skipping turn: "<<it<<std::endl;
      }

      if (m_Reader->AdvanceStream())
	break;
  }
  
  m_Reader->CloseStream();

  return true;
}

//
//
//
/*
VTKBeamBeam3D::VTKBeamBeam3D(MPI_Comm comm, const char* name)
{  
  
}


int VTKBeamBeam3D::bar(int c) const{
    return a + c;
}

int VTKBeamBeam3D::baz(double* d) const{
  cout << "input[0] = "<<*d << endl;;
  return 14;
}
*/










