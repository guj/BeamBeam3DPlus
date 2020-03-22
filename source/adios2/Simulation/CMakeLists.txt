#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
enable_language(Fortran)

set(base_source_dir ../../base/)
set(base_sources 		
		 ${base_source_dir}/Beambeam.f90
		 ${base_source_dir}/Commun.f90 	
		 ${base_source_dir}/CompDom.f90	 
		 ${base_source_dir}/Dampfluc.f90			
		 ${base_source_dir}/DepoScat.f90 	
		 ${base_source_dir}/Distribution.f90 
		 ${base_source_dir}/Extmap2nd.f90
		 ${base_source_dir}/Extmap3rd.f90	
		 ${base_source_dir}/Extmap4th.f90	
		 ${base_source_dir}/Feedback.f90  
		 ${base_source_dir}/FFT.f90
		 ${base_source_dir}/FieldSolver.f90		
		 ${base_source_dir}/Input.f90  
		 ${base_source_dir}/Linearmap.f90  
		 ${base_source_dir}/Orbit.f90  
		 ${base_source_dir}/Pgrid.f90  
		 ${base_source_dir}/Timer.f90  
		 ${base_source_dir}/Transfer.f90
		 ${base_source_dir}/Transpose.f90 
		 ${base_source_dir}/Utility.f90)
   		 	     		

#include_directories(${ADIOS2_INCLUDE_DIR})


list(APPEND base_sources AccSimulator.f90)
list(APPEND base_sources IoAdios2.f90)
list(APPEND base_sources ${base_source_dir}/main.f90)
list(APPEND base_sources Output.f90)

message(STATUS "..${ADIOS2_LIBRARIES}")

set(exe xmain_adios2)

add_executable(${exe} ${base_sources} )
#target_include_directories(${exe} PRIVATE ${ADIOS2_INCLUDE_DIR}/fortran)
target_include_directories(${exe} PRIVATE ${ADIOS2_INCLUDE_DIR_FORTRAN})

target_link_libraries(${exe} PRIVATE ${ADIOS2_LIBRARIES} ${MPI_mpi_mpifh_LIBRARY})
SET_TARGET_PROPERTIES(${exe}  PROPERTIES LINKER_LANGUAGE Fortran)

#!target_link_libraries(${exe} PRIVATE ${ADIOS2_LIBRARIES} MPI::MPI_Fortran) 
#!SET_TARGET_PROPERTIES(${exe}  PROPERTIES LINKER_LANGUAGE Fortran)


		  

