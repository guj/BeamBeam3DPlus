cmake_minimum_required(VERSION 3.10)
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
   		 	     		


list(APPEND base_sources AccSimulator.f90)
list(APPEND base_sources ${base_source_dir}/main.f90)
list(APPEND base_sources Output.f90)

message(STATUS "source: ${base_sources}")
set(exe xmain_sensei)

add_executable(${exe} ${base_sources} )
#target_include_directories(${exe} PRIVATE ${SRC_BEAM})
target_include_directories(${exe} PRIVATE ${SRC_FORTRAN})
message(STATUS  "... fortran include: ${SRC_FORTRAN}")
#
# not using MPI::MPI_Fortran but only use the MPI_mpi_mpifh_LIBRARY
# because the -flat_namespace would cause link error
# and I do not know how not to set it.
# 
#target_link_libraries(${exe} PRIVATE ${LIB_BEAM_W} ${MPI_mpi_mpifh_LIBRARY}) 
target_link_libraries(${exe} PRIVATE ${LIB_FORTRAN} ${MPI_mpi_mpifh_LIBRARY}) 
SET_TARGET_PROPERTIES(${exe}  PROPERTIES LINKER_LANGUAGE Fortran)


		  

