#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


set(base_sources 		
		 Beambeam.f90
		 Commun.f90 	
		 CompDom.f90	 
		 Dampfluc.f90			
		 DepoScat.f90 	
		 Distribution.f90 
		 Extmap2nd.f90
		 Extmap3rd.f90	
		 Extmap4th.f90	
		 Feedback.f90  
		 FFT.f90
		 FieldSolver.f90		
		 Input.f90  
		 Linearmap.f90  
		 Orbit.f90  
		 Pgrid.f90  
		 Timer.f90  
		 Transfer.f90
		 Transpose.f90 
		 Utility.f90)
   		 	     		


#add_library(beambeamCompute ${base_sources})
#SET_TARGET_PROPERTIES(beambeamCompute PROPERTIES LINKER_LANGUAGE Fortran)
#target_link_libraries(beambeamCompute PRIVATE MPI::MPI_Fortran)

list(APPEND base_sources AccSimulator.f90)
list(APPEND base_sources main.f90)
list(APPEND base_sources Output.f90)

add_executable(xmain ${base_sources} )
target_link_libraries(xmain PRIVATE MPI::MPI_Fortran)
SET_TARGET_PROPERTIES(xmain PROPERTIES LINKER_LANGUAGE Fortran)

#install(TARGETS xmain RUNTIME DESTINATION bin)
		  

