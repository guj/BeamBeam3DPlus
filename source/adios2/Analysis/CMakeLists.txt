#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
enable_language(Fortran)


set(sources 
	    main.f90
	    mod_tunefoot.f90)

message(STATUS "..${ADIOS2_LIBRARIES}")

set(exe tunefoot_adios2)
add_executable(${exe} ${sources} )
target_include_directories(${exe} PRIVATE ${ADIOS2_INCLUDE_DIR_FORTRAN})
#!target_link_libraries(${exe} PRIVATE ${ADIOS2_LIBRARIES} MPI::MPI_Fortran) 
#!SET_TARGET_PROPERTIES(${exe} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${exe} PRIVATE ${ADIOS2_LIBRARIES} ${MPI_mpi_mpifh_LIBRARY})
SET_TARGET_PROPERTIES(${exe}  PROPERTIES LINKER_LANGUAGE Fortran)



		  

