enable_language(Fortran)

set(sources 
	    mod_tunefoot.f90 main.f90)


set(exe tunefoot_sensei)
add_executable(${exe} ${sources} )


target_include_directories(${exe} PRIVATE ${SRC_FORTRAN})

target_link_libraries(${exe} PRIVATE  ${LIB_FORTRAN}  ${MPI_mpi_mpifh_LIBRARY}) 

SET_TARGET_PROPERTIES(${exe} PROPERTIES LINKER_LANGUAGE Fortran)



		  

