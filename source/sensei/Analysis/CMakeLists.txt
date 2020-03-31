enable_language(Fortran)

set(sources 
	    mod_tunefoot.f90 main.f90)

set(tgt sensei.tunefoot)
add_executable(${tgt} ${sources})
set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${tgt} ${flib}) 

target_compile_definitions(${tgt} PRIVATE USE_MPI)
target_link_libraries(${tgt} MPI::MPI_Fortran)








		  

