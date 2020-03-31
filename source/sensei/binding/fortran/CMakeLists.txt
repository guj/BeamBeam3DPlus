cmake_minimum_required(VERSION 3.10)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX QUIET)

set(F2C
     ${CMAKE_CURRENT_SOURCE_DIR}/f2c/sensei_f2c_beam.cpp
)


set(MODULES
     ${CMAKE_CURRENT_SOURCE_DIR}/modules/sensei_def.f90	
     ${CMAKE_CURRENT_SOURCE_DIR}/modules/sensei_w_mod.f90	
     ${CMAKE_CURRENT_SOURCE_DIR}/modules/sensei_r_mod.f90	
)

#set(flib senseiBeam_f)

add_library(${flib} SHARED ${MODULES} ${F2C})


target_include_directories(${flib}
  PUBLIC
    ${SENSEI_C_INCLUDE}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fortran>
)


#target_link_libraries(${flib} PRIVATE ${senseiC})
#target_link_libraries(${flib} PRIVATE MPI::MPI_C ${MPI_mpi_mpifh_LIBRARY})




set_target_properties(${flib} PROPERTIES LINKER_LANGUAGE Fortran)
#  target_link_libraries(${tgt} sADIOS2) 

    target_compile_definitions(${flib} PRIVATE USE_MPI)
#target_link_libraries(${flib} PRIVATE MPI::MPI_C MPI::MPI_Fortran sSENSEI)
target_link_libraries(${flib} PRIVATE MPI::MPI_C MPI::MPI_Fortran ${senseiC})




#########################
### for testing purpose:
#########################

set(simpleTest_source test/simpleTest.f90)


macro(fortran_add_test_helper testname)
  set(tgt sensei.${testname})
  add_executable(${tgt} ${simpleTest_source})
  set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${tgt} ${flib}) 

    target_compile_definitions(${tgt} PRIVATE USE_MPI)
    target_link_libraries(${tgt} MPI::MPI_Fortran)
endmacro()

fortran_add_test_helper(simpleTest)
