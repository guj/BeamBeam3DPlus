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


target_link_libraries(${flib} PRIVATE ${senseiC})

#target_link_libraries(${flib} PRIVATE MPI::MPI_Fortran MPI::MPI_C)
target_link_libraries(${flib} PRIVATE MPI::MPI_C ${MPI_mpi_mpifh_LIBRARY})



#########################
### for testing purpose:
#########################

set(test_source test/simpleTest.f90)

set(test_exe simpleTest)

add_executable(${test_exe}  ${test_source})

target_link_libraries(${test_exe} PRIVATE ${flib} ${MPI_mpi_mpifh_LIBRARY}) 
