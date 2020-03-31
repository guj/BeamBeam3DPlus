#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

macro(fortran_add_test_helper testname mpi)
  set(tgt adios2.${testname})
  add_executable(${tgt} ${testname}.f90 main.f90)
  set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${tgt} sADIOS2)
  if(mpi)
    target_compile_definitions(${tgt} PRIVATE USE_MPI)
    target_link_libraries(${tgt} MPI::MPI_Fortran)
  else()
    add_test(NAME ${pfx}${testname} COMMAND ${tgt})
  endif()
endmacro()


set(mpi TRUE)

fortran_add_test_helper(mod_tunefoot ${mpi})



