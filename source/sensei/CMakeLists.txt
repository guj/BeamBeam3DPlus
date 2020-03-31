cmake_minimum_required(VERSION 3.6)

find_package(HDF5 REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(SENSEI REQUIRED)

  add_library(sSENSEI INTERFACE)
  target_link_libraries(sSENSEI INTERFACE senseiCore)
  #install(TARGETS sADIOS2 EXPORT sADIOS2)
  #install(EXPORT sADIOS2 DESTINATION lib/cmake EXPORT_LINK_INTERFACE_LIBRARIES)



message(STATUS "sensei added")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(flib senseiBeam_f)
set(SRC_FORTRAN "${CMAKE_CURRENT_BINARY_DIR}/binding/fortran")
set(LIB_FORTRAN ${flib})

add_subdirectory(binding)
add_subdirectory(Simulation)
add_subdirectory(Analysis)






