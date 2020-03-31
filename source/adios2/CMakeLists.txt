cmake_minimum_required(VERSION 3.6)

find_package(ADIOS2 REQUIRED)

  add_library(sADIOS2 INTERFACE)
  target_link_libraries(sADIOS2 INTERFACE adios2::adios2_f)
  install(TARGETS sADIOS2 EXPORT sADIOS2)
  install(EXPORT sADIOS2 DESTINATION lib/cmake EXPORT_LINK_INTERFACE_LIBRARIES)



message(STATUS "adios2 found.")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(Simulation)
add_subdirectory(Analysis)
