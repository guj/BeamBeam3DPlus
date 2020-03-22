cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SENSEI_ROOT "$ENV{SENSEI_ROOT}")
set(VTK_DIR "$ENV{VTK_DIR}")

#find_package(HDF5 REQUIRED)
set(ENABLE_HDF5 true)
include(hdf5)

set(ENABLE_ADIOS2 true)
include(adios2)

if (SENSEI_ROOT)
   set(SENSEI_INCLUDE_OPTS HINTS ${SENSEI_ROOT}/include NO_DEFAULT_PATHS)
   set(SENSEI_LIBRARY_OPTS
      HINTS ${SENSEI_ROOT}/lib ${SENSEI_ROOT}/lib64
      NO_DEFAULT_PATHS
    )
endif()

if (VTK_DIR)
   include(vtk)
   set(VTK_INCLUDE_OPTS HINTS ${VTK_DIR}/include ${VTK_DIR}/include/vtk-8.1  NO_DEFAULT_PATHS)
endif()


find_path(SENSEI_INCLUDE_PATH senseiConfig.h ${SENSEI_INCLUDE_OPTS})
find_library(SENSEI_LIBRARY_CPP senseiCore ${SENSEI_LIBRARY_OPTS})


set(SENSEI_INCLUDE_DIR ${SENSEI_INCLUDE_PATH}
		       ${SENSEI_INCLUDE_PATH}/pugixml
		       )


set(SENSEI_LIBRARIES 
		     ${SENSEI_LIBRARY_CPP}
)

message(STATUS "... using sensei libs: ${SENSEI_LIBRARIES}  ${VTK_INCLUDE_DIRS}")
#add_library(sensei::sensei UNKNOWN IMPORTED)

#set_target_properties(sensei::sensei PROPERTIES
#        			     IMPORTED_LOCATION             "${ADIOS2_LIBRARY}"
#        			     INTERFACE_INCLUDE_DIRECTORIES "${ADIOS2_INCLUDE_DIR}"
#      )	



set(flib senseiBeam_f)
set(SRC_FORTRAN "${CMAKE_CURRENT_BINARY_DIR}/binding/fortran")
set(LIB_FORTRAN ${flib})

add_subdirectory(binding)
add_subdirectory(Analysis)
add_subdirectory(Simulation)

