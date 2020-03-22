cmake_minimum_required(VERSION 3.10)

set(senseiC senseiBeam_c)

set(senseiC_sources 	
		    ${CMAKE_CURRENT_SOURCE_DIR}/c/sensei_beam_c_api.cpp
		    ${CMAKE_CURRENT_SOURCE_DIR}/c/VTKBeamBeam3D.cpp			
)

add_library(${senseiC} SHARED ${senseiC_sources} )

target_include_directories(${senseiC} PRIVATE
                                      ${SENSEI_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
find_package(MPI)
target_link_libraries(${senseiC} PRIVATE
                                 ${SENSEI_LIBRARIES} MPI::MPI_C sVTK  sHDF5 sADIOS2)

set(SENSEI_C_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/c)

add_subdirectory(fortran)



