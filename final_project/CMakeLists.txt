cmake_minimum_required(VERSION 3.10)

PROJECT(final_project)
SET(VTK_DIR c:/MyProjects/VTK-bin)
SET(CMAKE_VERBOSE_MAKEFILE ON)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE} ${VTK_DIR}/Rendering)

add_executable(final_project final_project)


target_link_libraries(final_project glu32)
target_link_libraries(final_project opengl32)
if(VTK_LIBRARIES)
target_link_libraries(final_project ${VTK_LIBRARIES})
else()
target_link_libraries(final_project vtkHybrid)
endif()
