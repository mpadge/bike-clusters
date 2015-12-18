cmake_minimum_required(VERSION 2.8.7)
project(bike-clusters CXX)
include_directories(include)
ADD_DEFINITIONS(
    -std=c++11 
)
set(SOURCES_NEUTRAL src/Utils.c++ src/ClusterData.c++ src/ClusterCalculations src/mainNeutral.c++)
set(SOURCES_ACTUAL src/Utils.c++ src/ClusterData.c++ src/ClusterCalculations src/mainActual.c++)
set(SOURCES_RANDOM src/Utils.c++ src/randomHierarchy.c++)

add_executable(clustersNeutral ${SOURCES_NEUTRAL})
add_executable(clustersActual ${SOURCES_ACTUAL})
add_executable(rhier ${SOURCES_RANDOM})
target_link_libraries(clustersNeutral boost_program_options CGAL gmp)
target_link_libraries(clustersActual boost_program_options CGAL gmp)
target_link_libraries(rhier boost_program_options CGAL gmp)
