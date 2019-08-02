#Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> 
cmake_minimum_required(VERSION 2.8)
project(mgm)


## SETUP ##

# SETUP A MAKE CLEAN TO CLEAN EVERYTHING
set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
      "CMakeCache.txt;CMakeFiles;Makefile;cmake_install.cmake")

# DEFAULT BUILD: RELEASE
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()

if (NOT TARGET IIOLIB)
   include_directories(extern/iio)
   add_subdirectory(extern/iio)
endif() 
SET( LIBS ${LIBS} IIOLIB)

FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
   message("[OPENMP FOUND]")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


SET( LIBS ${LIBS} fftw3)
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
else ()
  set (CMAKE_C_STANDARD 99)
endif ()


# Build as ROS package also, if possible
find_package(catkin QUIET)
if (catkin_FOUND)
  message("[CATKIN FOUND] Configuring for ROS library compatible compilation!")
  include_directories(
    ${catkin_INCLUDE_DIRS}
  )

  # Specify export settings
  catkin_package(
    INCLUDE_DIRS include src/impl extern/iio
    LIBRARIES mgm iio
  )
else()
  message("[CATKIN NOT FOUND] Skipping ROS relevant compile steps...")
  message("(This is ok if you're not using MGM for ROS!)")
endif()

# INCLUDES AND COMPILATION ##

# Specify build include directories
include_directories(
  include/mgm
  include
  src/impl
)

add_library(mgm
  src/impl/img.cc
  src/impl/point.cc
  src/impl/mgm_costvolume.cc
  src/impl/mgm_core.cc
  src/impl/census_tools.cc
  src/impl/shear.c
  src/impl/stereo_utils.cc
  src/impl/mgm_multiscale.cc
  src/impl/smartparameter.cc
)
target_link_libraries(mgm ${LIBS})

add_library(iio extern/iio/iio.c)
target_link_libraries(iio IIOLIB)

add_executable(run_mgm_multi src/main_mgm_multi.cc)
target_link_libraries(run_mgm_multi ${LIBS} mgm)

add_executable(run_mgm src/main_mgm.cc)
target_link_libraries(run_mgm ${LIBS} mgm)

if (catkin_FOUND)
  install(DIRECTORY include/${PROJECT_NAME} src/impl extern/iio
    DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})
endif()
