#Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> 
cmake_minimum_required(VERSION 2.8)


# SETUP A MAKE CLEAN TO CLEAN EVERYTHING
#
#
set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
      "CMakeCache.txt;CMakeFiles;Makefile;cmake_install.cmake")

# DEFAULT BUILD: RELEASE
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()

if (NOT TARGET IIOLIB)
   include_directories("iio")
   add_subdirectory(iio)
endif() 
SET( LIBS ${LIBS} IIOLIB)

FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
   message("OPENMP FOUND")
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

#
add_executable (mgm_multi main_mgm_multi.cc img.cc point.cc mgm_costvolume.cc mgm_core.cc shear.c census_tools.cc stereo_utils.cc stereo_utils.h mgm_multiscale.cc mgm_multiscale.h)
target_link_libraries(mgm_multi ${LIBS})

add_executable (mgm main_mgm.cc mgm_multiscale.cc img.cc point.cc mgm_costvolume.cc mgm_core.cc census_tools.cc shear.c stereo_utils.cc stereo_utils.h)
target_link_libraries(mgm ${LIBS})

