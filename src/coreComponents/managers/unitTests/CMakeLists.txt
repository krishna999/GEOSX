#
# Specify list of tests
#

set( gtest_geosx_tests
     testRecursiveFieldApplication.cpp
     testMeshGeneration.cpp
     testFunctions.cpp
   )


set( dependencyList gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core conduit)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} conduit)
endif()

if ( ENABLE_MPI )
  set ( dependencyList ${dependencyList} mpi )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${dependencyList}
            )

    blt_add_test( NAME ${test_name}
            COMMAND ${test_name}
            )

endforeach()
