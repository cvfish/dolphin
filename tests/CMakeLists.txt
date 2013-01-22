# The CMakeLists file for Light-Matrix (Testing)

cmake_minimum_required(VERSION 2.8)
project(LightMatrix)

enable_testing()

set(LIGHTTEST_HOME $ENV{LIGHT_TEST_HOME})
if (LIGHTTEST_HOME)
message(STATUS "[Dolphin] Light-Test: ${LIGHTTEST_HOME}")
else (LIGHTTEST_HOME)
message(FATAL_ERROR "LIGHT_TEST_HOME is not set")
endif (LIGHTTEST_HOME)

set(LIGHTMAT_HOME $ENV{LIGHT_MATRIX_HOME})
if (LIGHTMAT_HOME)
message(STATUS "[Dolphin] Light-Matrix: ${LIGHTMAT_HOME}")
else (LIGHTMAT_HOME)
message(FATAL_ERROR "LIGHT_MATRIX_HOME is not set")
endif (LIGHTMAT_HOME)

include_directories(..)
include_directories(${LIGHTTEST_HOME})
include_directories(${LIGHTMAT_HOME})

include("../cmake_modules/CompilerConfig.cmake")

#==========================================================
#
#    third-party library dependencies
#
#==========================================================

set(CMAKE_MODULE_PATH "../cmake_modules" ${CMAKE_MODULE_PATH})

# Intel library 

find_package(ICCLib)
if (ICCLIB_FOUND)
message(STATUS "[Dolphin] ICC Library found: ${ICCLIB_PATH}")
else (ICCLIB_FOUND)
message(FATAL_ERROR "[Dolphin] ICC Library not found")
endif (ICCLIB_FOUND)

# BLAS and Lapack

find_package(MKL)
if (MKL_FOUND)
message(STATUS "[LMAT] Intel MKL found: ${MKLROOT_PATH}")
else (MKL_FOUND)
message(FATAL_ERROR "[LMAT] Intel MKL not found")
endif (MKL_FOUND)


#==========================================================
#
#    Header file groups (to be used as dependencies)
#
#==========================================================

set(INC ../dolphin)
set(EXECUTABLE_OUTPUT_PATH bin)

set(COMMON_BASE_HS
    ${INC}/common/common_base.h
    ${INC}/common/import_lmat.h
    ${INC}/common/properties.h)
    
set(COMMON_TOOLS_HS
    ${INC}/common/dpaccum.h
    ${INC}/common/common_calc.h)

set(COMMON_HS
    ${COMMON_BASE_HS}
    ${COMMON_TOOLS_HS})
    
    
#==========================================================
#
#    Unit testing targets
#
#==========================================================    
    
# main entry

add_library(test_main STATIC test_main.cpp)

# common module

set(COMMON_TEST_HS
    ${COMMON_HS})

add_executable(test_dpaccum ${COMMON_TEST_HS} common/test_dpaccum.cpp)
add_executable(test_common_calc ${COMMON_TEST_HS} common/test_common_calc.cpp)
add_executable(test_metrics ${COMMON_TEST_HS} common/test_metrics.cpp)

set(COMMON_TESTS
    test_dpaccum
    test_common_calc
    test_metrics)

# all tests

set(DOLPHIN_TESTS_USING_LINALG
    test_metrics)

set(DOLPHIN_ALL_TESTS
    ${COMMON_TESTS})


#==========================================================
#
#    treatments to all tests
#
#==========================================================

# Link to test_main
	
foreach(tname ${DOLPHIN_TESTS_USING_LINALG})
    target_link_libraries(${tname} ${MKL_LIBRARY})
endforeach(tname)	
	
foreach(tname ${DOLPHIN_ALL_TESTS})
	target_link_libraries(${tname} test_main)
    set_target_properties(${tname}
        PROPERTIES
        COMPILE_FLAGS "-DLMAT_USE_INTEL_SVML")
    target_link_libraries(${tname} ${SVML_LIBRARY})
endforeach(tname)	

# Add tests

foreach(tname ${DOLPHIN_ALL_TESTS})
    
	add_test(NAME ${tname} COMMAND ${tname})
endforeach(tname)



    