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
message(STATUS "[Dolphin] ICC Library not found")
endif (ICCLIB_FOUND)

set(SVML_FOUND ICCLIB_FOUND)

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
    ${INC}/common/dpaccum.h)

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

set(COMMON_TESTS
    test_dpaccum)

# all tests

set(DOLPHIN_ALL_TESTS
    ${COMMON_TESTS})


#==========================================================
#
#    treatments to all tests
#
#==========================================================

# Link to test_main
	
foreach(tname ${DOLPHIN_ALL_TESTS})
	target_link_libraries(${tname} test_main)
endforeach(tname)	


# Add tests

foreach(tname ${DOLPHIN_ALL_TESTS})
	add_test(NAME ${tname} COMMAND ${tname})
endforeach(tname)



    